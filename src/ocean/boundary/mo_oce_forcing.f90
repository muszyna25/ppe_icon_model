!>
!! Provide an implementation of the ocean forcing.
!!
!! Provide an implementation of the parameters used for surface forcing
!! of the hydrostatic ocean model.
!!
!! @author Peter Korn, MPI
!! @author Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
!!  Modification by Stephan Lorenz, MPI-M (2010-06):
!!   - renaming and adjustment to ocean domain and patch_oce
!!  Modified by Stephan Lorenz,     MPI-M (2010-07)
!!    adapted to structures discussed in 2010-01.
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_oce_forcing
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_ocean_nml,           ONLY: itestcase_oce, iforc_oce, analyt_forc, &
    & wstress_coeff, iforc_stat_oce, basin_height_deg, no_tracer, &
    & forcing_windstress_zonal_waveno, forcing_windstress_meridional_waveno, analytic_wind_amplitude
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_exception,           ONLY: finish, message
  USE mo_math_constants,      ONLY: pi, deg2rad
  USE mo_impl_constants,      ONLY: max_char_length, sea_boundary, success
  USE mo_math_utilities,      ONLY: gvec2cvec
  USE mo_sea_ice_types,       ONLY: t_sfc_flx
  USE mo_oce_state,           ONLY: set_oce_tracer_info

  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_var_list,            ONLY: add_var, add_ref, groups
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants
  
  IMPLICIT NONE
  PRIVATE
  
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  CHARACTER(LEN=12)           :: str_module    = 'oceForcing  '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  ! Public interface
  PUBLIC :: construct_ocean_forcing, destruct_ocean_forcing
  PUBLIC :: init_ocean_forcing
  
CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !! Constructor of surface fluxes for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  SUBROUTINE construct_ocean_forcing(p_patch, p_sfc_flx, var_list)
    !
    TYPE(t_patch),   INTENT(IN)    :: p_patch
    TYPE(t_sfc_flx), INTENT(INOUT) :: p_sfc_flx
    TYPE(t_var_list),INTENT(INOUT) :: var_list

    ! Local variables
    INTEGER                        :: alloc_cell_blocks, ist, jtrc, i
    CHARACTER(len=max_char_length) :: oce_tracer_names(no_tracer),&
    &                                 oce_tracer_units(no_tracer),&
    &                                 oce_tracer_longnames(no_tracer)
    INTEGER                        :: oce_tracer_codes(no_tracer)

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_ocean_forcing'

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    alloc_cell_blocks = p_patch%alloc_cell_blocks

    CALL add_var(var_list, 'forc_wind_u', p_sfc_flx%forc_wind_u , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_wind_u', 'm/s', 'forc_wind_u', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_wind_v', p_sfc_flx%forc_wind_v , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_wind_v', 'm/s', 'forc_wind_v', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_hflx', p_sfc_flx%forc_hflx , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_hflx', 'm/s', 'forc_hflx', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_fw_tot', p_sfc_flx%forc_fw_tot , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_fw_tot', 'm/s', 'forc_fw_tot', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_swflx', p_sfc_flx%forc_swflx , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_swflx', 'm/s', 'forc_swflx', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_lwflx', p_sfc_flx%forc_lwflx , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_lwflx', 'm/s', 'forc_lwflx', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_ssflx', p_sfc_flx%forc_ssflx , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_ssflx', 'm/s', 'forc_ssflx', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_slflx', p_sfc_flx%forc_slflx , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_slflx', 'm/s', 'forc_slflx', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_precip', p_sfc_flx%forc_precip , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_precip', 'm/s', 'forc_precip', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_evap', p_sfc_flx%forc_evap , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_evap', 'm/s', 'forc_evap', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_runoff', p_sfc_flx%forc_runoff , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_runoff', 'm/s', 'forc_runoff', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_fw_bc', p_sfc_flx%forc_fw_bc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_fw_bc', 'm/s', 'forc_fw_bc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_fw_bc_oce', p_sfc_flx%forc_fw_bc_oce , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_fw_bc_oce', 'm/s', 'forc_fw_bc_oce', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_fw_bc_ice', p_sfc_flx%forc_fw_bc_ice , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_fw_bc_ice', 'm/s', 'forc_fw_bc_ice', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_hfrelax', p_sfc_flx%forc_hfrelax , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_hfrelax', 'm/s', 'forc_hfrelax', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_fwrelax', p_sfc_flx%forc_fwrelax , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_fwrelax', 'm/s', 'forc_fwrelax', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_fw_ice_vol', p_sfc_flx%forc_fw_ice_vol, &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_fw_ice_vol', 'm/s', 'forc_fw_ice_vol', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    IF(no_tracer>=1) THEN
      ! there are four tracer related fields: tracer focing, tracer relaxation
      ! and both accumulated
      CALL add_var(var_list, 'forc_tracer', p_sfc_flx%forc_tracer , &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('forc_tracer', 'm/s', 'forc_tracer', DATATYPE_FLT32),&
        &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks,no_tracer/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      CALL add_var(var_list, 'forc_tracer_relax', p_sfc_flx%forc_tracer_relax , &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('forc_tracer_relax', 'm/s', 'forc_tracer_relax', DATATYPE_FLT32),&
        &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks,no_tracer/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      CALL add_var(var_list, 'forc_tracer_acc', p_sfc_flx%forc_tracer_acc , &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('forc_tracer_acc', 'm/s', 'forc_tracer_acc', DATATYPE_FLT32),&
        &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks,no_tracer/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      CALL add_var(var_list, 'forc_tracer_relax_acc', p_sfc_flx%forc_tracer_relax_acc , &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('forc_tracer_relax_acc', 'm/s', 'forc_tracer_relax_acc', DATATYPE_FLT32),&
        &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks,no_tracer/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ALLOCATE(p_sfc_flx%tracer_ptr(no_tracer*4))

      CALL set_oce_tracer_info(no_tracer           , &
        &                      oce_tracer_names    , &
        &                      oce_tracer_longnames, &
        &                      oce_tracer_codes    , &
        &                      oce_tracer_units)
      DO jtrc = 1,no_tracer
        CALL add_ref( var_list, 'forc_tracer', &
          &           'forc_tracer_'//TRIM(oce_tracer_names(jtrc)),          &
          &           p_sfc_flx%tracer_ptr(jtrc)%p,    &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,&
          &           t_cf_var('forc_tracer'//TRIM(oce_tracer_names(jtrc)), &
          &                    oce_tracer_units(jtrc), &
          &                    'forcing: '//TRIM(oce_tracer_longnames(jtrc)), DATATYPE_FLT32), &
          &           t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
          &           ldims=(/nproma,alloc_cell_blocks/))
      END DO
      DO jtrc = no_tracer+1,2*no_tracer
        i = jtrc - no_tracer
        CALL add_ref( var_list, 'forc_tracer_relax', &
          &           'forc_tracer_relax_'//TRIM(oce_tracer_names(i)),          &
          &           p_sfc_flx%tracer_ptr(jtrc)%p,    &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,&
          &           t_cf_var('forc_tracer_relax'//TRIM(oce_tracer_names(i)), &
          &                    oce_tracer_units(i), &
          &                    'forcing relaxation accumulated: '//TRIM(oce_tracer_longnames(i)), DATATYPE_FLT32), &
          &           t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
          &           ldims=(/nproma,alloc_cell_blocks/))
      END DO
      DO jtrc = (2*no_tracer)+1,3*no_tracer
        i = jtrc - 2*no_tracer
        CALL add_ref( var_list, 'forc_tracer_acc', &
          &           'forc_tracer_acc_'//TRIM(oce_tracer_names(i)),          &
          &           p_sfc_flx%tracer_ptr(jtrc)%p,    &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,&
          &           t_cf_var('forc_tracer_acc'//TRIM(oce_tracer_names(i)), &
          &                    oce_tracer_units(i), &
          &                    'forcing accumulated: '//TRIM(oce_tracer_longnames(i)), DATATYPE_FLT32), &
          &           t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
          &           ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
      END DO
      DO jtrc = (3*no_tracer)+1,4*no_tracer
        i = jtrc - 3*no_tracer
        CALL add_ref( var_list, 'forc_tracer_relax_acc', &
          &           'forc_tracer_relax_acc_'//TRIM(oce_tracer_names(i)),          &
          &           p_sfc_flx%tracer_ptr(jtrc)%p,    &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,&
          &           t_cf_var('forc_tracer_relax_acc'//TRIM(oce_tracer_names(i)), &
          &                    oce_tracer_units(i), &
          &                    'forcing relaxation accumulated: '//TRIM(oce_tracer_longnames(i)), DATATYPE_FLT32), &
          &           t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
          &           ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
      END DO
    ENDIF
    CALL add_var(var_list, 'forc_wind_u_acc', p_sfc_flx%forc_wind_u_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_wind_u_acc', 'm/s', 'forc_wind_u_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'forc_wind_v_acc', p_sfc_flx%forc_wind_v_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_wind_v_acc', 'm/s', 'forc_wind_v_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'forc_hflx_acc', p_sfc_flx%forc_hflx_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_hflx_acc', 'm/s', 'forc_hflx_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'forc_fw_tot_acc', p_sfc_flx%forc_fw_tot_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_fw_tot_acc', 'm/s', 'forc_fw_tot_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'forc_swflx_acc', p_sfc_flx%forc_swflx_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_swflx_acc', 'm/s', 'forc_swflx_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'forc_lwflx_acc', p_sfc_flx%forc_lwflx_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_lwflx_acc', 'm/s', 'forc_lwflx_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'forc_ssflx_acc', p_sfc_flx%forc_ssflx_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_ssflx_acc', 'm/s', 'forc_ssflx_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'forc_slflx_acc', p_sfc_flx%forc_slflx_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_slflx_acc', 'm/s', 'forc_slflx_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'forc_precip_acc', p_sfc_flx%forc_precip_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_precip_acc', 'm/s', 'forc_precip_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'forc_evap_acc', p_sfc_flx%forc_evap_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_evap_acc', 'm/s', 'forc_evap_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'forc_runoff_acc', p_sfc_flx%forc_runoff_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_runoff_acc', 'm/s', 'forc_runoff_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'forc_fw_bc_acc', p_sfc_flx%forc_fw_bc_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_fw_bc_acc', 'm/s', 'forc_fw_bc_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'forc_fw_bc_oce_acc', p_sfc_flx%forc_fw_bc_oce_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_fw_bc_oce_acc', 'm/s', 'forc_fw_bc_oce_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'forc_fw_bc_ice_acc', p_sfc_flx%forc_fw_bc_ice_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_fw_bc_ice_acc', 'm/s', 'forc_fw_bc_ice_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'forc_hfrelax_acc', p_sfc_flx%forc_hfrelax_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_hfrelax_acc', 'm/s', 'forc_hfrelax_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'forc_fwrelax_acc', p_sfc_flx%forc_fwrelax_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_fwrelax_acc', 'm/s', 'forc_fwrelax_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'forc_fw_ice_vol_acc', p_sfc_flx%forc_fw_ice_vol_acc, &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_fw_ice_vol_acc', 'm/s', 'forc_fw_ice_vol_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    ! cartesians
    ALLOCATE(p_sfc_flx%forc_wind_cc(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for forcing wind_cc  failed')
    END IF

    ! init of cartesian coordinates:
    p_sfc_flx%forc_wind_cc(:,:)%x(1) = 0.0_wp
    p_sfc_flx%forc_wind_cc(:,:)%x(2) = 0.0_wp
    p_sfc_flx%forc_wind_cc(:,:)%x(3) = 0.0_wp

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE construct_ocean_forcing
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !>
  !! Destructor surface flux forcing for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  SUBROUTINE destruct_ocean_forcing(p_sfc_flx)
    TYPE(t_sfc_flx), INTENT(INOUT) :: p_sfc_flx
    !
    ! Local variables

    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_ocean_forcing'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    ! forcing fields are handled by the ocean_default_list

    DEALLOCATE(p_sfc_flx%forc_wind_cc, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for forcing wind cc failed')
    END IF
    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE destruct_ocean_forcing
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !
  !>
  !! Initialization of stationary surface fluxes for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2011-09)
  !
  SUBROUTINE init_ocean_forcing(p_patch_3d, p_sfc_flx)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: p_patch_3d
    TYPE(t_sfc_flx)                             :: p_sfc_flx
    
    ! Local variables
    INTEGER :: jc, jb
    INTEGER :: i_startidx_c, i_endidx_c
    
    REAL(wp) :: z_lat, z_lon, y_length, y_center
    
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_forcing:init_ho_sfcflx'
    
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d   => p_patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )
    
    all_cells => patch_2d%cells%ALL
    
    ! analytical forcing
    IF (iforc_oce == analyt_forc) THEN
      
      IF (itestcase_oce == 27 .OR. itestcase_oce == 29) iforc_stat_oce = 1
      
      SELECT CASE (iforc_stat_oce)
      
      CASE (0)
        CALL message(TRIM(routine), &
          & 'iforc_stat_oce=0: no stationary wind forcing applied' )
        
      CASE (1)
        CALL message(TRIM(routine), 'Testcase (27,29): Apply stationary wind forcing' )
        y_length                        = basin_height_deg * deg2rad
        forcing_windstress_zonal_waveno = 1.0_wp
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
          
          DO jc = i_startidx_c, i_endidx_c
            
            IF(p_patch_3d%lsm_c(jc,1,jb)<=sea_boundary)THEN
              
              z_lat = patch_2d%cells%center(jc,jb)%lat
              z_lon = patch_2d%cells%center(jc,jb)%lon
              
              p_sfc_flx%forc_wind_u(jc,jb) = wstress_coeff * COS(forcing_windstress_zonal_waveno*pi*(z_lat-y_length)/y_length)
            ELSE
              p_sfc_flx%forc_wind_u(jc,jb) = 0.0_wp
            ENDIF
            p_sfc_flx%forc_wind_v(jc,jb) = 0.0_wp
            
            !Init cartesian wind
            CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),&
              & p_sfc_flx%forc_wind_v(jc,jb),&
              & patch_2d%cells%center(jc,jb)%lon,&
              & patch_2d%cells%center(jc,jb)%lat,&
              & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
              & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
              & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
            
          END DO
        END DO
        
      CASE (2)
        CALL message(TRIM(routine), &
          & 'iforc_stat_oce=2: stationary wind forcing over basin - u=cos(n*(lat-lat_0)/lat_0)')
        
        ! Latitudes vary from -pi/2 to pi/2
        y_length                        = basin_height_deg * deg2rad
        forcing_windstress_zonal_waveno = 1.0_wp
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
          DO jc = i_startidx_c, i_endidx_c
            z_lat = patch_2d%cells%center(jc,jb)%lat
            z_lon = patch_2d%cells%center(jc,jb)%lon
            IF (p_patch_3d%lsm_c(jc,1,jb)<=sea_boundary) THEN
              p_sfc_flx%forc_wind_u(jc,jb) =  wstress_coeff * COS(forcing_windstress_zonal_waveno*pi*(z_lat-y_length)/y_length)
            ELSE
              p_sfc_flx%forc_wind_u(jc,jb) = 0.0_wp
            ENDIF
            p_sfc_flx%forc_wind_v(jc,jb) = 0.0_wp
            
            !Init cartesian wind
            IF(p_patch_3d%lsm_c(jc,1,jb)<=sea_boundary)THEN
              CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),      &
                & p_sfc_flx%forc_wind_v(jc,jb),      &
                & patch_2d%cells%center(jc,jb)%lon,   &
                & patch_2d%cells%center(jc,jb)%lat,   &
                & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
            ELSE
              p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
            ENDIF
          END DO
        END DO
        
      CASE (3)
        CALL message(TRIM(routine), &
          & 'iforc_stat_oce=3: apply stationary wind forcing globally - u=cos(n*lat/lat_0); v=0')
        
        ! Use here global scale:
        y_length                        = 180.0_wp * deg2rad
        y_center                        = -60.0_wp * deg2rad
        forcing_windstress_zonal_waveno = 3.0_wp
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
          DO jc = i_startidx_c, i_endidx_c
            z_lat = patch_2d%cells%center(jc,jb)%lat
            z_lon = patch_2d%cells%center(jc,jb)%lon
            IF (p_patch_3d%lsm_c(jc,1,jb)<=sea_boundary) THEN
              p_sfc_flx%forc_wind_u(jc,jb) = analytic_wind_amplitude * &
                & COS(forcing_windstress_zonal_waveno*pi*(z_lat-y_center)/y_length)
            ELSE
              p_sfc_flx%forc_wind_u(jc,jb) = 0.0_wp
            ENDIF
            p_sfc_flx%forc_wind_v(jc,jb) = 0.0_wp
            
            !Init cartesian wind
            IF(p_patch_3d%lsm_c(jc,1,jb)<=sea_boundary)THEN
              CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),      &
                & p_sfc_flx%forc_wind_v(jc,jb),      &
                & patch_2d%cells%center(jc,jb)%lon,   &
                & patch_2d%cells%center(jc,jb)%lat,   &
                & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
            ELSE
              p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
            ENDIF
          END DO
        END DO
        
      CASE (4)
        CALL message(TRIM(routine), &
          & 'iforc_stat_oce=4: stationary wind forcing: u=cos(n*lat)*cos(lat) for APE (still does not work)')
        
        ! Forcing for ape
        forcing_windstress_zonal_waveno      = 3.0_wp
        forcing_windstress_meridional_waveno = 3.0_wp
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
          DO jc = i_startidx_c, i_endidx_c
            z_lat = patch_2d%cells%center(jc,jb)%lat
            z_lon = patch_2d%cells%center(jc,jb)%lon
            IF (p_patch_3d%lsm_c(jc,1,jb)<=sea_boundary) THEN
              
              p_sfc_flx%forc_wind_u(jc,jb) =  wstress_coeff * analytic_wind_amplitude * &
                & COS(z_lat) * COS(forcing_windstress_zonal_waveno * z_lat) * COS(z_lon)
              
              p_sfc_flx%forc_wind_v(jc,jb) = - wstress_coeff * analytic_wind_amplitude * &
                & COS(z_lat) * COS(forcing_windstress_zonal_waveno * z_lat) * SIN(forcing_windstress_meridional_waveno * z_lon)
              
            ELSE
              p_sfc_flx%forc_wind_u(jc,jb) = 0.0_wp
              p_sfc_flx%forc_wind_v(jc,jb) = 0.0_wp
            ENDIF
            
            !Init cartesian wind
            IF(p_patch_3d%lsm_c(jc,1,jb)<=sea_boundary)THEN
              CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),      &
                & p_sfc_flx%forc_wind_v(jc,jb),      &
                & patch_2d%cells%center(jc,jb)%lon,   &
                & patch_2d%cells%center(jc,jb)%lat,   &
                & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
            ELSE
              p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
            ENDIF
          END DO
        END DO
        
      CASE (5)
        CALL message(TRIM(routine), &
          & 'iforc_stat_oce=5: stationary wind forcing: u=cos(n*lat)*cos(lat) for APE (still does not work)')
        
        ! Forcing for ape
        forcing_windstress_zonal_waveno      = 3.0_wp
        forcing_windstress_meridional_waveno = 3.0_wp
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
          DO jc = i_startidx_c, i_endidx_c
            z_lat = patch_2d%cells%center(jc,jb)%lat
            z_lon = patch_2d%cells%center(jc,jb)%lon
            IF (p_patch_3d%lsm_c(jc,1,jb)<=sea_boundary) THEN
              
              p_sfc_flx%forc_wind_u(jc,jb) =  wstress_coeff * analytic_wind_amplitude * &
                & COS(z_lat) * COS(forcing_windstress_zonal_waveno * z_lat)
              
              p_sfc_flx%forc_wind_v(jc,jb) = 0
            ELSE
              p_sfc_flx%forc_wind_u(jc,jb) = 0.0_wp
              p_sfc_flx%forc_wind_v(jc,jb) = 0.0_wp
            ENDIF
            
            !Init cartesian wind
            IF(p_patch_3d%lsm_c(jc,1,jb)<=sea_boundary)THEN
              CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),      &
                & p_sfc_flx%forc_wind_v(jc,jb),      &
                & patch_2d%cells%center(jc,jb)%lon,   &
                & patch_2d%cells%center(jc,jb)%lat,   &
                & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
            ELSE
              p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
            ENDIF
          END DO
        END DO
        
      CASE default
        
        CALL message(TRIM(routine), 'STOP: Stationary Analytical Forcing not implemented' )
        CALL finish(TRIM(routine), 'CHOSEN STATIONARY FORCING OPTION NOT SUPPORTED - TERMINATE')
        
      END SELECT
      
      !---------Debug Diagnostics-------------------------------------------
      idt_src=0  ! output print level - 0: print in any case
      CALL dbg_print('analytical forcing u'      ,p_sfc_flx%forc_wind_u, str_module,idt_src, in_subset=patch_2d%cells%owned)
      CALL dbg_print('analytical forcing v'      ,p_sfc_flx%forc_wind_v, str_module,idt_src, in_subset=patch_2d%cells%owned)
      !---------------------------------------------------------------------
      
    END IF
    
  END SUBROUTINE init_ocean_forcing
  
  
  
END MODULE mo_oce_forcing
