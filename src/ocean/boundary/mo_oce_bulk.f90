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
!!  Modification by Stephan Lorenz, MPI-M:
!!   - renaming and adjustment to ocean domain and patch_oce (2010-06)
!!   - for parallel ocean: 3-dim ocean grid in v_base        (2011-07)
!!   - adding OMIP fluxes for sea ice                        (2011-09)
!!   - restructuring code                                    (2014-03)
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
MODULE mo_oce_bulk
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2007
!
!-------------------------------------------------------------------------
!
USE mo_kind,                ONLY: wp
USE mo_parallel_config,     ONLY: nproma
USE mo_run_config,          ONLY: dtime, ltimer
USE mo_sync,                ONLY: sync_c, sync_patch_array, global_sum_array
USE mo_timer,               ONLY: timer_start, timer_stop, timer_coupling
USE mo_io_units,            ONLY: filename_max
USE mo_mpi,                 ONLY: my_process_is_stdio, p_io, p_bcast,                   &
  &                               p_comm_work_test, p_comm_work
USE mo_parallel_config,     ONLY: p_test_run
USE mo_netcdf_read,         ONLY: read_netcdf_data
USE mo_datetime,            ONLY: t_datetime
USE mo_time_config,         ONLY: time_config
USE mo_ext_data_types,      ONLY: t_external_data
USE mo_ocean_ext_data,      ONLY: ext_data
USE mo_grid_config,         ONLY: nroot
USE mo_ocean_nml,           ONLY: iforc_oce, forcing_timescale, relax_analytical_type,  &
  &                               no_tracer, n_zlev, basin_center_lat,                  &
  &                               basin_center_lon, basin_width_deg, basin_height_deg,  &
  &                               para_surfRelax_Temp, i_apply_surface_hflux,           &
  &                               para_surfRelax_Salt,  type_surfRelax_Temp, type_surfRelax_Salt,  &
  &                               No_Forcing, Analytical_Forcing, OMIP_FluxFromFile,    &
  &                               Atmo_FluxFromFile, Coupled_FluxFromAtmo, Coupled_FluxFromFile, &
  &                               i_sea_ice, forcing_enable_freshwater, &
  &                               forcing_set_runoff_to_zero,           &
  &                               forcing_windstress_u_type,            &
  &                               forcing_windstress_v_type,            &
  &                               forcing_fluxes_type,                  &
  &                               limit_elevation, seaice_limit, l_relaxsal_ice
USE mo_dynamics_config,     ONLY: nold
USE mo_model_domain,        ONLY: t_patch, t_patch_3D
USE mo_util_dbg_prnt,       ONLY: dbg_print
USE mo_dbg_nml,             ONLY: idbg_mxmn
USE mo_oce_types,           ONLY: t_hydro_ocean_state
USE mo_exception,           ONLY: finish, message, message_text
USE mo_math_constants,      ONLY: pi, deg2rad, rad2deg
USE mo_physical_constants,  ONLY: rho_ref, als, alv, tmelt, tf, mu, clw, albedoW_sim
USE mo_impl_constants,      ONLY: max_char_length, sea_boundary, MIN_DOLIC
USE mo_math_utilities,      ONLY: gvec2cvec, cvec2gvec
USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean
USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
USE mo_sea_ice,             ONLY: calc_bulk_flux_ice, calc_bulk_flux_oce,                  &
  &                               ice_slow, ice_fast

USE mo_ocean_coupling,      ONLY: couple_ocean_toatmo_fluxes

IMPLICIT NONE

! required for reading netcdf files
INCLUDE 'netcdf.inc'

PRIVATE

! Public interface
PUBLIC  :: update_surface_flux

! private routines
PRIVATE :: update_flux_analytical
PRIVATE :: update_flux_fromFile
PRIVATE :: update_flux_from_atm_flx
PRIVATE :: update_relaxation_flux
PRIVATE :: read_forc_data_oce
PRIVATE :: balance_elevation

CHARACTER(len=*), PARAMETER :: version = '$Id$'
CHARACTER(len=12)           :: str_module    = 'oceBulk     '  ! Output of module for 1 line debug
INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug

CONTAINS

  !-------------------------------------------------------------------------
  !
  !>
  !! Update surface flux forcing for hydrostatic ocean
  !!
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
!<OptimizeEssential>
  SUBROUTINE update_surface_flux(p_patch_3D, p_os, p_as, p_ice, Qatm, p_sfc_flx, jstep, datetime, p_op_coeff)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_atmos_for_ocean)                     :: p_as
    TYPE(t_atmos_fluxes)                        :: Qatm
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_sfc_flx)                             :: p_sfc_flx
    INTEGER, INTENT(IN)                         :: jstep
    TYPE(t_datetime), INTENT(INOUT)             :: datetime
    TYPE(t_operator_coeff),   INTENT(IN)        :: p_op_coeff
    !
    ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:update_surface_flux'
    INTEGER  :: jc, jb, trac_no
    INTEGER  :: i_startidx_c, i_endidx_c
    REAL(wp) :: dsec
    REAL(wp) ::  z_c2(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), POINTER     :: t_top(:,:), s_top(:,:)

    TYPE(t_patch), POINTER:: p_patch 
    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain

    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------

    all_cells       => p_patch%cells%all
    cells_in_domain => p_patch%cells%in_domain

    t_top =>p_os%p_prog(nold(1))%tracer(:,1,:,1)
    s_top =>p_os%p_prog(nold(1))%tracer(:,1,:,2)

    !  calculate time
    dsec  = datetime%daysec        ! real seconds since begin of day

    SELECT CASE (iforc_oce)

    CASE (No_Forcing)                !  10

      ! CALL message(TRIM(routine), 'No  forcing applied' )
      CONTINUE

    CASE (Analytical_Forcing)        !  11

      !  Driving the ocean with analytically calculated fluxes
      CALL update_flux_analytical(p_patch_3D, p_os, p_sfc_flx)

    CASE (OMIP_FluxFromFile)         !  12

      !  Driving the ocean with OMIP (no NCEP activated anymore):
      !   1) read OMIP data (read relaxation data, type_surfRelax_Temp=2)
      !   2) call bulk formula (with or without sea ice)
      !   3) call ice model (todo: move below after case construct to unify call of sea ice)
      CALL update_flux_fromFile(p_patch_3D, p_os, p_as, p_ice, Qatm, p_sfc_flx, jstep, datetime, p_op_coeff)

    CASE (Atmo_FluxFromFile)                                          !  13

      ! 1) Read field data from file
      ! 2) CALL calc_atm_fluxes_from_bulk (p_patch, p_as, p_os, p_ice, Qatm)
      ! 3) CALL update_flux_from_atm_flx(p_patch, p_as, p_os, p_ice, Qatm, p_sfc_flx)

      CALL message(TRIM(routine), 'STOP: Ocean Forcing option 13 not implemented yet' )
      CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION NOT SUPPORTED - TERMINATE')

    CASE (Coupled_FluxFromAtmo)                                       !  14

      !  Driving the ocean in a coupled mode:
      !  atmospheric fluxes drive the ocean; fluxes are calculated by atmospheric model
      !  use atmospheric fluxes directly, i.e. avoid call to "calc_atm_fluxes_from_bulk"
      !  and do a direct assignment of atmospheric state to surface fluxes.
      !  3) call ice model (todo: move below after case construct to unify call of sea ice)
      !
      CALL couple_ocean_toatmo_fluxes(p_patch_3D, p_os, p_as, p_ice, Qatm, p_sfc_flx, jstep, datetime)

      ! call of sea ice model
      IF (i_sea_ice >= 1) THEN

        Qatm%SWnetw (:,:)   = p_sfc_flx%HeatFlux_ShortWave(:,:)
        Qatm%LWnetw (:,:)   = p_sfc_flx%HeatFlux_LongWave (:,:)
        Qatm%sensw  (:,:)   = p_sfc_flx%HeatFlux_Sensible (:,:)
        Qatm%latw   (:,:)   = p_sfc_flx%HeatFlux_Latent   (:,:)
     
      ! Precipitation on ice is snow when we're below the freezing point
            WHERE ( ALL( p_ice%Tsurf(:,:,:) < 0._wp, 2 ) )
                Qatm%rpreci(:,:) = p_sfc_flx%FrshFlux_SnowFall(:,:)
                Qatm%rprecw(:,:) = p_sfc_flx%FrshFlux_Precipitation(:,:)
            ELSEWHERE
                Qatm%rpreci(:,:) = 0._wp
                Qatm%rprecw(:,:) = p_sfc_flx%FrshFlux_Precipitation(:,:) + p_sfc_flx%FrshFlux_SnowFall(:,:)
            ENDWHERE
   
        CALL ice_slow(p_patch_3D, p_os, p_as, p_ice, Qatm, p_sfc_flx, p_op_coeff)
        
        ! Sum of freshwater flux (for salt) F = P - E + Sn (no runoff yet included??)
        p_sfc_flx%FrshFlux_TotalOcean(:,:) = p_patch_3d%wet_c(:,1,:)*( 1.0_wp-p_ice%concSum(:,:) ) * &
          &                                  (p_sfc_flx%FrshFlux_Precipitation(:,:) +                &
          &                                   p_sfc_flx%FrshFlux_Evaporation(:,:) +                  &
          &                                   p_sfc_flx%FrshFlux_SnowFall(:,:))
        p_sfc_flx%FrshFlux_TotalSalt(:,:)     =  p_sfc_flx%FrshFlux_TotalOcean(:,:) + p_sfc_flx%FrshFlux_TotalIce(:,:)
        ! sum of flux from sea ice to the ocean is stored in p_sfc_flx%HeatFlux_Total
        !  done in mo_sea_ice:upper_ocean_TS

        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src=1  ! output print level (1-5, fix)
        CALL dbg_print('UpdSfc: hi after slow'     ,p_ice%hi  ,str_module, idt_src, in_subset=p_patch%cells%owned)
        idt_src=3  ! output print level (1-5, fix)
        CALL dbg_print('UpdSfc: Conc. after slow'  ,p_ice%conc,str_module, idt_src, in_subset=p_patch%cells%owned)
        !---------------------------------------------------------------------

      ELSE

        p_sfc_flx%topBoundCond_windStress_u(:,:) = Qatm%stress_xw(:,:)
        p_sfc_flx%topBoundCond_windStress_v(:,:) = Qatm%stress_yw(:,:)

        p_sfc_flx%FrshFlux_TotalOcean(:,:) = p_patch_3d%wet_c(:,1,:)*( 1.0_wp-p_ice%concSum(:,:) ) * &
          &                                 (p_sfc_flx%FrshFlux_Precipitation(:,:) +                 &
          &                                  p_sfc_flx%FrshFlux_Evaporation(:,:) +                   &
          &                                  p_sfc_flx%FrshFlux_SnowFall(:,:))
        p_sfc_flx%FrshFlux_TotalSalt(:,:)  = p_sfc_flx%FrshFlux_TotalOcean(:,:)

      ENDIF


    CASE (Coupled_FluxFromFile)                                       !  15
      !1) bulk formula to atmospheric state and proceed as above, the only distinction
      !   to OMIP is that atmospheric info is coming from model rather than file

      CALL message(TRIM(routine), 'STOP: Ocean Forcing option 15 not implemented yet' )
      CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION NOT SUPPORTED - TERMINATE')

    CASE DEFAULT

      CALL message(TRIM(routine), 'STOP: Ocean Forcing option not implemented' )
      CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION DOES NOT EXIST - TERMINATE')

    END SELECT

    !
    ! After final updating of flux forcing for ocean (from file, bulk formula, or coupling) sea ice module is called
    !   - not yet possible, trying to disentangle update fluxes before and after sea ice (surface flux module)
    !
  ! IF (i_sea_ice >= 1) THEN

  !   Qatm%SWnetw (:,:)   = p_sfc_flx%HeatFlux_ShortWave(:,:)
  !   Qatm%LWnetw (:,:)   = p_sfc_flx%HeatFlux_LongWave(:,:)

  !   CALL ice_slow(p_patch_3D, p_os, p_as, p_ice, Qatm, p_sfc_flx, p_op_coeff)

  !   ! sum of flux from sea ice to the ocean is stored in p_sfc_flx%HeatFlux_Total
  !   !  done in mo_sea_ice:upper_ocean_TS

  !   !---------DEBUG DIAGNOSTICS-------------------------------------------
  !   idt_src=1  ! output print level (1-5, fix)
  !   CALL dbg_print('UpdSfc: hi after slow'     ,p_ice%hi       ,str_module,idt_src, in_subset=p_patch%cells%owned)
  !   idt_src=3  ! output print level (1-5, fix)
  !   CALL dbg_print('UpdSfc: Conc. after slow'  ,p_ice%conc     ,str_module,idt_src, in_subset=p_patch%cells%owned)
  !   !---------------------------------------------------------------------

    !
    ! After final updating of zonal and merdional components (from file, bulk formula, or coupling)
    ! cartesian coordinates are calculated
    !
    IF (iforc_oce > NO_FORCING) THEN
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF(p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
            CALL gvec2cvec(  p_sfc_flx%topBoundCond_windStress_u(jc,jb),&
                           & p_sfc_flx%topBoundCond_windStress_v(jc,jb),&
                           & p_patch%cells%center(jc,jb)%lon,&
                           & p_patch%cells%center(jc,jb)%lat,&
                           & p_sfc_flx%topBoundCond_windStress_cc(jc,jb)%x(1),&
                           & p_sfc_flx%topBoundCond_windStress_cc(jc,jb)%x(2),&
                           & p_sfc_flx%topBoundCond_windStress_cc(jc,jb)%x(3))
          ELSE
            p_sfc_flx%topBoundCond_windStress_u(jc,jb)         = 0.0_wp
            p_sfc_flx%topBoundCond_windStress_v(jc,jb)         = 0.0_wp
            p_sfc_flx%topBoundCond_windStress_cc(jc,jb)%x      = 0.0_wp
          ENDIF
        END DO
      END DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=1  ! output print level (1-5, fix)
      CALL dbg_print('UpdSfc: windStr u'       ,p_sfc_flx%topBoundCond_windStress_u      ,str_module,idt_src, &
        &  in_subset=p_patch%cells%owned)
      idt_src=2  ! output print level (1-5, fix)
      CALL dbg_print('UpdSfc: windStr v'       ,p_sfc_flx%topBoundCond_windStress_v      ,str_module,idt_src, &
        &  in_subset=p_patch%cells%owned)
      CALL dbg_print('UpdSfc: windStr cc%x(1)' ,p_sfc_flx%topBoundCond_windStress_cc%x(1),str_module,idt_src, &
        &  in_subset=p_patch%cells%owned)
      CALL dbg_print('UpdSfc: windStr cc%x(2)' ,p_sfc_flx%topBoundCond_windStress_cc%x(2),str_module,idt_src, &
        &  in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

    END IF

    !-------------------------------------------------------------------------
    ! Set surface boundary conditions to zero
    !p_sfc_flx%topBoundCond_Temp_vdiff(:,:) = 0.0_wp  ! heat flux BC not yet checked
    If (no_tracer>1) p_sfc_flx%topBoundCond_Salt_vdiff(:,:) = 0.0_wp

    !-------------------------------------------------------------------------
    ! Apply temperature relaxation to surface boundary condition
    !  - 2011-12: this is alternative to forcing by fluxes, not in addition

    IF (type_surfRelax_Temp >= 1) THEN

      trac_no = 1   !  tracer no 1: temperature
      CALL update_relaxation_flux(p_patch_3D, p_as, p_os, p_ice, p_sfc_flx, trac_no)

    END IF

    !-------------------------------------------------------------------------
    ! Apply net surface heat flux to boundary condition
    !  - heat flux is applied alternatively to temperature relaxation for coupling
    !  - also done if sea ice model is used since HeatFlux_Total is set in mo_sea_ice
    !  - with OMIP-forcing and sea_ice=0 we need type_surfRelax_Temp=-1
    !    since there is no HeatFlux_Total over open water when using OMIP-forcing
    !  - i_apply_surface_hflux=1 provides net surface heat flux globally
    !
    IF (no_tracer > 0) THEN
      IF (type_surfRelax_Temp == -1 .OR. i_sea_ice >= 1 .OR. i_apply_surface_hflux == 1) THEN

        ! Heat flux boundary condition for diffusion
        !   D = d/dz(K_v*dT/dz)  where
        ! Boundary condition at surface (upper bound of D at center of first layer)
        !   is calculated from net surface heat flux Q_s [W/m2]
        !   which is calculated by the atmosphere (coupled) or read from flux file (see above)
        !   Q_s = Rho*Cp*Q_T  with density Rho and Cp specific heat capacity
        !   K_v*dT/dz(surf) = Q_T = Q_s/Rho/Cp  [K*m/s]
        ! discretized:
        !   top_bc_tracer = topBoundCond_Temp_vdiff = HeatFlux_Total / (rho_ref*clw)

        p_sfc_flx%topBoundCond_Temp_vdiff(:,:) = p_sfc_flx%HeatFlux_Total(:,:) / (rho_ref*clw)

        !---------DEBUG DIAGNOSTICS-------------------------------------------
        CALL dbg_print('UpdSfc: HeatFlxTotal[W/m2]',p_sfc_flx%HeatFlux_Total         ,str_module,1, &
          &  in_subset=p_patch%cells%owned)
        CALL dbg_print('UpdSfc: topBC_T_vd[K*m/s]', p_sfc_flx%topBoundCond_Temp_vdiff,str_module,3, &
          &  in_subset=p_patch%cells%owned)
        !---------------------------------------------------------------------

      END IF
    END IF

    !-------------------------------------------------------------------------
    ! Apply salinity relaxation to surface boundary condition

    IF (type_surfRelax_Salt >= 1 .AND. no_tracer >1) THEN

      trac_no = 2   !  tracer no 2: salinity
      CALL update_relaxation_flux(p_patch_3D, p_as, p_os, p_ice, p_sfc_flx, trac_no)

    ENDIF  !  type_surfRelax_Salt >=1  salinity relaxation

    !-------------------------------------------------------------------------
    ! Apply freshwater forcing to surface boundary condition, independent of salinity relaxation

    ! Freshwater forcing activated as boundary condition in vertical Diffusion D, see above
    ! Vertical diffusion term for salinity Q_S in tracer equation and freshwater forcing W_s is
    !   Q_S = K_v*dS/dz(surf) = -W_s*S(nold)  [psu*m/s]

    IF (forcing_enable_freshwater) THEN

      p_sfc_flx%topBoundCond_Salt_vdiff(:,:) = p_sfc_flx%topBoundCond_Salt_vdiff(:,:) &
        &                                    - p_sfc_flx%FrshFlux_TotalSalt(:,:)*s_top(:,:)*p_patch_3d%wet_c(:,1,:)

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('UpdSfc: topBC_S_vd[psu*m/s]', p_sfc_flx%topBoundCond_Salt_vdiff,str_module,3, &
        &  in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

    ENDIF

    ! Sum of freshwater volume flux F = P - E + R + F_relax in [m/s] (independent of l_forc_frehsw)
    !  - add implicit freshwater flux due to relaxation to volume forcing term
    IF (no_tracer >1) THEN
      !old formulation <r14213:
      !p_sfc_flx%FrshFlux_VolumeTotal(:,:) = (p_sfc_flx%FrshFlux_TotalSalt(:,:) + p_sfc_flx%forc_fwrelax(:,:))
      p_sfc_flx%FrshFlux_VolumeTotal(:,:) = p_sfc_flx%FrshFlux_Runoff(:,:) + &
        &                          p_sfc_flx%FrshFlux_VolumeIce(:,:) + &
        &                          p_sfc_flx%FrshFlux_TotalOcean(:,:)  + &
        &                          p_sfc_flx%forc_fwrelax(:,:)

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('UpdSfc:VolumeTotal[m/s]',p_sfc_flx%FrshFlux_VolumeTotal,str_module,1, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------
    END IF
    
    ! apply additional volume flux to surface elevation
    !  - add to h_old before explicit term
    !  - no change in salt concentration
    !  - volume flux is considered for forcing_enable_freshwater=true only
    !    i.e. for salinity relaxation only, no volume flux is applied
    IF (forcing_enable_freshwater) THEN
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          p_os%p_prog(nold(1))%h(jc,jb) = p_os%p_prog(nold(1))%h(jc,jb) + p_sfc_flx%FrshFlux_VolumeTotal(jc,jb)*dtime
        END DO
      END DO
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('UpdSfc: h-old+fwfVol ',p_os%p_prog(nold(1))%h  ,str_module, 1, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------
    END IF
    
    ! apply volume flux correction: 
    !  - sea level is balanced to zero over ocean surface
    !  - correction applied daily
    IF (limit_elevation .AND. dsec < dtime) THEN
      CALL balance_elevation(p_patch_3D, p_os%p_prog(nold(1))%h)
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('UpdSfc: h-old+BalElev',p_os%p_prog(nold(1))%h  ,str_module, 2, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------
    END IF

  END SUBROUTINE update_surface_flux

  !-------------------------------------------------------------------------
  !
  !>
  !! Update surface flux forcing from file
  !!
  !! Provides surface forcing fluxes for ocean model from file.
  !!  Reads OMIP/NCEP fluxes via netcdf for bulk formula
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010/2014)
  !
  SUBROUTINE update_flux_fromFile(p_patch_3D, p_os, p_as, p_ice, Qatm, p_sfc_flx, jstep, datetime, p_op_coeff)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_atmos_for_ocean)                     :: p_as
    TYPE(t_atmos_fluxes)                        :: Qatm
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_sfc_flx)                             :: p_sfc_flx
    INTEGER, INTENT(IN)                         :: jstep
    TYPE(t_datetime), INTENT(INOUT)             :: datetime
    TYPE(t_operator_coeff),   INTENT(IN)        :: p_op_coeff
    !
    ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:update_flux_fromFile'
    INTEGER  :: jmon, jdmon, jmon1, jmon2, ylen, yday
    INTEGER  :: iniyear, curyear, offset
    INTEGER  :: jc, jb, no_set
    INTEGER  :: i_startidx_c, i_endidx_c
    REAL(wp) :: rday1, rday2, dtm1, dsec, z_smax
    REAL(wp) ::  z_c2(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) ::   Tfw(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), POINTER     :: t_top(:,:), s_top(:,:)

    TYPE(t_patch), POINTER:: p_patch 
    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain

    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------

    all_cells       => p_patch%cells%all
    cells_in_domain => p_patch%cells%in_domain

    t_top =>p_os%p_prog(nold(1))%tracer(:,1,:,1)
    s_top =>p_os%p_prog(nold(1))%tracer(:,1,:,2)

    !  calculate day and month
    jmon  = datetime%month         ! integer current month
    jdmon = datetime%day           ! integer day in month
    yday  = datetime%yeaday        ! integer current day in year
    ylen  = datetime%yealen        ! integer days in year (365 or 366)
    dsec  = datetime%daysec        ! real seconds since begin of day
   !ytim  = datetime%yeatim        ! real time since begin of year

    !-------------------------------------------------------------------------
    ! Applying annual forcing read from file in mo_ext_data:
    !  - stepping daily in monthly data (preliminary solution)

    !jdmon = mod(jdays+1,30)-1     ! no of days in month

    ! To Do: use fraction of month for interpolation
    !frcmon= datetime%monfrc       ! fraction of month
    !rday1 = frcmon+0.5_wp
    !rday2 = 1.0_wp-rday1
    !IF (rday1 > 1.0_wp)  THEN
    !  rday2=rday1
    !  rday1=1.0_wp-rday1
    !END IF

    !njday = int(86400._wp/dtime)  ! no of timesteps per day

    ! iforc_type: read time varying OMIP or NCEP flux forcing from file:
         ! 1: read wind stress (records 1, 2) and temperature (record 3)
         ! 2: read full OMIP dataset for bulk formula in mo_oce_bulk (12 records)
         ! 3: as 1; read surface heat (record 4) and freshwater flux (record 5) add.
         ! 4: as 1; read 4 parts of heat flux, precip/evap flux additionally
         ! 5: read full NCEP datasets; read monthly mean data of consecutive years

    ! Read forcing file in chunks of one year length fixed
    !  - #slo# 2012-02-17: first quick solution for reading NCEP data
    !  - ext_data has rank n_dom due to grid refinement in the atmosphere but not in the ocean

    ! Check if file should be read:
    !   - for iforc_type=5 only - NCEP type forcing
    !   - read annual data at Jan, 1st: seconds of year are less than a timestep
    !   - or at begin of each run (must not be first of january)
    !IF (iforc_type == 5) THEN
    IF (forcing_windstress_u_type == 5 .AND. forcing_windstress_v_type == 5 .AND. forcing_fluxes_type == 5) THEN
      dtm1 = dtime - 1.0_wp

      IF ( (jmon == 1 .AND. jdmon == 1 .AND. dsec < dtm1) .OR. (jstep == 1) ) THEN

        ! use initial date to define correct set (year) of reading NCEP data
        !  - with offset=0 always the first year of NCEP data is used
        iniyear = time_config%ini_datetime%year
        !curyear = time_config%cur_datetime%year  ! not updated each timestep
        curyear = datetime%year
        offset = 0
        no_set = offset + curyear-iniyear + 1 

        idt_src=2  ! output print level (1-5, fix)
     !  IF (idbg_mxmn >= idt_src) THEN
     !    WRITE(message_text,'(a,i2,a,i2,a,e15.5))') 'Read NCEP data: month=', &
     !      &  jmon,' day=',jdmon,' seconds=',dsec
     !    CALL message(TRIM(routine), message_text) 
        WRITE(message_text,'(a,3i5)') 'Read NCEP data: init. year, current year, no. of set:', &
          &                            iniyear, curyear, no_set
        CALL message(TRIM(routine), message_text) 
     !  END IF

        CALL read_forc_data_oce(p_patch, ext_data, no_set)

      END IF

    END IF

    !
    ! use annual forcing-data:
    !
    IF (forcing_timescale == 1)  THEN

      jmon1=1
      jmon2=1
      rday1=0.5_wp
      rday2=0.5_wp

    !
    ! interpolate monthly forcing-data daily:
    !
    ELSE IF (forcing_timescale == 12)  THEN

      jmon1=jmon-1
      jmon2=jmon
      rday1=REAL(15-jdmon,wp)/30.0_wp
      rday2=REAL(15+jdmon,wp)/30.0_wp
      IF (jdmon > 15)  THEN
        jmon1=jmon
        jmon2=jmon+1
        rday1=REAL(45-jdmon,wp)/30.0_wp
        rday2=REAL(jdmon-15,wp)/30.0_wp
      END IF

      IF (jmon1 ==  0) jmon1=12
      IF (jmon1 == 13) jmon1=1
      IF (jmon2 ==  0) jmon2=12
      IF (jmon2 == 13) jmon2=1

    !
    ! apply daily forcing-data directly:
    !
    ELSE

      ! - now daily data sets are read in mo_ext_data
      ! - use rday1, rday2, jmon1 = jmon2 = yday for controling correct day in year
      ! - no interpolation applied, 
      jmon1 = yday
      jmon2 = jmon1
      rday1 = 1.0_wp
      rday2 = 0.0_wp

      ! Leap year: read Feb, 28 twice since only 365 data-sets are available
      IF (ylen == 366) then
        IF (yday>59) jmon1=yday-1
        jmon2=jmon1
      ENDIF

    END IF

    !
    ! OMIP data read in mo_ext_data into variable ext_data
    !
    ! IF (iforc_type >= 1)  THEN
    IF (forcing_windstress_u_type > 0 .AND. forcing_windstress_u_type < 101 ) THEN ! file based forcing

      ! provide OMIP fluxes for wind stress forcing
      ! 1:  wind_u(:,:)   !  'stress_x': zonal wind stress       [Pa]
      ! 2:  wind_v(:,:)   !  'stress_y': meridional wind stress  [Pa]

      ! ext_data has rank n_dom due to grid refinement in the atmosphere but not in the ocean
      p_sfc_flx%topBoundCond_windStress_u(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,1) + &
        &                                        rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,1)

     ! Wind stress boundary condition for vertical diffusion D:
     !   D = d/dz(K_v*du/dz)  where
     ! Boundary condition at surface (upper bound of D at center of first layer)
     !   derived from wind-stress boundary condition Tau (in Pascal Pa=N/m2) read from OMIP data (or elsewhere)
     !   K_v*du/dz(surf) = F_D = Tau/Rho [ m2/s2 ]
     ! discretized:
     !   top_bc_u_c = topBoundCond_windStress_u / rho_ref
     !
     ! This is equivalent to an additonal forcing term F_u in the velocity equation, i.e. outside
     ! the vertical diffusion, following MITGCM:
     !   F_u = F_D/dz = Tau / (Rho*dz)  [ m/s2 ]

     ! The devision by rho_ref is done in top_bound_cond_horz_veloc (z_scale)

    END IF

    IF (forcing_windstress_v_type > 0 .AND. forcing_windstress_v_type < 101 ) THEN
      p_sfc_flx%topBoundCond_windStress_v(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,2) + &
        &                                        rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,2)
    END IF

    !IF (iforc_type == 2 .OR. iforc_type == 5) THEN
    IF (forcing_fluxes_type > 0 .AND. forcing_fluxes_type < 101 ) THEN

      !-------------------------------------------------------------------------
      ! provide OMIP fluxes for sea ice (interface to ocean)
      ! 4:  tafo(:,:),   &  ! 2 m air temperature                              [C]
      ! 5:  ftdew(:,:),  &  ! 2 m dew-point temperature                        [K]
      ! 6:  fu10(:,:) ,  &  ! 10 m wind speed                                  [m/s]
      ! 7:  fclou(:,:),  &  ! Fractional cloud cover
      ! 8:  pao(:,:),    &  ! Surface atmospheric pressure                     [hPa]
      ! 9:  fswr(:,:),   &  ! Incoming surface solar radiation                 [W/m]
      ! 10:  precip(:,:), &  ! precipitation rate                              [m/s]
      ! 11:  evap  (:,:), &  ! evaporation   rate                              [m/s]
      ! 12:  runoff(:,:)     ! river runoff  rate                              [m/s]
      ! 13: u(:,:),      &  ! 10m zonal wind speed                             [m/s]
      ! 14: v(:,:),      &  ! 10m meridional wind speed                        [m/s]

      p_as%tafo(:,:)  = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,4) + &
        &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,4)
      !  - change units to deg C, subtract tmelt (0 deg C, 273.15)
      p_as%tafo(:,:)  = p_as%tafo(:,:) - tmelt
      p_as%ftdew(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,5) + &
        &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,5)
      p_as%fu10(:,:)  = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,6) + &
        &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,6)
      p_as%fclou(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,7) + &
        &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,7)
      p_as%pao(:,:)   = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,8) + &
        &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,8)
      !  don't - change units to mb/hPa
      !p_as%pao(:,:)   = p_as%pao(:,:) !* 0.01
      p_as%fswr(:,:)  = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,9) + &
        &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,9)
      p_as%u(:,:)     = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,13) + &
        &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,13)
      p_as%v(:,:)     = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,14) + &
        &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,14)

      ! provide precipitation, evaporation, runoff flux data for freshwater forcing of ocean 
      !  - not changed via bulk formula, stored in surface flux data
      !  - Attention: as in MPIOM evaporation is calculated from latent heat flux (which is depentent on current SST)
      !               therefore not applied here
      p_sfc_flx%FrshFlux_Precipitation(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,10) + &
        &                                     rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,10)
      !p_sfc_flx%FrshFlux_Evaporation  (:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,11) + &
      !  &                                     rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,11)
      IF (forcing_set_runoff_to_zero) THEN
        p_sfc_flx%FrshFlux_Runoff(:,:) = 0.0_wp
      ELSE
        p_sfc_flx%FrshFlux_Runoff(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,12) + &
          &                              rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,12)
      ENDIF

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,4)
      CALL dbg_print('UpdSfc: Ext data4-ta/mon1' ,z_c2        ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,4)
      CALL dbg_print('UpdSfc: Ext data4-ta/mon2' ,z_c2        ,str_module,idt_src, in_subset=p_patch%cells%owned)
      CALL dbg_print('UpdSfc: p_as%tafo'         ,p_as%tafo   ,str_module,idt_src, in_subset=p_patch%cells%owned)

      IF (forcing_enable_freshwater) THEN
        idt_src=3  ! output print level (1-5, fix)
        CALL dbg_print('UpdSfc: p_sfc_flx%FrshFlux_Precipitation',p_sfc_flx%FrshFlux_Precipitation   ,str_module,idt_src, &
          & in_subset=p_patch%cells%owned)
        CALL dbg_print('UpdSfc: p_sfc_flx%FrshFlux_Runoff'       ,p_sfc_flx%FrshFlux_Runoff          ,str_module,idt_src, &
          & in_subset=p_patch%cells%owned)
      ENDIF
      !---------------------------------------------------------------------

    END IF  !  iforc_type=2 or 5

    !IF (iforc_type == 3) THEN

      !-------------------------------------------------------------------------
      ! #slo# This is a first try for "simple flux coupling" - not used anymore >r14000
      ! Apply surface heat and freshwater fluxes (records 4 and 5)
      ! 4:  hflx(:,:)   !  net surface heat flux               [W/m2]
      ! 5:  fwbc(:,:)   !  net freshwater flux                 [m/s]

    !END IF
    
    ! this is used for "intermediate complexity flux forcing" - not used anymore >r14000
    !IF (iforc_type == 4) THEN

      !-------------------------------------------------------------------------
      ! Apply 4 parts of surface heat and 2 parts of freshwater fluxes (records 4 to 9)
      ! 4:  swflx(:,:)   !  surface short wave heat flux        [W/m2]
      ! 5:  lwflx(:,:)   !  surface long  wave heat flux        [W/m2]
      ! 6:  ssflx(:,:)   !  surface sensible   heat flux        [W/m2]
      ! 7:  slflx(:,:)   !  surface latent     heat flux        [W/m2]
      ! 8:  precip(:,:)  !  total precipitation flux            [m/s]
      ! 9:  evap(:,:)    !  evaporation flux                    [m/s]


    !ENDIF  ! i_forc_type == 4

    IF (type_surfRelax_Temp == 2)  THEN

      !-------------------------------------------------------------------------
      ! Apply temperature relaxation data (record 3) from stationary forcing
      !  - change units to deg C, subtract tmelt (0 deg C, 273.15)
      !  - this is not done for type_surfRelax_Temp=3, since init-data is in Celsius

       p_sfc_flx%data_surfRelax_Temp(:,:) = &
         &  rday1*(ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,3)-tmelt) + &
         &  rday2*(ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,3)-tmelt)

    END IF

    IF (type_surfRelax_Salt == 2 .AND. no_tracer >1) THEN

      !-------------------------------------------------------------------------
      ! Apply salinity relaxation data (record ??) from stationary forcing

    !  p_sfc_flx%data_surfRelax_Salt(:,:) = &
    !    &  rday1*(ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,x)-tmelt) + &
    !    &  rday2*(ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,x)-tmelt)
      CALL finish(TRIM(ROUTINE),' type_surfRelax_Salt=2 (reading from flux file) not yet implemented')

    END IF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    IF (idbg_mxmn >= idt_src) THEN
      WRITE(message_text,'(a,i6,2(a,i4),2(a,f12.8))') 'FLUX time interpolation: jt=',jstep, &
        &  ' mon1=',jmon1,' mon2=',jmon2,' day1=',rday1,' day2=',rday2
      CALL message (' ', message_text)
    END IF
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,1)
    CALL dbg_print('FlxFil: Ext data1-u/mon1'  ,z_c2 ,str_module,idt_src, in_subset=p_patch%cells%owned)
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,1)
    CALL dbg_print('FlxFil: Ext data1-u/mon2'  ,z_c2 ,str_module,idt_src, in_subset=p_patch%cells%owned)
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,2)
    CALL dbg_print('FlxFil: Ext data2-v/mon1'  ,z_c2 ,str_module,idt_src, in_subset=p_patch%cells%owned)
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,2)
    CALL dbg_print('FlxFil: Ext data2-v/mon2'  ,z_c2 ,str_module,idt_src, in_subset=p_patch%cells%owned)
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,3)
    CALL dbg_print('FlxFil: Ext data3-t/mon1'  ,z_c2 ,str_module,idt_src, in_subset=p_patch%cells%owned)
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,3)
    CALL dbg_print('FlxFil: Ext data3-t/mon2'  ,z_c2 ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    IF (i_sea_ice >= 1) THEN

      !IF (iforc_type == 2 .OR. iforc_type == 5) THEN
      IF (forcing_fluxes_type > 0 .AND. forcing_fluxes_type < 101 ) THEN

        ! bulk formula are calculated globally using specific OMIP or NCEP fluxes
        CALL calc_bulk_flux_oce(p_patch, p_as, p_os , Qatm, datetime)
        CALL calc_bulk_flux_ice(p_patch, p_as, p_ice, Qatm, datetime)

        ! evaporation results from latent heat flux, as provided by bulk formula using OMIP/NCEP fluxes
        IF (forcing_enable_freshwater) THEN
          ! under sea ice evaporation is neglected, Qatm%latw is flux in the absence of sea ice
          ! TODO: evaporation of ice and snow must be implemented
          p_sfc_flx%FrshFlux_Evaporation(:,:) = Qatm%latw(:,:) / (alv*rho_ref)
          p_sfc_flx%FrshFlux_TotalOcean(:,:) = p_patch_3d%wet_c(:,1,:)*( 1.0_wp-p_ice%concSum(:,:) ) * & 
            &                                  ( p_sfc_flx%FrshFlux_Precipitation(:,:) + p_sfc_flx%FrshFlux_Evaporation(:,:) )
          ! Precipitation on ice is snow when we're below the freezing point
          ! TODO: Use 10 m temperature, not Tsurf - Also, do this in calc_bulk_flux_oce and
          ! calc_bulk_flux_ice
          WHERE ( ALL( p_ice%Tsurf(:,:,:) < 0._wp, 2 ) )
            Qatm%rpreci(:,:) = p_sfc_flx%FrshFlux_Precipitation(:,:)
            Qatm%rprecw(:,:) = 0._wp
          ELSEWHERE
            Qatm%rpreci(:,:) = 0._wp
            Qatm%rprecw(:,:) = p_sfc_flx%FrshFlux_Precipitation(:,:)
          ENDWHERE
        ENDIF

        ! TODO:
        !  - specify evaporation over snow/ice/water differently - currently only over open water is considered

      ENDIF
      
      IF ( no_tracer >= 2 ) THEN
        Tfw(:,:) = -mu*s_top(:,:)
      ELSE
        Tfw = Tf
      ENDIF
      CALL dbg_print('FlxFil: i-alb (bef ifast)'  ,Qatm%albvisdir ,str_module,4, in_subset=p_patch%cells%owned)

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        CALL ice_fast(i_startidx_c, i_endidx_c, nproma, p_ice%kice, dtime, &
          &   p_ice% Tsurf(:,:,jb),   &
          &   p_ice% T1   (:,:,jb),   &
          &   p_ice% T2   (:,:,jb),   &
          &   p_ice% hi   (:,:,jb),   &
          &   p_ice% hs   (:,:,jb),   &
          &   p_ice% Qtop (:,:,jb),   &
          &   p_ice% Qbot (:,:,jb),   & 
          &   Qatm%SWnet  (:,:,jb),   &
          &   Qatm%lat(:,:,jb) + Qatm%sens(:,:,jb) + Qatm%LWnet(:,:,jb),   & 
          &   Qatm%dlatdT(:,:,jb) + Qatm%dsensdT(:,:,jb) + Qatm%dLWdT(:,:,jb),   & 
          &   Tfw         (:,  jb),   &
          &   Qatm%albvisdir(:,:,jb), &
          &   Qatm%albvisdif(:,:,jb), &
          &   Qatm%albnirdir(:,:,jb), &
          &   Qatm%albnirdif(:,:,jb), &
          &   doy=datetime%yeaday)
      ENDDO

      ! Ocean albedo model
      Qatm%albvisdirw = albedoW_sim
      Qatm%albvisdifw = albedoW_sim
      Qatm%albnirdirw = albedoW_sim
      Qatm%albnirdifw = albedoW_sim

      ! #slo# 2012-12:
      ! sum of flux from sea ice to the ocean is stored in p_sfc_flx%HeatFlux_Total
      ! diagnosis of 4 parts is stored in p_sfc_flx%HeatFlux_ShortWave/LongWave/Sensible/Latent
      ! this diagnosis is done in mo_sea_ice:upper_ocean_TS
      ! 
      ! under ice the conductive heat flux is not yet stored specifically
      ! the sum HeatFlux_Total is aggregated and stored accordingly which cannot be done here

      ! ATTENTION
      !   ice_slow sets the fluxes in Qatm to zero for a new accumulation in ice_fast
      !   this should be done by the coupler if ice_fast is moved to the atmosphere

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=4  ! output print level (1-5, fix)
      CALL dbg_print('FlxFil: hi before slow'    ,p_ice%hi       ,str_module,idt_src, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: Conc. before slow' ,p_ice%conc     ,str_module,idt_src, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: ConcSum. bef slow' ,p_ice%concSum  ,str_module,idt_src, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: T1 before slow'    ,p_ice%t1       ,str_module,idt_src, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: T2 before slow'    ,p_ice%t2       ,str_module,idt_src, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: TSurf before slow' ,p_ice%tsurf    ,str_module,idt_src, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: Qtop  before slow' ,p_ice%Qtop     ,str_module,idt_src, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: Qbot  before slow' ,p_ice%Qbot     ,str_module,idt_src, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: i-alb before slow' ,Qatm%albvisdir ,str_module,idt_src, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

      CALL ice_slow(p_patch_3D, p_os, p_as, p_ice, Qatm, p_sfc_flx, p_op_coeff)

      !IF ( forcing_enable_freshwater .AND. (iforc_type == 2 .OR. iforc_type == 5) ) THEN
      IF ( forcing_enable_freshwater .AND. (forcing_fluxes_type > 0 .AND. forcing_fluxes_type < 101 ) ) THEN

        p_sfc_flx%FrshFlux_TotalSalt(:,:) = p_sfc_flx%FrshFlux_Runoff(:,:)                        &
          &           + p_sfc_flx%FrshFlux_TotalIce(:,:) + p_sfc_flx%FrshFlux_TotalOcean(:,:)

        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src=2  ! output print level (1-5, fix)
        CALL dbg_print('FlxFil:OMIP/NCEP:TotalSalt',p_sfc_flx%FrshFlux_TotalSalt,  str_module,idt_src,in_subset=p_patch%cells%owned)
        idt_src=3  ! output print level (1-5, fix)
        CALL dbg_print('FlxFil:OMIP/NCEP:Evap  '   ,p_sfc_flx%FrshFlux_Evaporation,str_module,idt_src,in_subset=p_patch%cells%owned)
        CALL dbg_print('FlxFil:OMIP/NCEP:FlxTotIce',p_sfc_flx%FrshFlux_TotalIce,   str_module,idt_src,in_subset=p_patch%cells%owned)
        CALL dbg_print('FlxFil:OMIP/NCEP:FlxTotOce',p_sfc_flx%FrshFlux_TotalOcean, str_module,idt_src,in_subset=p_patch%cells%owned)
        !---------------------------------------------------------------------

      ENDIF

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=4  ! output print level (1-5, fix)
      CALL dbg_print('FlxFil: hi after slow'     ,p_ice%hi       ,str_module, 4, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: Conc. after slow'  ,p_ice%conc     ,str_module, 4, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: ConcSum after slow',p_ice%concSum  ,str_module, 4, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: T1 after slow'     ,p_ice%t1       ,str_module, 4, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: T2 after slow'     ,p_ice%t2       ,str_module, 4, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: TSurf before slow' ,p_ice%tsurf    ,str_module, 4, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

      ! limit sea ice thickness to seaice_limit of surface layer depth, without elevation
      !   - no energy balance correction
      !   - number of ice classes currently kice=1 - sum of classes must be limited
      !   - only sea ice, no snow is considered
      IF (seaice_limit < 0.999999_wp) THEN
        z_smax = seaice_limit*p_patch_3D%p_patch_1D(1)%del_zlev_m(1)
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
          DO jc = i_startidx_c, i_endidx_c
            p_ice%hi(jc,:,jb) = MIN(p_ice%hi(jc,:,jb), z_smax)
          END DO
        END DO
      END IF

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('FlxFil: hi aft. limiter'     ,p_ice%hi       ,str_module, 4, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: Conc. aft. limiter'  ,p_ice%conc     ,str_module, 4, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: ConcSum aft. limiter',p_ice%concSum  ,str_module, 4, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: T1 aft. limiter'     ,p_ice%t1       ,str_module, 4, in_subset=p_patch%cells%owned)
      CALL dbg_print('FlxFil: T2 aft. limiter'     ,p_ice%t2       ,str_module, 4, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

    ELSE   !  no sea ice

      ! bulk formula applied to boundary forcing for ocean model:
      !  - no sea ice and no temperature relaxation
      !  - apply net surface heat flux in W/m2
      IF (forcing_fluxes_type > 0 .AND. forcing_fluxes_type < 101 ) CALL calc_bulk_flux_oce(p_patch,p_as,p_os,Qatm,datetime)
      !IF (iforc_type == 2 .OR. iforc_type == 5) CALL calc_bulk_flux_oce(p_patch, p_as, p_os, Qatm, datetime)
      p_sfc_flx%topBoundCond_windStress_u(:,:) = Qatm%stress_xw(:,:)
      p_sfc_flx%topBoundCond_windStress_v(:,:) = Qatm%stress_yw(:,:)

      type_surfRelax_Temp = 0   !  hack

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c

          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN
            p_sfc_flx%HeatFlux_ShortWave(jc,jb) = Qatm%SWnetw(jc,jb) ! net SW radiation flux over water
            p_sfc_flx%HeatFlux_LongWave (jc,jb) = Qatm%LWnetw(jc,jb) ! net LW radiation flux over water
            p_sfc_flx%HeatFlux_Sensible (jc,jb) = Qatm%sensw (jc,jb) ! Sensible heat flux over water
            p_sfc_flx%HeatFlux_Latent   (jc,jb) = Qatm%latw  (jc,jb) ! Latent heat flux over water
          ELSE
            p_sfc_flx%HeatFlux_ShortWave(jc,jb) = 0.0_wp
            p_sfc_flx%HeatFlux_LongWave (jc,jb) = 0.0_wp
            p_sfc_flx%HeatFlux_Sensible (jc,jb) = 0.0_wp
            p_sfc_flx%HeatFlux_Latent   (jc,jb) = 0.0_wp
          END IF

        ENDDO
      ENDDO

      ! for the setup with bulk and without sea ice the threshold for temperature is set to tf
      WHERE (t_top(:,:) .LT. Tf)
        t_top(:,:) = Tf
      ENDWHERE

      ! sum of fluxes for ocean boundary condition
      p_sfc_flx%HeatFlux_Total(:,:) = p_sfc_flx%HeatFlux_ShortWave(:,:) + p_sfc_flx%HeatFlux_LongWave(:,:) &
        &                      + p_sfc_flx%HeatFlux_Sensible(:,:) + p_sfc_flx%HeatFlux_Latent(:,:)

    ENDIF  !  sea ice

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('FlxFil: Bulk SW-flux'      ,p_sfc_flx%HeatFlux_ShortWave,str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('FlxFil: Bulk LW-flux'      ,p_sfc_flx%HeatFlux_LongWave ,str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('FlxFil: Bulk Sens.  HF'    ,p_sfc_flx%HeatFlux_Sensible ,str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('FlxFil: Bulk Latent HF'    ,p_sfc_flx%HeatFlux_Latent   ,str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('FlxFil: Bulk Total  HF'    ,p_sfc_flx%HeatFlux_Total    ,str_module, 2, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE update_flux_fromFile

  !-------------------------------------------------------------------------
  !
  !> Calculates the temperature (and salinity) relaxation term for vertical diffusion boundary condition
  !!  
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2014)
  !
!<OptimizeEssential>
  SUBROUTINE update_relaxation_flux(p_patch_3D, p_as, p_os, p_ice, p_sfc_flx, tracer_no)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_atmos_for_ocean),      INTENT(IN)    :: p_as
    TYPE(t_hydro_ocean_state),    INTENT(IN)    :: p_os
    TYPE (t_sea_ice),             INTENT (IN)   :: p_ice
    TYPE(t_sfc_flx)                             :: p_sfc_flx
    INTEGER,                       INTENT(IN)   :: tracer_no       !  no of tracer: 1=temperature, 2=salinity

    !Local variables 
    INTEGER :: jc, jb
    INTEGER :: i_startidx_c, i_endidx_c
    REAL(wp) :: z_tmin, z_relax, z_topBCSalt_old
    REAL(wp) :: z_c        (nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), PARAMETER   :: seconds_per_month = 2.592e6_wp     ! TODO: use real month length
    TYPE(t_patch), POINTER :: p_patch
    REAL(wp),      POINTER :: t_top(:,:), s_top(:,:)
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------  
    p_patch         => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------

    all_cells => p_patch%cells%all

    t_top =>p_os%p_prog(nold(1))%tracer(:,1,:,1)
    s_top =>p_os%p_prog(nold(1))%tracer(:,1,:,2)


    IF (tracer_no == 1) THEN  ! temperature relaxation

      !  - set minimum temperature to tf (-1.9 deg C) for simple temp-relax
      !  - set to zero on land points
    
      !z_tmin = -1.0_wp
      z_tmin = tf  !  -1.9 deg C
    
      ! LL: this is not the proper check in this point, should be removed - #SLO#: after restructuring
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN
            p_sfc_flx%data_surfRelax_Temp(jc,jb) &
              & = max(p_sfc_flx%data_surfRelax_Temp(jc,jb), z_tmin)
          ELSE
            p_sfc_flx%data_surfRelax_Temp(jc,jb) = 0.0_wp
          END IF
        END DO
      END DO
    
      ! Temperature relaxation activated as boundary condition in vertical Diffusion D:
      !   D = d/dz(K_v*dT/dz)  where
      ! Boundary condition at surface (upper bound of D at center of first layer)
      !   is relaxation to temperature (tau = relaxation constant [1/s] ):
      !   K_v*dT/dz(surf) = Q_T = -dz/tau*(T-T*) [ K*m/s ]
      ! discretized: temperature-relaxation-data T* = T_data = data_surfRelax_Temp
      !   top_bc_tracer = topBoundCond_Temp_vdiff = -(del_zlev_m+h) / relax_param[s] * (tracer - data_surfRelax_Temp)
      !
      ! This is equivalent to an additonal forcing term in the tracer equation, i.e. outside
      ! the vertical diffusion, following MITGCM:
      !    F_T  = Q_T/dz = -1/tau * (T-T*) [ K/s ]
      ! when using the sign convention
      !   dT/dt = Operators + F_T
      ! i.e. F_T <0 for  T-T* >0 (i.e. decreasing temperature if it is warmer than relaxation data) 
      ! 
      ! Extended boundary condition (relaxation term plus heat flux) is not yet implemented
    
      ! EFFECTIVE RESTORING PARAMETER: 1.0_wp/(para_surfRelax_Temp*seconds_per_month)
    
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
    
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN
            z_relax = (p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) + p_os%p_prog(nold(1))%h(jc,jb)) / &
              &       (para_surfRelax_Temp*seconds_per_month)
            p_sfc_flx%topBoundCond_Temp_vdiff(jc,jb) = -z_relax*(t_top(jc,jb)-p_sfc_flx%data_surfRelax_Temp(jc,jb))
          ELSE
            p_sfc_flx%topBoundCond_Temp_vdiff(jc,jb) = 0.0_wp
          ENDIF
    
        END DO
      END DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=2  ! output print level (1-5, fix)
      z_c(:,:) = p_sfc_flx%data_surfRelax_Temp(:,:)
      CALL dbg_print('UpdRlx: T-relax: T*'       ,z_c, str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:) = p_sfc_flx%data_surfRelax_Temp(:,:)-t_top(:,:)
      CALL dbg_print('UpdRlx: T-relax: T*-T'     ,z_c, str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:) = p_sfc_flx%topBoundCond_Temp_vdiff(:,:)
      CALL dbg_print('UpdRlx: T-relax: T [K*m/s]',z_c, str_module,idt_src, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

      ! Heat flux diagnosed for all ocean only relaxation cases
      ! TODO: discriminate hflx and hfrelax

      ! Heat flux diagnosed for relaxation cases, see above
      !   Q_s = Rho*Cp*Q_T  [W/m2]  with density Rho and Cp specific heat capacity
      ! where
      !   Q_T = K_v*dT/dz(surf) = Q_s/Rho/Cp  [K*m/s]

      p_sfc_flx%HeatFlux_Total(:,:) = p_sfc_flx%topBoundCond_Temp_vdiff(:,:) * rho_ref * clw

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=1  ! output print level (1-5, fix)
      CALL dbg_print('UpdSfc:T-relax-hflx [W/m2]',p_sfc_flx%HeatFlux_Total,str_module,idt_src, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

    ELSE IF (tracer_no == 2) THEN  ! salinity relaxation

      ! Salinity relaxation activated as boundary condition in vertical Diffusion D:
      !   D = d/dz(K_v*dS/dz)  where
      ! Boundary condition at surface (upper bound of D at center of first layer)
      !   is relaxation to salinity (tau = relaxation constant [1/s] ):
      !   K_v*dS/dz(surf) = Q_S = -dz/tau*(S-S*) [ psu*m/s ]
      ! discretized: salinity-relaxation-data S* = S_data = data_surfRelax_Salt
      !   top_bc_tracer = topBoundCond_Salt_vdiff = -(del_zlev_m+h) / relax_param[s] * (tracer - data_surfRelax_Salt)
      !
      ! This is equivalent to an additonal forcing term in the tracer equation, i.e. outside
      ! the vertical diffusion, following MITGCM:
      !    F_S  = Q_S/dz = -1/tau * (S-S*) [ psu/s ]
      ! when using the sign convention
      !   dS/dt = Operators + F_S
      ! i.e. F_S <0 for  S-S* >0 (i.e. decreasing salinity if it is saltier than relaxation data) 
      ! note that the freshwater flux is opposite in sign to F_S, see below,
      ! i.e. fwf >0 for  S-S* >0 (i.e. increasing freshwater flux to decrease the salinity)

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN

            !z_relax = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)&
            !          &/(para_surfRelax_Temp*seconds_per_month)
            z_relax = (p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)+p_os%p_prog(nold(1))%h(jc,jb)) / &
              &       (para_surfRelax_Salt*seconds_per_month)
            ! 
            ! If sea ice is present (and l_relaxsal_ice), salinity relaxation is proportional to open water,
            !   under sea ice, no relaxation is applied, according to the procedure in MPIOM
            IF (l_relaxsal_ice .AND. i_sea_ice >=1) z_relax = (1.0_wp-p_ice%concsum(jc,jb))*z_relax

            z_topBCSalt_old              = p_sfc_flx%topBoundCond_Salt_vdiff(jc,jb)
            p_sfc_flx%topBoundCond_Salt_vdiff(jc,jb) = p_sfc_flx%topBoundCond_Salt_vdiff(jc,jb) &
              &                              -z_relax*(s_top(jc,jb)-p_sfc_flx%data_surfRelax_Salt(jc,jb))

            ! Diagnosed freshwater flux due to relaxation [m/s]
            ! this flux is applied as volume forcing in surface equation in fill_rhs4surface_eq_ab
            p_sfc_flx%forc_fwrelax(jc,jb) = (z_topBCSalt_old-p_sfc_flx%topBoundCond_Salt_vdiff(jc,jb)) / s_top(jc,jb)

          ELSE
            p_sfc_flx%topBoundCond_Salt_vdiff(jc,jb) = 0.0_wp
            p_sfc_flx%forc_fwrelax(jc,jb)  = 0.0_wp
          ENDIF
        END DO
      END DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=2  ! output print level (1-5, fix)
      CALL dbg_print('UpdRlx:forc-fwrelax[m/s]'  ,p_sfc_flx%forc_fwrelax  ,str_module,idt_src, in_subset=p_patch%cells%owned)
      idt_src=3  ! output print level (1-5, fix)
      z_c(:,:) = p_sfc_flx%data_surfRelax_Salt(:,:)
      CALL dbg_print('UpdRlx: S-relax: S*'       ,z_c                     ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:) = p_sfc_flx%data_surfRelax_Salt(:,:)-s_top(:,:)
      CALL dbg_print('UpdRlx: S-relax: S*-S'     ,z_c                     ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:) = p_sfc_flx%topBoundCond_Salt_vdiff(:,:)
      CALL dbg_print('UpdRlx: S-relax:S[psu*m/s]',z_c                     ,str_module,idt_src, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

    END IF  ! tracer_no

  END SUBROUTINE update_relaxation_flux

  !-------------------------------------------------------------------------
  !
  !> Takes thermal calc_atm_fluxes_from_bulk to calculate surface fluxes for ocean forcing:
  !!  heat, freshwater and momentum.
  !!  not active or tested yet (2012/08)
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011). Originally written by D. Notz.
  !
  SUBROUTINE update_flux_from_atm_flx(p_patch_3D, p_as, p_os, p_ice, Qatm, p_sfc_flx)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_atmos_for_ocean),      INTENT(IN)    :: p_as
    TYPE(t_hydro_ocean_state),    INTENT(IN)    :: p_os
    TYPE (t_sea_ice),             INTENT (IN)   :: p_ice
    TYPE (t_atmos_fluxes),        INTENT (INOUT):: Qatm
    TYPE(t_sfc_flx)                             :: p_sfc_flx

    !Local variables 
    REAL(wp) :: z_rho_w = 1.22_wp  !near surface air density [kg/m^3] cf. Large/Yeager, sect 4.1, p.17
    REAL(wp) :: z_C_d0, z_C_d1, z_C_d
    REAL(wp) :: z_norm, z_v, z_relax

    INTEGER :: jc, jb, i
    INTEGER :: i_startidx_c, i_endidx_c
    REAL(wp):: z_evap        (nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp):: z_Q_freshwater(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:update_flux_from_atm_flx'
    TYPE(t_patch), POINTER :: p_patch
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------  
    p_patch         => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    all_cells => p_patch%cells%all

    !Relaxation parameter from namelist for salinity.
    z_relax = para_surfRelax_Temp/(30.0_wp*24.0_wp*3600.0_wp)

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        DO i = 1, p_ice%kice
          !surface heat forcing as sum of sensible, latent, longwave and shortwave heat fluxes
          IF (p_ice% hi(jc,jb,i) > 0._wp)THEN

            p_sfc_flx%HeatFlux_Total(jc,jb)              &
              & =  Qatm%sens(jc,jb,i) + Qatm%lat(jc,jb,i)& ! Sensible + latent heat flux at ice surface
              & +  Qatm%LWnet(jc,jb,i)                   & ! net LW radiation flux over ice surface
              & +  Qatm%bot(jc,jb,i)                       ! Ocean heat flux at ice bottom 
                                                           ! liquid/solid  precipitation rate
            !                                                are zero

            !This prepares freshwater flux calculation below; eq. (64) in Marsland et al.
            z_evap(jc,jb) = Qatm%lat(jc,jb,i)/(als*z_rho_w)

          ELSE

            p_sfc_flx%HeatFlux_Total(jc,jb)            &
            & =  Qatm%sensw(jc,jb) + Qatm%latw(jc,jb)  & ! Sensible + latent heat flux over water
            & +  Qatm%LWnetw(jc,jb)                    & ! net LW radiation flux over water
            & +  Qatm%SWnetw(jc,jb)                      ! net SW radiation flux ove water
                                                         ! liquid/solid  precipitation rate are zero

           !This prepares freshwater flux calculation below; eq. (64) in Marsland et al.
            z_evap(jc,jb) = Qatm%latw(jc,jb)/(alv*z_rho_w)
          ENDIF
        END DO

        !calculate surface freshwater flux       
        !following MPI-OM as described in Marsland et al, formula (63)-(65)

        !calculate evaporation from latent heat flux and latent heat of vaporisation
        !This is (63) in Marsland et al.
        !+River runoff +glacial meltwater
        z_Q_freshwater(jc,jb) = (Qatm%rpreci(jc,jb) + Qatm%rprecw(jc,jb)) -  z_evap(jc,jb)  

        !Now the freshwater flux calculation is finished; this is (65) in Marsland et al.
        !Relaxation of top layer salinity to observed salinity
        !
        !  Attention, check consistency in the model:
        !   - salinity relaxation is here in addition to the formulation at the end of update_surface_flux
        !   - also, according to (65) of Marsland, there is a bug below:
        !     multiplication with S1 (tracer(2)) is missing
        !   - has to be checked and merged with salinity boundary condition in update_surface_flux
        !
        p_sfc_flx%topBoundCond_Salt_vdiff(jc,jb) =                         &
          & (p_patch_3D%p_patch_1D(1)%del_zlev_m(1)+z_Q_freshwater(jc,jb)) &
          & /p_patch_3D%p_patch_1D(1)%del_zlev_m(1)                        &  !  * tracer(jc,1,jb,2)
          & +z_relax*(p_os%p_prog(nold(1))%tracer(jc,1,jb,2)-p_sfc_flx%data_surfRelax_Salt(jc,jb))


        !calculate wind stress    
        z_norm = sqrt(p_as%u(jc,jb)*p_as%u(jc,jb)+p_as%v(jc,jb)*p_as%v(jc,jb))

        !calculate drag coefficient for wind following 
        ! Kara, Rochford, Hurlburt, Air-Sea Flux Estimates And the 1997-1998 Enso Event
        ! Boundary-Layer Meteorology, 103, 439-458 (2002)
        !
        z_v = MAX(2.5_wp, MIN(p_as%fu10(jc,jb),32.5_wp))

        z_C_d0 = 1.0E-3_wp*(0.692_wp+0.071_wp*z_v-0.00070_wp*z_norm)
        z_C_d1 = 1.0E-3_wp*(0.083_wp-0.0054_wp*z_v-0.000093_wp*z_norm)
        z_C_d  = z_C_d0 + z_C_d1*(p_as%tafo(jc,jb)-p_os%p_prog(nold(1))%tracer(jc,1,jb,1))

        !write(*,*)'final wind stress coeff',z_C_d
        p_sfc_flx%topBoundCond_windStress_u(jc,jb) = z_rho_w*z_C_d*z_norm &
          &  *(p_as%u(jc,jb)- p_os%p_diag%u(jc,1,jb))

        p_sfc_flx%topBoundCond_windStress_v(jc,jb) = z_rho_w*z_C_d*z_norm &
          &  *(p_as%v(jc,jb) - p_os%p_diag%v(jc,1,jb))
   
      END DO
    END DO

    IF (type_surfRelax_Temp==1) THEN

       p_sfc_flx%topBoundCond_Temp_vdiff(:,:)=  z_relax * &
         &                                      (p_sfc_flx%data_surfRelax_Temp(:,:)-p_os%p_prog(nold(1))%tracer(:,1,:,1))

    ENDIF

  END SUBROUTINE update_flux_from_atm_flx
  !-------------------------------------------------------------------------
  !
  !>
  !! Update surface flux forcing for hydrostatic ocean
  !!
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
!<OptimizeEssential>
  SUBROUTINE update_flux_analytical(p_patch_3D, p_os, p_sfc_flx)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)    :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_sfc_flx)                             :: p_sfc_flx
    !
    ! local variables
    INTEGER :: jc, jb
    INTEGER :: i_startidx_c, i_endidx_c
    !INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
    !INTEGER :: rl_start_c, rl_end_c

    REAL(wp) :: z_lat, z_lon, z_lat_deg
    REAL(wp) :: y_length               !basin extension in y direction in degrees
    REAL(wp) :: z_T_init(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: z_perlat, z_perlon, z_permax, z_perwid, z_relax, z_dst
    INTEGER  :: z_dolic
    REAL(wp) :: z_temp_max, z_temp_min, z_temp_incr
    !CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:update_flux_analytical'
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------  
    p_patch         => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------
    all_cells => p_patch%cells%all

    SELECT CASE (relax_analytical_type)

    CASE(30,32,27)

     IF(no_tracer>=1.AND.type_surfRelax_Temp/=0)THEN

        y_length = basin_height_deg * deg2rad
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
          DO jc = i_startidx_c, i_endidx_c

            IF(p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary)THEN

              z_T_init(jc,jb) = 20.0_wp- p_patch_3D%p_patch_1D(1)%zlev_m(1)*15.0_wp/4000.0_wp

              z_lat = p_patch%cells%center(jc,jb)%lat
              z_lon = p_patch%cells%center(jc,jb)%lon
 
              ! Add temperature perturbation at new values
              z_perlat = basin_center_lat + 0.1_wp*basin_height_deg
              z_perlon = basin_center_lon + 0.1_wp*basin_width_deg 
              z_permax  = 0.1_wp
              z_perwid  =  10.0_wp

              z_relax = para_surfRelax_Temp/(30.0_wp*24.0_wp*3600.0_wp)

             z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
             IF (z_dolic > MIN_DOLIC) THEN

               z_dst=sqrt((z_lat-z_perlat*deg2rad)**2+(z_lon-z_perlon*deg2rad)**2)

               IF(z_dst<=5.0_wp*deg2rad)THEN
                 z_T_init = z_T_init &
                 &        + z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2) &
                 &        * sin(pi*p_patch_3D%p_patch_1D(1)%zlev_m(1)/4000.0_wp)
               ENDIF
               ! up to here z_init is identically initialized than temperature

               !add local cold perturbation 
               IF(z_dst<=10.5_wp*deg2rad)THEN
                 z_T_init(jc,jb)= z_T_init(jc,jb) - exp(-(z_dst/(z_perwid*deg2rad))**2)
               ENDIF

               p_sfc_flx%data_surfRelax_Temp(jc,jb)=z_T_init(jc,jb)

               p_sfc_flx%topBoundCond_Temp_vdiff(jc,jb) =  z_relax * &          
                 &  ( p_sfc_flx%data_surfRelax_Temp(jc,jb)-p_os%p_prog(nold(1))%tracer(jc,1,jb,1) )

             END IF
           ELSE
             p_sfc_flx%topBoundCond_windStress_cc(jc,jb)%x(:) = 0.0_wp
             p_sfc_flx%topBoundCond_windStress_u(jc,jb)       = 0.0_wp
             p_sfc_flx%topBoundCond_windStress_v(jc,jb)       = 0.0_wp
           ENDIF 
        END DO
      END DO

    ENDIF

    CASE (33)
      IF(type_surfRelax_Temp>=1)THEN
        z_relax = para_surfRelax_Temp/(30.0_wp*24.0_wp*3600.0_wp)

        p_sfc_flx%topBoundCond_Temp_vdiff(:,:) = z_relax*( p_sfc_flx%data_surfRelax_Temp(:,:) &
          &                                               -p_os%p_prog(nold(1))%tracer(:,1,:,1) )

      END IF

    CASE(51)

      IF(type_surfRelax_Temp>=1)THEN

        z_relax = para_surfRelax_Temp/(30.0_wp*24.0_wp*3600.0_wp)

        z_temp_max  = 30.5_wp
        z_temp_min  = 0.5_wp
        z_temp_incr = (z_temp_max-z_temp_min)/(n_zlev-1.0_wp)

      !Add horizontal variation
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          z_lat = p_patch%cells%center(jc,jb)%lat
          z_lat_deg = z_lat*rad2deg

            IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN

              z_temp_max     =0.01_wp*(z_lat_deg-basin_center_lat)*(z_lat_deg-basin_center_lat)
              z_T_init(jc,jb)=30.5_wp

              z_T_init(jc,jb)&
              &=z_T_init(jc,jb)*exp(-z_temp_max/basin_height_deg)
            ELSE
              z_T_init(jc,jb)=0.0_wp
            ENDIF
        END DO
      END DO
      p_sfc_flx%data_surfRelax_Temp(:,:)=z_T_init(:,:)

      p_sfc_flx%topBoundCond_Temp_vdiff(:,:) = z_relax*( p_sfc_flx%data_surfRelax_Temp(:,:) &
        &                                               -p_os%p_prog(nold(1))%tracer(:,1,:,1) )

      END IF

    END SELECT

  END SUBROUTINE update_flux_analytical

  !-------------------------------------------------------------------------
  !>
  !! Balance sea level to zero over global ocean
  !!
  !! Balance sea level to zero over global ocean
  !! This routine uses parts of mo_oce_diagnostics
  !!
  !! @par Revision History
  !! Initial revision by Stephan Lorenz, MPI (2013-04)
  !!
  !!
  SUBROUTINE balance_elevation (p_patch_3D, h_old)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)    :: p_patch_3D
    REAL(wp), INTENT(INOUT)                 :: h_old(1:nproma,1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    TYPE(t_patch), POINTER                  :: p_patch
    TYPE(t_subset_range), POINTER           :: all_cells

    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: jc, jb
    REAL(wp) :: ocean_are, glob_slev, corr_slev

    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
 
    ! parallelize correctly
    ocean_are = p_patch_3D%p_patch_1D(1)%ocean_area(1)
    glob_slev = global_sum_array(p_patch%cells%area(:,:)*h_old(:,:)*p_patch_3D%wet_halo_zero_c(:,1,:))
    corr_slev = glob_slev/ocean_are

    idt_src=2
    IF ((my_process_is_stdio()) .AND. (idbg_mxmn >= idt_src)) &
      & write(0,*)' BALANCE_ELEVATION(Dom): ocean_are, glob_slev, corr_slev =',ocean_are, glob_slev, glob_slev/ocean_are

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc =  i_startidx_c, i_endidx_c
        IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN
          ! subtract or scale?
          h_old(jc,jb) = h_old(jc,jb) - corr_slev
          !h_old(jc,jb) = h_old(jc,jb) * (1.0_wp - corr_slev)
          !h_old(jc,jb) = h_old(jc,jb) - h_old(jc,jb)*corr_slev
        END IF
      END DO
    END DO

  END SUBROUTINE balance_elevation

  !-------------------------------------------------------------------------
  !>
  !! Read ocean forcing data from netcdf
  !!
  !! Read ocean forcing data for NCEP or other forcing
  !! This routine reads annual data sets of length forcing_timescale
  !!
  !! @par Revision History
  !! Initial revision by Stephan Lorenz, MPI (2012-02-17)
  !!
  !!
  SUBROUTINE read_forc_data_oce (p_patch, ext_data, no_set)

    TYPE(t_patch), INTENT(IN)            :: p_patch
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)
    INTEGER,       INTENT(IN)            :: no_set          !  no of set in file to be read

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_oce_bulk:read_forc_data_oce'

    CHARACTER(filename_max) :: ncep_file   !< file name for reading in

    LOGICAL :: l_exist
    INTEGER :: jg, i_lev, i_cell_type, no_cells, no_tst, jtime, jt !, jc, jb
    INTEGER :: ncid, dimid,mpi_comm
    INTEGER :: i_start(2),i_count(2), jcells

    REAL(wp):: z_flux(nproma,p_patch%alloc_cell_blocks,forcing_timescale)  ! set length is forcing_timescale, 3rd dimension
    REAL(wp):: z_c   (nproma,forcing_timescale,p_patch%alloc_cell_blocks)  ! 2nd dimension is forcing_timescale
    !TYPE (t_keyword_list), POINTER :: keywords => NULL()

    !-------------------------------------------------------------------------

    !  READ NCEP FORCING

    !-------------------------------------------------------------------------

    !CALL message (TRIM(routine), 'start')

    IF (iforc_oce == 12) THEN

    !DO jg = 1,n_dom
      jg = 1

      i_lev       = p_patch%level
      i_cell_type = p_patch%cell_type

      IF(my_process_is_stdio()) THEN
        !
        WRITE (ncep_file,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',i_lev, '-flux.nc'

        !ncep_file=TRIM('/pool/data/ICON/external/iconR2B04-flux.nc')
        CALL message( TRIM(routine),'Ocean NCEP forcing flux file is: '//TRIM(ncep_file) )
        INQUIRE (FILE=ncep_file, EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'NCEP forcing flux file is not found.')
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(ncep_file), NF_NOWRITE, ncid))
        !CALL message( TRIM(routine),'Ocean NCEP flux file opened for read' )

        !
        ! get and check number of cells in ncep data
        !
        CALL nf(nf_inq_dimid(ncid, 'ncells', dimid))
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))

        IF(p_patch%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in NCEP flux file do not match.')
        ENDIF

        !
        ! get number of timesteps
        !
        CALL nf(nf_inq_dimid(ncid, 'time', dimid))
        CALL nf(nf_inq_dimlen(ncid, dimid, no_tst))
        !
        ! check - s.b.

      ENDIF
      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF
      CALL p_bcast(no_tst, p_io, mpi_comm)

      !-------------------------------------------------------
      !
      ! Read 12 monthly NCEP data sets for triangle centers using 4-dim routine
      !
      !-------------------------------------------------------

      jcells = p_patch%n_patch_cells  !  global dimension
      jtime  = forcing_timescale              !  time period to read (not yet)


      ! provide NCEP fluxes for sea ice (interface to ocean)
      ! 1:  'stress_x': zonal wind stress       [Pa]
      ! 2:  'stress_y': meridional wind stress  [Pa]
      ! 3:  'SST"     : sea surface temperature [K]

      ! zonal wind stress
      !write(0,*) ' ncep set 1: dimensions:',p_patch%n_patch_cells_g, p_patch%n_patch_cells, &
      ! &  forcing_timescale, nproma, p_patch%nblks_c
      !CALL read_netcdf_data (ncid, 'stress_x', p_patch%n_patch_cells_g,      &
      !  &                    p_patch%n_patch_cells, p_patch%cells%decomp_info%glb_index, &
      !  &                    forcing_timescale, z_flx2(:,:,:))
      !write(0,*) ' READ_FORC, READ 1: first data sets: stress-x, block=5, index=1,5:'
      !do jt=1,jtime
      !  write(0,*) 'jt=',jt,' val:',(z_flx2(jc,jt,5),jc=1,5)
      !enddo

      ! start-pointer and length of pointer for reading data:
      ! start: first set (1,1); second year (1,jtime+1)
      i_start(1) = 1
      i_start(2) = jtime*(no_set-1) + 1  ! position pointer to set no_set
      i_count(1) = jcells                ! length of pointer, dim 1 of z_dummy_array
      i_count(2) = jtime                 ! length of pointer, dim 2 of z_dummy_array

      idt_src=2  ! output print level (1-5, fix)
      IF (idbg_mxmn >= idt_src) THEN
        !
        WRITE(message_text,'(A,I6,A)')  'Ocean NCEP flux file contains',no_tst,' data sets'
        CALL message( TRIM(routine), TRIM(message_text) )

        WRITE(message_text,'(4(A,I4))')  'NCEP data set: length = ',jtime, &
          &   '; no. of set =',no_set,                                     &
          &   '; pos. of ptr =', i_start(2)
        CALL message( TRIM(routine), TRIM(message_text) )
      END IF

      CALL read_netcdf_data (ncid, 'stress_x', p_patch%n_patch_cells_g,      &
        &                    p_patch%n_patch_cells, p_patch%cells%decomp_info%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))


      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,1) = z_flux(:,:,jt)
      END DO

      ! meridional wind stress
      CALL read_netcdf_data (ncid, 'stress_y', p_patch%n_patch_cells_g,      &
        &                    p_patch%n_patch_cells, p_patch%cells%decomp_info%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,2) = z_flux(:,:,jt)
      END DO

      ! SST
      CALL read_netcdf_data (ncid, 'SST', p_patch%n_patch_cells_g,           &
        &                    p_patch%n_patch_cells, p_patch%cells%decomp_info%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,3) = z_flux(:,:,jt)
      END DO

 !    ! Read complete NCEP data sets for focing ocean model (iforc_type=5)
 !    ! 4:  tafo(:,:),   &  ! 2 m air temperature                              [C]
 !    ! 5:  ftdew(:,:),  &  ! 2 m dew-point temperature                        [K]
 !    ! 6:  fu10(:,:) ,  &  ! 10 m wind speed                                  [m/s]
 !    ! 7:  fclou(:,:),  &  ! Fractional cloud cover
 !    ! 8:  pao(:,:),    &  ! Surface atmospheric pressure                     [hPa]
 !    ! 9:  fswr(:,:),   &  ! Incoming surface solar radiation                 [W/m]

      ! 2m-temperature
      CALL read_netcdf_data (ncid, 'temp_2m', p_patch%n_patch_cells_g,       &
        &                    p_patch%n_patch_cells, p_patch%cells%decomp_info%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,4) = z_flux(:,:,jt)
      END DO

      ! 2m dewpoint temperature
      CALL read_netcdf_data (ncid, 'dpt_temp_2m', p_patch%n_patch_cells_g,   &
        &                    p_patch%n_patch_cells, p_patch%cells%decomp_info%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,5) = z_flux(:,:,jt)
      END DO

      ! Scalar wind
      CALL read_netcdf_data (ncid, 'scalar_wind', p_patch%n_patch_cells_g,   &
        &                    p_patch%n_patch_cells, p_patch%cells%decomp_info%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,6) = z_flux(:,:,jt)
      END DO

      ! cloud cover
      CALL read_netcdf_data (ncid, 'cloud', p_patch%n_patch_cells_g,         &
        &                    p_patch%n_patch_cells, p_patch%cells%decomp_info%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,7) = z_flux(:,:,jt)
      END DO

      ! sea level pressure
      CALL read_netcdf_data (ncid, 'pressure', p_patch%n_patch_cells_g,      &
        &                    p_patch%n_patch_cells, p_patch%cells%decomp_info%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,8) = z_flux(:,:,jt)
      END DO

      ! total solar radiation
      CALL read_netcdf_data (ncid, 'tot_solar', p_patch%n_patch_cells_g,     &
        &                    p_patch%n_patch_cells, p_patch%cells%decomp_info%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,9) = z_flux(:,:,jt)
      END DO

      ! precipitation
  !   CALL read_netcdf_data (ncid, 'precip', p_patch%n_patch_cells_g,        &
  !     &                    p_patch%n_patch_cells, p_patch%cells%decomp_info%glb_index, &
  !     &                    jtime, i_start, i_count, z_flux(:,:,:))
  !   DO jt = 1, jtime
  !     ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,10) = z_flux(:,:,jt)
  !   END DO

      ! evaporation or downward surface LW flux
  !   CALL read_netcdf_data (ncid, 'evap', p_patch%n_patch_cells_g,          &
  !     &                    p_patch%n_patch_cells, p_patch%cells%decomp_info%glb_index, &
  !     &                    jtime, i_start, i_count, z_flux(:,:,:))
  !   DO jt = 1, jtime
  !     ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,11) = z_flux(:,:,jt)
  !   END DO
  !   CALL read_netcdf_data (ncid, 'dlwrf', p_patch%n_patch_cells_g,         &
  !     &                    p_patch%n_patch_cells, p_patch%cells%decomp_info%glb_index, &
  !     &                    jtime, i_start, i_count, z_flux(:,:,:))
  !   DO jt = 1, jtime
  !     ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,11) = z_flux(:,:,jt)
  !   END DO

      ! runoff
  !   CALL read_netcdf_data (ncid, 'runoff', p_patch%n_patch_cells_g,        &
  !     &                    p_patch%n_patch_cells, p_patch%cells%decomp_info%glb_index, &
  !     &                    jtime, i_start, i_count, z_flux(:,:,:))
  !   DO jt = 1, jtime
  !     ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,12) = z_flux(:,:,jt)
  !   END DO


      !
      ! close file
      !
      IF(my_process_is_stdio()) CALL nf(nf_close(ncid))

    !ENDDO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,1)
      CALL dbg_print('ReadFc: NCEP: stress-x'    ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,2)
      CALL dbg_print('ReadFc: NCEP: stress-y'    ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,3)
      CALL dbg_print('ReadFc: NCEP: SST'         ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      idt_src=4  ! output print level (1-5, fix)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,4)
      CALL dbg_print('ReadFc: NCEP: temp_2m'     ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,5)
      CALL dbg_print('ReadFc: NCEP: dpt_temp_2m' ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,6)
      CALL dbg_print('ReadFc: NCEP: scalar_wind' ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,7)
      CALL dbg_print('ReadFc: NCEP: cloudiness'  ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,8)
      CALL dbg_print('ReadFc: NCEP: pressure'    ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,9)
      CALL dbg_print('ReadFc: NCEP: total solar' ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
    ! z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,10)
    ! CALL dbg_print('ReadFc: NCEP: precip.'     ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
    ! z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,11)
    ! CALL dbg_print('ReadFc: NCEP: evaporation' ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
    ! z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,12)
    ! CALL dbg_print('ReadFc: NCEP: runoff'      ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

      idt_src=2  ! output print level (1-5, fix)
      IF (idbg_mxmn >= idt_src) &
        & CALL message( TRIM(routine),'Ocean NCEP fluxes for external data read' )

    END IF ! iforc_oce=12

  END SUBROUTINE read_forc_data_oce

  !-------------------------------------------------------------------------

!<OptimizeEssential>
  SUBROUTINE nf(status)

    INTEGER, INTENT(in) :: status

    IF (status /= nf_noerr) THEN
      CALL finish('mo_oce_bulk netCDF error', nf_strerror(status))
    ENDIF

  END SUBROUTINE nf


END MODULE mo_oce_bulk
