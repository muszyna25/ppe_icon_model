!>
!! Provide an implementation of the ocean surface module.
!!
!! Provide an implementation of the parameters used for surface forcing
!! of the hydrostatic ocean model.
!!
!! @author Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
!!  Restructured code by Stephan Lorenz, MPI-M: (2015-04)
!!  Restructured code by Vladimir Lapin, MPI-M: (2017-01)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ocean_bulk_forcing
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2007
!
!-------------------------------------------------------------------------
!
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_sync,                ONLY: global_sum_array
  USE mo_io_units,            ONLY: filename_max
  USE mo_mpi,                 ONLY: my_process_is_stdio, p_io, p_bcast, p_comm_work_test, p_comm_work
  USE mo_parallel_config,     ONLY: p_test_run
  USE mo_read_interface,      ONLY: openInputFile, closeFile, t_stream_id, &
    &                               on_cells, read_2D_time  !, read_3D
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_ocean_ext_data,      ONLY: ext_data
  USE mo_dynamics_config,     ONLY: nold
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_dbg_nml,             ONLY: idbg_mxmn

  USE mo_ocean_nml,           ONLY: iforc_oce, forcing_timescale,  forcing_frequency, &
    &  no_tracer, para_surfRelax_Temp, type_surfRelax_Temp,             &
    &  para_surfRelax_Salt, type_surfRelax_Salt,                                &
    &  i_sea_ice, l_relaxsal_ice, forcing_enable_freshwater,                    &
    &  forcing_set_runoff_to_zero, OMIP_FluxFromFile, OceanReferenceDensity
  USE mo_sea_ice_nml,         ONLY: use_calculated_ocean_stress, stress_ice_zero

  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_math_constants,      ONLY: rad2deg !, deg2rad
  USE mo_physical_constants,  ONLY: alv, tmelt, clw, stbo, zemiss_def
  USE mo_physical_constants,  ONLY: rd, cpd, fr_fac, alf, cd_ia, Cd_io, rho_ref
  USE mo_impl_constants,      ONLY: max_char_length, sea_boundary
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_ocean_surface_types, ONLY: t_ocean_surface, t_atmos_for_ocean

  USE mo_math_utilities,      ONLY: gvec2cvec
  USE mtime,                  ONLY: datetime, getDayOfYearFromDateTime, getNoOfDaysInYearDateTime

  
  IMPLICIT NONE
  
  ! required for reading netcdf files
  INCLUDE 'netcdf.inc'

  PRIVATE

  ! Public interface
  PUBLIC :: update_surface_relaxation
  PUBLIC :: apply_surface_relaxation

  PUBLIC :: update_flux_fromFile
  PUBLIC :: calc_omip_budgets_ice
  PUBLIC :: calc_omip_budgets_oce

  PUBLIC :: update_ocean_surface_stress

  PUBLIC :: balance_elevation


  CHARACTER(len=12)           :: str_module    = 'OceanBulkForcing'  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug
  REAL(wp), PARAMETER         :: seconds_per_month = 2.592e6_wp  ! TODO: use real month length

CONTAINS

!**********************************************************************
!----------------------------- Relaxation -----------------------------
!**********************************************************************

  !-------------------------------------------------------------------------
  !
  !> Calculates surface temperature and salinity tracer relaxation
  !!   relaxation terms for tracer equation and surface fluxes are calculated
  !!   in addition to other surface tracer fluxes
  !!   surface tracer restoring is applied either in apply_surface_relaxation
  !!   or in adding surface relaxation fluxes to total forcing fluxes
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2014)
  !
  SUBROUTINE update_surface_relaxation(p_patch_3D, p_os, p_ice, p_oce_sfc, tracer_no)

    TYPE (t_patch_3D ),    TARGET, INTENT(IN) :: p_patch_3D
    TYPE (t_hydro_ocean_state), INTENT(INOUT) :: p_os
    TYPE (t_sea_ice),              INTENT(IN) :: p_ice
    TYPE (t_ocean_surface)                    :: p_oce_sfc
    INTEGER,                       INTENT(IN) :: tracer_no       !  no of tracer: 1=temperature, 2=salinity

    !Local variables
    INTEGER                       :: jc, jb
    INTEGER                       :: i_startidx_c, i_endidx_c
    REAL(wp)                      :: relax_strength, thick
    TYPE(t_patch), POINTER        :: p_patch
    REAL(wp),      POINTER        :: t_top(:,:), s_top(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    !-------------------------------------------------------------------------

    t_top => p_os%p_prog(nold(1))%tracer(:,1,:,1)


    IF (tracer_no == 1) THEN  ! surface temperature relaxation
      !
      ! Temperature relaxation activated as additonal forcing term in the tracer equation
      ! implemented as simple time-dependent relaxation (time needed to restore tracer completely back to T*)
      !    F_T  = Q_T/dz = -1/tau * (T-T*) [ K/s ]  (where Q_T is boundary condition for vertical diffusion in [K*m/s])
      ! when using the sign convention
      !   dT/dt = Operators + F_T
      ! i.e. F_T <0 for  T-T* >0 (i.e. decreasing temperature T if T is warmer than relaxation data T*)

      ! EFFECTIVE RESTORING PARAMETER: 1.0_wp/(para_surfRelax_Temp*seconds_per_month)

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN

            !relax_strength = (p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) + p_os%p_prog(nold(1))%h(jc,jb)) / &
            !  &       (para_surfRelax_Temp*seconds_per_month)
!           p_oce_sfc%TopBC_Temp_vdiff(jc,jb) = -relax_strength*(t_top(jc,jb)-p_oce_sfc%data_surfRelax_Temp(jc,jb))
            relax_strength = 1.0_wp / (para_surfRelax_Temp*seconds_per_month)

            ! calculate additional temperature restoring rate F_T due to relaxation [K/s]
            p_oce_sfc%TempFlux_Relax(jc,jb) = -relax_strength*(t_top(jc,jb)-p_oce_sfc%data_surfRelax_Temp(jc,jb))

            ! Diagnosed heat flux Q_surf due to relaxation
            !  Q_surf = F_T*dz * (rho*Cp) = -dz/tau*(T-T*) * (rho*Cp)  [W/m2]
            !  HeatFlux_Relax = thick * TempFlux_Relax * (OceanReferenceDensity*clw)
            ! this heat flux is negative if relaxation flux is negative, i.e. heat is released if temperature decreases
            ! this flux is for diagnosis only and not added to tracer forcing

            thick = (p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)+p_os%p_prog(nold(1))%h(jc,jb))
            p_oce_sfc%HeatFlux_Relax(jc,jb) = p_oce_sfc%TempFlux_Relax(jc,jb) * thick * OceanReferenceDensity*clw

          ENDIF

        END DO
      END DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('UpdSfcRlx:HeatFlx_Rlx[W/m2]',p_oce_sfc%HeatFlux_Relax     ,str_module,2, in_subset=p_patch%cells%owned)
      CALL dbg_print('UpdSfcRlx: T* to relax to'  ,p_oce_sfc%data_surfRelax_Temp,str_module,4, in_subset=p_patch%cells%owned)
      CALL dbg_print('UpdSfcRlx: 1/tau*(T*-T)'    ,p_oce_sfc%TempFlux_Relax     ,str_module,3, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

    ELSE IF (tracer_no == 2) THEN  ! surface salinity relaxation
      !
      ! Salinity relaxation activated as additonal forcing term in the tracer equation
      ! implemented as simple time-dependent relaxation (time needed to restore tracer completely back to S*)
      !    F_S  = -1/tau * (S-S*) [ psu/s ]
      ! when using the sign convention
      !   dS/dt = Operators + F_S
      ! i.e. F_S <0 for  S-S* >0 (i.e. decreasing salinity S if S is saltier than relaxation data S*)
      ! note that the freshwater flux is opposite in sign to F_S, see below,
      ! i.e. fwf >0 for  S-S* >0 (i.e. increasing freshwater flux to decrease salinity)

      s_top => p_os%p_prog(nold(1))%tracer(:,1,:,2)

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN

            relax_strength = 1.0_wp / (para_surfRelax_Salt*seconds_per_month)
            !
            ! If sea ice is present (and l_relaxsal_ice), salinity relaxation is proportional to open water,
            !   under sea ice, no relaxation is applied, according to the procedure in MPIOM
            IF (l_relaxsal_ice .AND. i_sea_ice >=1) relax_strength = (1.0_wp-p_ice%concsum(jc,jb))*relax_strength

            ! calculate additional salt restoring rate F_S due to relaxation [psu/s]
            p_oce_sfc%SaltFlux_Relax(jc,jb) = -relax_strength*(s_top(jc,jb)-p_oce_sfc%data_surfRelax_Salt(jc,jb))

            ! Diagnosed freshwater flux due to relaxation (equivalent to heat flux Q)
            !  Fw_S = F_S*dz/S = dz/tau * (S-S*)/S  [m/s]
            ! this flux is applied as volume forcing in surface equation in fill_rhs4surface_eq_ab
            thick = (p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)+p_os%p_prog(nold(1))%h(jc,jb))
            p_oce_sfc%FrshFlux_Relax(jc,jb) = -p_oce_sfc%SaltFlux_Relax(jc,jb) * thick / s_top(jc,jb)

          ENDIF
        END DO
      END DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('UpdSfcRlx:FrshFlxRelax[m/s]',p_oce_sfc%FrshFlux_Relax     ,str_module,2, in_subset=p_patch%cells%owned)
      CALL dbg_print('UpdSfcRlx: S* to relax to'  ,p_oce_sfc%data_surfRelax_Salt,str_module,4, in_subset=p_patch%cells%owned)
      CALL dbg_print('UpdSfcRlx: 1/tau*(S*-S)'    ,p_oce_sfc%SaltFlux_Relax     ,str_module,3, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

    END IF  ! tracer_no

  END SUBROUTINE update_surface_relaxation

  !-------------------------------------------------------------------------
  !
  !> Calculates surface temperature and salinity tracer relaxation
  !!   relaxation terms for tracer equation and surface fluxes are calculated
  !!   in addition to other surface tracer fluxes
  !!   surface tracer restoring is applied either in apply_surface_relaxation
  !!   or in adding surface relaxation fluxes to total forcing fluxes
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2014)
  !
  SUBROUTINE apply_surface_relaxation(p_patch_3D, p_os, p_oce_sfc, tracer_no)

    TYPE (t_patch_3D ),    TARGET, INTENT(IN)    :: p_patch_3D
    TYPE (t_hydro_ocean_state),    INTENT(INOUT) :: p_os
    TYPE (t_ocean_surface), INTENT(IN)           :: p_oce_sfc
    INTEGER,                      INTENT(IN)     :: tracer_no       !  no of tracer: 1=temperature, 2=salinity

    !Local variables
    INTEGER :: jc, jb
    INTEGER :: i_startidx_c, i_endidx_c
    REAL(wp) :: t_top_old  (nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: s_top_old  (nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_patch), POINTER :: p_patch
    REAL(wp),      POINTER :: t_top(:,:), s_top(:,:)
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------

    all_cells => p_patch%cells%all

    t_top =>p_os%p_prog(nold(1))%tracer(:,1,:,1)
    t_top_old(:,:) = t_top(:,:)


    ! add relaxation term to temperature tracer
    IF (tracer_no == 1) THEN

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c

          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN
            t_top_old(jc,jb) = t_top(jc,jb)
            t_top(jc,jb)     = t_top_old(jc,jb) + p_oce_sfc%TempFlux_Relax(jc,jb)*dtime
          ENDIF

        END DO
      END DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('AppTrcRlx: TempFluxRelax'  , p_oce_sfc%TempFlux_Relax, str_module, 3, in_subset=p_patch%cells%owned)
      CALL dbg_print('AppTrcRlx: Old Temperature', t_top_old                  , str_module, 3, in_subset=p_patch%cells%owned)
      CALL dbg_print('AppTrcRlx: New Temperature', t_top                      , str_module, 2, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

    ! add relaxation term to salinity tracer
    ELSE IF (tracer_no == 2) THEN

      s_top =>p_os%p_prog(nold(1))%tracer(:,1,:,2)
      s_top_old(:,:) = s_top(:,:)

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN
            s_top(jc,jb)     = s_top_old(jc,jb) + p_oce_sfc%SaltFlux_Relax(jc,jb)*dtime
          ENDIF
        END DO
      END DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('AppTrcRlx: SaltFluxRelax', p_oce_sfc%SaltFlux_Relax, str_module, 3, in_subset=p_patch%cells%owned)
      CALL dbg_print('AppTrcRlx: Old Salt'     , s_top_old                  , str_module, 3, in_subset=p_patch%cells%owned)
      CALL dbg_print('AppTrcRlx: New Salt'     , s_top                      , str_module, 2, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

    END IF  ! tracer_no

  END SUBROUTINE apply_surface_relaxation

!**********************************************************************
!-------------------------------- OMIP --------------------------------
!**********************************************************************

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
  SUBROUTINE update_flux_fromFile(p_patch_3D, p_as, this_datetime)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_atmos_for_ocean)                     :: p_as
    TYPE(datetime), POINTER                     :: this_datetime
    !
    ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_bulk_forcing:update_flux_fromFile'
    INTEGER  :: jmon, jdmon, jmon1, jmon2, ylen, yday
    REAL(wp) :: rday1, rday2
    REAL(wp) ::  z_c2(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: sodt


    TYPE(t_patch), POINTER:: p_patch 
    !TYPE(t_subset_range), POINTER :: all_cells

    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------

    !all_cells       => p_patch%cells%all

    !  calculate day and month
    jmon  = this_datetime%date%month
    jdmon = this_datetime%date%day
    yday  = getDayOfYearFromDateTime(this_datetime)
    ylen  = getNoOfDaysInYearDateTime(this_datetime)

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
      ! - forcing data sets are read in mo_ext_data , forcing_timescale allocates and reads the no. of forcing steps
      ! - forcing_frequency is a namelist variable and controls how often the same forcing step is used  
      ! - jmon1 for controling correct forcing step
      ! - no time interpolation applied (jom2 = jmon1, rday2 = 0)
      ! - sotd = seconds of this day

      sodt=REAL(this_datetime%time%hour*3600._wp + &
           this_datetime%time%minute*60._wp + this_datetime%time%second)

      IF (forcing_timescale == 28 .OR. forcing_timescale == 29 .OR. forcing_timescale == 30 .OR. forcing_timescale == 31 )  THEN
        jmon1 = 1 + ((jdmon-1) * 86400.0_wp/forcing_frequency) + INT( sodt / forcing_frequency )
      ELSE IF (forcing_timescale == 28*24 .OR. forcing_timescale == 29*24  &
                               .OR. forcing_timescale == 30*24 .OR. forcing_timescale == 31*24 )  THEN
        jmon1 = 1 + ((jdmon-1) * 86400.0_wp/forcing_frequency) + INT( sodt / forcing_frequency )
      ELSE
        jmon1 = 1 + ((yday-1) * 86400.0_wp/forcing_frequency) + INT( sodt / forcing_frequency )
      ENDIF

      idt_src = 10
      IF ((my_process_is_stdio()) .AND. (idbg_mxmn >= idt_src)) &
      & write(0,"(a,i6,a,3i8)")' use forcing record ',jmon1,' at ', yday, this_datetime%time%hour,this_datetime%time%minute &
                  ,this_datetime%time%second

      jmon2 = jmon1
      rday1 = 1.0_wp
      rday2 = 0.0_wp


      ! Leap year in OMIP forcing: read Feb, 28 twice since only 365 data-sets are available
      IF (ylen == 366 .and. forcing_timescale == 365 ) then
        IF (yday>59) jmon1=yday-1
        jmon2=jmon1
      ENDIF




    END IF

    !
    ! OMIP data read in mo_ext_data into variable ext_data
    !

    ! file based wind forcing:
    ! provide OMIP fluxes for wind stress forcing
    ! data set 1:  wind_u(:,:)   !  'stress_x': zonal wind stress       [Pa]
    ! data set 2:  wind_v(:,:)   !  'stress_y': meridional wind stress  [Pa]
    !  - forcing_windstress_u_type and v_type not used anymore
    !  - full OMIP data read if iforc_oce=OMIP_FluxFromFile (=11)

    ! ext_data has rank n_dom due to grid refinement in the atmosphere but not in the ocean
    !IF (forcing_windstress_u_type == 1)
    p_as%topBoundCond_windStress_u(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,1) + &
      &                                   rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,1)

    !IF (forcing_windstress_v_type == 1) THEN
    p_as%topBoundCond_windStress_v(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,2) + &
      &                                   rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,2)


    !-------------------------------------------------------------------------
    ! provide OMIP fluxes for sea ice (interface to ocean)
    ! data set 4:  tafo(:,:),   &  ! 2 m air temperature                              [C]
    ! data set 5:  ftdew(:,:),  &  ! 2 m dew-point temperature                        [K]
    ! data set 6:  fu10(:,:) ,  &  ! 10 m wind speed                                  [m/s]
    ! data set 7:  fclou(:,:),  &  ! Fractional cloud cover
    ! data set 8:  pao(:,:),    &  ! Surface atmospheric pressure                     [hPa]
    ! data set 9:  fswr(:,:),   &  ! Incoming surface solar radiation                 [W/m]
    ! data set 10:  precip(:,:), &  ! precipitation rate                              [m/s]
    ! data set 11:  evap  (:,:), &  ! evaporation   rate                              [m/s]
    ! data set 12:  runoff(:,:)     ! river runoff  rate                              [m/s]
    ! data set 13: u(:,:),      &  ! 10m zonal wind speed                             [m/s]
    ! data set 14: v(:,:),      &  ! 10m meridional wind speed                        [m/s]

    !IF (iforc_type == 2 .OR. iforc_type == 5) THEN
    !IF (forcing_fluxes_type > 0 .AND. forcing_fluxes_type < 101 ) THEN
    !  - forcing_fluxes_type = 1 not used anymore,
    !  - full OMIP data read if iforc_oce=OMIP_FluxFromFile (=11)

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
    p_as%FrshFlux_Precipitation(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,10) + &
      &                                     rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,10)
    !p_as%FrshFlux_Evaporation  (:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,11) + &
    !  &                                     rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,11)
    IF (forcing_set_runoff_to_zero) THEN
      p_as%FrshFlux_Runoff(:,:) = 0.0_wp
    ELSE
      p_as%FrshFlux_Runoff(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,12) + &
        &                              rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,12)
    ENDIF

 !  ! for test only - introduced temporarily
 !  p_as%tafo(:,:)  = 292.9_wp
 !  !  - change units to deg C, subtract tmelt (0 deg C, 273.15)
 !  p_as%tafo(:,:)  = p_as%tafo(:,:) - 273.15
 !  p_as%ftdew(:,:) = 289.877
 !  p_as%fu10(:,:)  = 7.84831
 !  p_as%fclou(:,:) = 0.897972
 !  p_as%fswr(:,:)  = 289.489
 !  p_as%u(:,:)     = 0.0_wp
 !  p_as%v(:,:)     = 0.0_wp
 !  p_as%topBoundCond_windStress_u(:,:) = 0.0_wp
 !  p_as%topBoundCond_windStress_v(:,:) = 0.0_wp
 !  p_as%FrshFlux_Precipitation(:,:) = 1.04634e-8
 !  p_as%FrshFlux_Runoff(:,:) = 0.0_wp
 !  p_as%pao(:,:)   = 101300.0_wp

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,4)
    CALL dbg_print('FlxFil: Ext data4-ta/mon1' ,z_c2        ,str_module,3, in_subset=p_patch%cells%owned)
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,4)
    CALL dbg_print('FlxFil: Ext data4-ta/mon2' ,z_c2        ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('FlxFil: p_as%tafo'         ,p_as%tafo   ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('FlxFil: p_as%windStr-u',p_as%topBoundCond_windStress_u, str_module,3,in_subset=p_patch%cells%owned)
    CALL dbg_print('FlxFil: p_as%windStr-v',p_as%topBoundCond_windStress_v, str_module,4,in_subset=p_patch%cells%owned)
    CALL dbg_print('FlxFil: Precipitation' ,p_as%FrshFlux_Precipitation,str_module,3,in_subset=p_patch%cells%owned)
    CALL dbg_print('FlxFil: Runoff'        ,p_as%FrshFlux_Runoff       ,str_module,3,in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    IF (type_surfRelax_Temp == 2)  THEN

      !-------------------------------------------------------------------------
      ! Apply temperature relaxation data (record 3) from stationary forcing
      !  - change units to deg C, subtract tmelt (0 deg C, 273.15)
      !  - this is not done for type_surfRelax_Temp=3, since init-data is in Celsius

       p_as%data_surfRelax_Temp(:,:) = &
         &  rday1*(ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,3)-tmelt) + &
         &  rday2*(ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,3)-tmelt)

    END IF

    IF (type_surfRelax_Salt == 2 .AND. no_tracer >1) THEN

      !-------------------------------------------------------------------------
      ! Apply salinity relaxation data (record ??) from stationary forcing
      CALL finish(TRIM(ROUTINE),' type_surfRelax_Salt=2 (reading from flux file) not yet implemented')

    END IF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    IF (idbg_mxmn >= idt_src) THEN
      WRITE(message_text,'(a,2(a,i4),2(a,f12.8))') 'FLUX time interpolation:', &
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

  END SUBROUTINE update_flux_fromFile


  !-------------------------------------------------------------------------
  !
  !> Calc_omip_budgets_ice equals sbr "Budget" in MPIOM.
  !! Sets the atmospheric fluxes over *SEA ICE ONLY* for the update of the ice
  !! temperature and ice growth rates for OMIP forcing
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-07). Originally code written by
  !! Dirk Notz, following MPIOM. Code transfered to ICON.
  !! Einar Olason, split calc_atm_fluxes_from_bulk into calc_bulk_flux_ice and calc_bulk_flux_oce
  !! so that the ocean model can be run without the ice model, but with OMIP fluxes.
  !!
  !! Rewritten by Stephan Lorenz, MPI-M (2015-06).
  !!  Using interface with parameters in order to call budget routine independent of ocean model

  SUBROUTINE calc_omip_budgets_ice(geolat, tafoC, ftdewC, fu10, fclou, pao, fswr,                &
    &                              kice, tice, hice, albvisdir, albvisdif, albnirdir, albnirdif, &
    &                              LWnetIce, SWnetIce, sensIce, latentIce,                       &
    &                              dLWdTIce, dsensdTIce, dlatdTIce)                              

 !  INPUT variables for OMIP via parameter:
    REAL(wp), INTENT(in)    :: geolat(:,:)      ! latitude                             [rad]
    REAL(wp), INTENT(in)    :: tafoC(:,:)       ! 2 m air temperature in Celsius       [C]
    REAL(wp), INTENT(in)    :: ftdewC(:,:)      ! 2 m dew point temperature in Celsius [C]
    REAL(wp), INTENT(in)    :: fu10(:,:)        ! 10 m wind speed                      [m/s]
    REAL(wp), INTENT(in)    :: fclou(:,:)       ! Fractional cloud cover               [frac]
    REAL(wp), INTENT(in)    :: pao(:,:)         ! Surface atmospheric pressure         [hPa]
    REAL(wp), INTENT(in)    :: fswr(:,:)        ! Incoming surface solar radiation     [W/m2]
    INTEGER,  INTENT(in)    :: kice             ! number of ice classes (currently 1)
    REAL(wp), INTENT(in)    :: tice(:,:,:)      ! surface ice temperature per class    [C]
    REAL(wp), INTENT(in)    :: hice(:,:,:)      ! ice thickness per class              [m]
    REAL(wp), INTENT(in)    :: albvisdir(:,:,:) ! direct ice albedo per class
    REAL(wp), INTENT(in)    :: albvisdif(:,:,:) ! diffuse ice albedo per class
    REAL(wp), INTENT(in)    :: albnirdir(:,:,:) ! direct near infrared ice albedo per class
    REAL(wp), INTENT(in)    :: albnirdif(:,:,:) ! diffuse near infrared ice albedo per class

 !  OUTPUT variables for sea ice model via parameter (inout since icefree part is not touched)
    REAL(wp), INTENT(inout) :: LWnetIce (:,:,:) ! net longwave heat flux over ice      [W/m2]
    REAL(wp), INTENT(inout) :: SWnetIce (:,:,:) ! net shortwave heat flux over ice     [W/m2]
    REAL(wp), INTENT(inout) :: sensIce  (:,:,:) ! sensible heat flux over ice          [W/m2]
    REAL(wp), INTENT(inout) :: latentIce(:,:,:) ! latent heat flux over ice            [W/m2]
    REAL(wp), INTENT(inout) :: dLWdTIce (:,:,:) ! derivitave of LWnetIce w.r.t temperature
    REAL(wp), INTENT(inout) :: dsensdTIce(:,:,:)! derivitave of sensIce w.r.t temperature
    REAL(wp), INTENT(inout) :: dlatdTIce(:,:,:) ! derivitave of latentIce w.r.t temperature

 !  Local variables
 !  REAL(wp), DIMENSION (nproma,p_patch%alloc_cell_blocks) :: &
    REAL(wp), DIMENSION (SIZE(tafoC,1), SIZE(tafoC,2)) ::  &
      & Tsurf,          &  ! Surface temperature in Celsius                  [C]
      & tafoK,          &  ! Air temperature at 2 m in Kelvin                [K]
      & fu10lim,        &  ! wind speed at 10 m height in range 2.5...32     [m/s]
      & esta,           &  ! water vapor pressure at 2 m height              [Pa]
      & esti,           &  ! water vapor pressure at ice surface             [Pa]
      & sphumida,       &  ! Specific humididty at 2 m height
      & sphumidi,       &  ! Specific humididty at ice surface
      & rhoair,         &  ! air density                                     [kg/m^3]
      & dragl0,         &  ! part of dragl
      & dragl1,         &  ! part of dragl
      & dragl,          &  ! Drag coefficient for latent   heat flux
      & drags,          &  ! Drag coefficient for sensible heat flux (=0.95 dragl)
      & fakts,          &  ! Effect of cloudiness on LW radiation
      & humi,           &  ! Effect of air humidity on LW radiation
      & fa, fi,         &  ! Enhancment factor for vapor pressure
      & dsphumididesti, &  ! Derivative of sphumidi w.r.t. esti
      & destidT,        &  ! Derivative of esti w.r.t. T
      & dfdT               ! Derivative of f w.r.t. T
 !    & wspeed             ! Wind speed                                      [m/s]

    INTEGER :: i
    REAL(wp) :: aw,bw,cw,dw,ai,bi,ci,di,AAw,BBw,CCw,AAi,BBi,CCi,alpha,beta
    REAL(wp) :: fvisdir, fvisdif, fnirdir, fnirdif, local_rad2deg

    tafoK(:,:)  = tafoC(:,:) + tmelt  ! Change units of tafo  to Kelvin

    ! set to zero for NAG
    sphumida(:,:)  = 0.0_wp
    fa      (:,:)  = 0.0_wp
    esta    (:,:)  = 0.0_wp
    rhoair  (:,:)  = 0.0_wp

    !-----------------------------------------------------------------------
    ! Compute water vapor pressure and specific humididty in 2m height (esta)
    ! and at water surface (estw) according to "Buck Research Manual (1996)
    ! (see manuals for instruments at http://www.buck-research.com/);
    ! updated from Buck, A. L., New equations for computing vapor pressure and
    ! enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981"
    !-----------------------------------------------------------------------
    ! #slo# 2015-03: the comment above is now valid 
    ! the values for ice are not changed in Buck (1996) in comparison to Buck (1981)

    ! the following commented values are from Buck (1981)
    ! aw=611.21_wp; bw=18.729_wp; cw=257.87_wp; dw=227.3_wp

    ! here are the updated values for open water according to Buck (1996)
    ai=611.15_wp; bi=23.036_wp; ci=279.82_wp; di=333.7_wp
    aw=611.21_wp; bw=18.678_wp; cw=257.14_wp; dw=234.5_wp

    AAw=7.2e-4_wp; BBw=3.20e-6_wp; CCw=5.9e-10_wp
    AAi=2.2e-4_wp; BBi=3.83e-6_wp; CCi=6.4e-10_wp

    alpha=0.62197_wp; beta=0.37803_wp

    ! #slo# correction: pressure in enhancement formula is in mb (hPa) according to Buck 1981 and 1996
    fa(:,:)        = 1.0_wp+AAw+pao*0.01_wp*(BBw+CCw*ftdewC**2)
    esta(:,:)      = fa * aw*EXP((bw-ftdewC/dw)*ftdewC/(ftdewC+cw))
    sphumida(:,:)  = alpha * esta/(pao-beta*esta)

    !-----------------------------------------------------------------------
    !  Compute longwave radiation according to
    !         Berliand, M. E., and T. G. Berliand, 1952: Determining the net
    !         long-wave radiation of the Earth with consideration of the effect
    !         of cloudiness. Izv. Akad. Nauk SSSR, Ser. Geofiz., 1, 6478.
    !         cited by: Budyko, Climate and Life, 1974.
    !         Note that for humi, esta is given in [mmHg] in the original
    !         publication. Therefore, 0.05*sqrt(esta/100) is used rather than
    !         0.058*sqrt(esta)
    !  This is the formula used in MPI-OM when using the QLOBERL preprocessing option (currently
    !  the default usage).
    !-----------------------------------------------------------------------

    ! Berliand & Berliand ('52) calculate only LWnet
    humi    = 0.39_wp - 0.05_wp*SQRT(esta/100._wp)

    ! icon-identical calculation of rad2deg
    local_rad2deg = 180.0_wp / 3.14159265358979323846264338327950288_wp
    fakts   =  1.0_wp - ( 0.5_wp + 0.4_wp/90._wp &
      &         *MIN(ABS(local_rad2deg*geolat(:,:)),60._wp) ) * fclou(:,:)**2
  !   &         *MIN(ABS(rad2deg*p_patch%cells%center(:,:)%lat),60._wp) ) * fclou(:,:)**2

    !-----------------------------------------------------------------------
    !  Calculate bulk equations according to
    !      Kara, B. A., P. A. Rochford, and H. E. Hurlburt, 2002:
    !      Air-Sea Flux Estimates And The 19971998 Enso Event,  Bound.-Lay.
    !      Met., 103(3), 439-458, doi: 10.1023/A:1014945408605.
    !-----------------------------------------------------------------------

    ! with nag there is floating invalid operation on rest of last nproma-block only due to pao=nan
    ! rhoair(:,:) = pao(:,:) / (rd*tafoK(:,:)*(1.0_wp+0.61_wp*sphumida(:,:)) ) !  error with nag
    WHERE (pao(:,:)>0.0_wp) rhoair(:,:) = pao(:,:) / (rd*tafoK(:,:)*(1.0_wp+0.61_wp*sphumida(:,:)) )

    fu10lim(:,:)    = MAX (2.5_wp, MIN(32.5_wp,fu10(:,:)) )
    dragl1(:,:)     = 1e-3_wp*(-0.0154_wp + 0.5698_wp/fu10lim(:,:) &
      &               - 0.6743_wp/(fu10lim(:,:) * fu10lim(:,:)))
    dragl0(:,:)     = 1e-3_wp*(0.8195_wp+0.0506_wp*fu10lim(:,:) &
      &               - 0.0009_wp*fu10lim(:,:)*fu10lim(:,:))

    ! Fractions of SWin in each band (from cice)
    fvisdir=0.28_wp; fvisdif=0.24_wp; fnirdir=0.31_wp; fnirdif=0.17_wp
    Tsurf(:,:) = 0.0_wp ! For debug output

    ! Over sea ice area only
    !  TODO: in case of no ice model, ice variables cannot be used here
    !  ice classes: currently one class (kice=1) is used, therefore formulation can be simplified to 2-dim variables as in mpiom
    DO i = 1, kice
      WHERE (hice(:,i,:)>0._wp)
        
        !  albedo model: atmos_fluxes%albvisdir, albvisdif, albnirdir, albnirdif
        !   - all 4 albedos are the same (i_ice_albedo = 1), they are calculated in ice_fast and should be stored in p_ice
        SWnetIce(:,i,:)  = ( 1._wp-albvisdir(:,i,:) )*fvisdir*fswr(:,:) +   &
          &                ( 1._wp-albvisdif(:,i,:) )*fvisdif*fswr(:,:) +   &
          &                ( 1._wp-albnirdir(:,i,:) )*fnirdir*fswr(:,:) +   &
          &                ( 1._wp-albnirdif(:,i,:) )*fnirdif*fswr(:,:)
      ! Tsurf(:,:)       = p_ice%Tsurf(:,i,:)
        Tsurf(:,:)       = tice(:,i,:)
        ! pressure in enhancement formula is in mb (hPa) according to Buck 1981 and 1996
        fi(:,:)          = 1.0_wp+AAi+pao(:,:)*0.01_wp*(BBi+CCi*Tsurf(:,:) **2)
        esti(:,:)        = fi(:,:)*ai*EXP((bi-Tsurf(:,:) /di)*Tsurf(:,:) /(Tsurf(:,:) +ci))
        sphumidi(:,:)    = alpha*esti(:,:)/(pao(:,:)-beta*esti(:,:))
        ! This may not be the best drag parametrisation to use over ice
        dragl(:,:)       = dragl0(:,:) + dragl1(:,:) * (Tsurf(:,:)-tafoC(:,:))
        ! A reasonableee maximum and minimum is needed for dragl in case there's a large difference
        ! between the 2-m and surface temperatures.
        dragl(:,:)       = MAX(0.5e-3_wp, MIN(3.0e-3_wp,dragl(:,:)))
        drags(:,:)       = 0.95_wp * dragl(:,:)

        LWnetIce(:,i,:)  = -fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4 &
           &               -4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - tafoC(:,:))
        ! same form as MPIOM:
        !atmos_fluxes%LWnet (:,i,:)  = - (fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4 &
        !  &         + 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:)))
        dLWdTIce(:,i,:)  = -4._wp*zemiss_def*stbo*tafoK(:,:)**3
        sensIce(:,i,:)   = drags(:,:) * rhoair(:,:)*cpd*fu10(:,:) * fr_fac * (tafoC(:,:) -Tsurf(:,:))
        latentIce(:,i,:) = dragl(:,:) * rhoair(:,:)* alf *fu10(:,:) * fr_fac &
          &                   * (sphumida(:,:)-sphumidi(:,:))

        dsensdTIce(:,i,:)   = 0.95_wp*cpd*rhoair(:,:)*fu10(:,:)&
          &                   *(dragl0(:,:) - 2.0_wp*dragl(:,:))
        dsphumididesti(:,:) = alpha/(pao(:,:)-beta*esti(:,:)) &
          &                   * (1.0_wp + beta*esti(:,:)/(pao(:,:)-beta*esti(:,:)))
        destidT(:,:)        = (bi*ci*di-Tsurf(:,:)*(2.0_wp*ci+Tsurf(:,:)))&
          &                   /(di*(ci+Tsurf(:,:))**2) * esti(:,:)
        dfdT(:,:)           = 2.0_wp*CCi*BBi*Tsurf(:,:)
        dlatdTIce(:,i,:)    = alf*rhoair(:,:)*fu10(:,:)* &
          &                  ( (sphumida(:,:)-sphumidi(:,:))*dragl1(:,:) &
          &                    - dragl(:,:)*dsphumididesti(:,:)*(fi(:,:)*destidT(:,:) &
          &                    + esti(:,:)*dfdT(:,:)) )
      ENDWHERE
    ENDDO


  END SUBROUTINE calc_omip_budgets_ice

  !-------------------------------------------------------------------------
  !
  !> Forcing_from_bulk equals sbr "Budget_omip" in MPIOM.
  !! Sets the atmospheric fluxes for the update of
  !! temperature of *OPEN WATER* for OMIP forcing.
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2012-08). Originally code written by
  !! Dirk Notz, following MPIOM. Code transfered to ICON.
   
  SUBROUTINE calc_omip_budgets_oce(p_patch, p_as, p_os, atmos_fluxes)
    TYPE(t_patch),            INTENT(IN), TARGET    :: p_patch
    TYPE(t_atmos_for_ocean),  INTENT(IN)    :: p_as
    TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os
    TYPE(t_atmos_fluxes),     INTENT(INOUT) :: atmos_fluxes

 !  INPUT variables:
 !  p_as%tafo(:,:),      : 2 m air temperature                              [C]
 !  p_as%ftdew(:,:),     : 2 m dew-point temperature                        [K]
 !  p_as%fu10(:,:) ,     : 10 m wind speed                                  [m/s]
 !  p_as%fclou(:,:),     : Fractional cloud cover
 !  p_as%pao(:,:),       : Surface atmospheric pressure                     [hPa]
 !  p_as%fswr(:,:),      : Incoming surface solar radiation                 [W/m]
 !  p_os%tracer(:,1,:,1) : SST
 !  atmos_fluxes%albvisdirw, albvisdifw, albnirdirw, albnirdifw
 !
 !  OUTPUT variables:  atmos_fluxes - heat fluxes and wind stress over open ocean
 !  atmos_fluxes%LWnetw   : long wave
 !  atmos_fluxes%SWnetw   : short wave
 !  atmos_fluxes%sensw    : sensible
 !  atmos_fluxes%latw     : latent
 !  atmos_fluxes%stress_xw: zonal stress, water
 !  atmos_fluxes%stress_yw: meridional stress, water
 !  atmos_fluxes%stress_x : zonal stress, ice
 !  atmos_fluxes%stress_y : meridional stress, ice

 !  Local variables
    REAL(wp), DIMENSION (nproma,p_patch%alloc_cell_blocks) ::           &
      & Tsurf,          &  ! Surface temperature                             [C]
      & tafoK,          &  ! Air temperature at 2 m in Kelvin                [K]
      & fu10lim,        &  ! wind speed at 10 m height in range 2.5...32     [m/s]
      & esta,           &  ! water vapor pressure at 2 m height              [Pa]
      & estw,           &  ! water vapor pressure at water surface           [Pa]
      & sphumida,       &  ! Specific humididty at 2 m height
      & sphumidw,       &  ! Specific humididty at water surface
      & ftdewC,         &  ! Dew point temperature in Celsius                [C]
      & rhoair,         &  ! air density                                     [kg/m^3]
      & dragl0,         &  ! part of dragl
      & dragl1,         &  ! part of dragl
      & dragl,          &  ! Drag coefficient for latent   heat flux
      & drags,          &  ! Drag coefficient for sensible heat flux (=0.95 dragl)
      & fakts,          &  ! Effect of cloudiness on LW radiation
      & humi,           &  ! Effect of air humidity on LW radiation
      & fa, fw,         &  ! Enhancment factor for vapor pressure
      & wspeed,         &  ! Wind speed                                      [m/s]
      & C_ao               ! Drag coefficient for atm-ocean stress           [m/s]

    INTEGER :: jb, jc, i_startidx_c, i_endidx_c
    REAL(wp) :: aw,bw,cw,dw,AAw,BBw,CCw,alpha,beta
    REAL(wp) :: fvisdir, fvisdif, fnirdir, fnirdif

    TYPE(t_subset_range), POINTER :: all_cells

    Tsurf(:,:)  = p_os%p_prog(nold(1))%tracer(:,1,:,1)  ! set surface temp = mixed layer temp
    tafoK(:,:)  = p_as%tafo(:,:)  + tmelt               ! Change units of tafo  to Kelvin
    ftdewC(:,:) = p_as%ftdew(:,:) - tmelt               ! Change units of ftdew to Celsius

    ! subset range pointer
    all_cells => p_patch%cells%all

    !-----------------------------------------------------------------------
    ! Compute water vapor pressure and specific humididty in 2m height (esta)
    ! and at water surface (estw) according to "Buck Research Manual (1996)
    ! (see manuals for instruments at http://www.buck-research.com/);
    ! updated from Buck, A. L., New equations for computing vapor pressure and
    ! enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981"
    !-----------------------------------------------------------------------
    AAw   = 7.2e-4_wp;  BBw  = 3.20e-6_wp; CCw = 5.9e-10_wp
    alpha = 0.62197_wp; beta = 0.37803_wp
    ! #slo# 2015-03: the comment above is now valid - the following commented values are from Buck (1981)
    ! aw    = 611.21_wp; bw    = 18.729_wp;  cw  = 257.87_wp; dw = 227.3_wp
    ! these are the updated values according to Buck (1996)
    aw    = 611.21_wp; bw    = 18.678_wp;  cw  = 257.14_wp; dw = 234.5_wp

    ! #slo# correction: pressure in enhancement formula is in mb (hPa) according to Buck 1981 and 1996
   !fa(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*(BBw+CCw*ftdewC(:,:)**2)
    fa(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*0.01_wp*(BBw+CCw*ftdewC(:,:)**2)
    esta(:,:) = fa(:,:) * aw*EXP((bw-ftdewC(:,:)/dw)*ftdewC(:,:)/(ftdewC(:,:)+cw))
   !esta(:,:) =           aw*EXP((bw-ftdewC(:,:)/dw)*ftdewC(:,:)/(ftdewC(:,:)+cw))
   !fw(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*(BBw+CCw*Tsurf(:,:) **2)
    fw(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*0.01_wp*(BBw+CCw*Tsurf(:,:) **2)
   !estw(:,:) = fw(:,:) *aw*EXP((bw-Tsurf(:,:) /dw)*Tsurf(:,:) /(Tsurf(:,:) +cw))
    ! For a given surface salinity we should multiply estw with  1 - 0.000537*S
    ! #slo# correction according to MPIOM: lowering of saturation vapor pressure over saline water
    !       is taken constant to 0.9815
    estw(:,:) = 0.9815_wp*fw(:,:)*aw*EXP((bw-Tsurf(:,:) /dw)*Tsurf(:,:) /(Tsurf(:,:) +cw))

    sphumida(:,:)  = alpha * esta(:,:)/(p_as%pao(:,:)-beta*esta(:,:))
    sphumidw(:,:)  = alpha * estw(:,:)/(p_as%pao(:,:)-beta*estw(:,:))

    !-----------------------------------------------------------------------
    !  Compute longwave radiation according to
    !         Berliand, M. E., and T. G. Berliand, 1952: Determining the net
    !         long-wave radiation of the Earth with consideration of the effect
    !         of cloudiness. Izv. Akad. Nauk SSSR, Ser. Geofiz., 1, 6478.
    !         cited by: Budyko, Climate and Life, 1974.
    !         Note that for humi, esta is given in [mmHg] in the original
    !         publication. Therefore, 0.05*sqrt(esta/100) is used rather than
    !         0.058*sqrt(esta)
    !  This is the formula used in MPI-OM when using the QLOBERL preprocessing option (currently
    !  the default usage).
    !-----------------------------------------------------------------------

    humi(:,:)    = 0.39_wp - 0.05_wp*SQRT(esta(:,:)/100._wp)
    fakts(:,:)   =  1.0_wp - ( 0.5_wp + 0.4_wp/90._wp &
      &         *MIN(ABS(rad2deg*p_patch%cells%center(:,:)%lat),60._wp) ) * p_as%fclou(:,:)**2
    ! Berliand & Berliand ('52) calculate only LWnetw

    ! #eoo# 2012-12-14: another bugfix
    ! #slo# #hha# 2012-12-13: bugfix, corrected form
    atmos_fluxes%LWnetw(:,:) = - fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4  &
      &                - 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:))
    ! same form as MPIOM:
    !atmos_fluxes%LWnetw(:,:) = - (fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4  &
    !  &         + 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:)))
    ! bug
    !atmos_fluxes%LWnetw(:,:) = fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4  &
    !  &         - 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:))

    ! Fractions of SWin in each band (from cice)
    fvisdir=0.28_wp; fvisdif=0.24_wp; fnirdir=0.31_wp; fnirdif=0.17_wp
    atmos_fluxes%SWnetw(:,:) = ( 1._wp-atmos_fluxes%albvisdirw(:,:) )*fvisdir*p_as%fswr(:,:) +   &
      &                ( 1._wp-atmos_fluxes%albvisdifw(:,:) )*fvisdif*p_as%fswr(:,:) +   &
      &                ( 1._wp-atmos_fluxes%albnirdirw(:,:) )*fnirdir*p_as%fswr(:,:) +   &
      &                ( 1._wp-atmos_fluxes%albnirdifw(:,:) )*fnirdif*p_as%fswr(:,:)

    !-----------------------------------------------------------------------
    !  Calculate bulk equations according to
    !      Kara, B. A., P. A. Rochford, and H. E. Hurlburt, 2002:
    !      Air-Sea Flux Estimates And The 19971998 Enso Event,  Bound.-Lay.
    !      Met., 103(3), 439-458, doi: 10.1023/A:1014945408605.
    !-----------------------------------------------------------------------

    rhoair(:,:) = 0._wp
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c,i_endidx_c

        rhoair(jc,jb) = p_as%pao(jc,jb)                &
          &            /(rd*tafoK(jc,jb)*(1.0_wp+0.61_wp*sphumida(jc,jb)) )

      END DO
    END DO

    fu10lim(:,:)    = MAX (2.5_wp, MIN(32.5_wp,p_as%fu10(:,:)) )
    dragl1(:,:)     = 1e-3_wp*(-0.0154_wp + 0.5698_wp/fu10lim(:,:) &
      &               - 0.6743_wp/(fu10lim(:,:) * fu10lim(:,:)))
    dragl0(:,:)     = 1e-3_wp*(0.8195_wp+0.0506_wp*fu10lim(:,:) &
      &               - 0.0009_wp*fu10lim(:,:)*fu10lim(:,:))
    dragl(:,:)      = dragl0(:,:) + dragl1(:,:) * (Tsurf(:,:)-p_as%tafo(:,:))
    ! A reasonable maximum and minimum is needed for dragl in case there's a large difference
    ! between the 2-m and surface temperatures.
    dragl(:,:)      = MAX(0.5e-3_wp, MIN(3.0e-3_wp,dragl(:,:)))
    drags(:,:)      = 0.95_wp * dragl(:,:)
    atmos_fluxes%sensw(:,:) = drags(:,:)*rhoair(:,:)*cpd*p_as%fu10(:,:) * fr_fac &
      &               * (p_as%tafo(:,:) -Tsurf(:,:))
    atmos_fluxes%latw(:,:)  = dragl(:,:)*rhoair(:,:)*alv*p_as%fu10(:,:) * fr_fac &
      &               * (sphumida(:,:)-sphumidw(:,:))

    ! wind stress over ice and open water
    IF (use_calculated_ocean_stress) THEN
      !-----------------------------------------------------------------------
      !  Calculate oceanic wind stress according to:
      !   Gill (Atmosphere-Ocean Dynamics, 1982, Academic Press) (see also Smith, 1980, J. Phys
      !   Oceanogr., 10, 709-726)
      !-----------------------------------------------------------------------

      wspeed(:,:) = SQRT( p_as%u**2 + p_as%v**2 )
      C_ao(:,:)   = MIN( 2._wp, MAX(1.1_wp, 0.61_wp+0.063_wp*wspeed ) )*1e-3_wp
      atmos_fluxes%stress_xw(:,:) = C_ao(:,:)*rhoair(:,:)*wspeed(:,:)*p_as%u(:,:)! over water
      atmos_fluxes%stress_yw(:,:) = C_ao(:,:)*rhoair(:,:)*wspeed(:,:)*p_as%v(:,:)! over water

      atmos_fluxes%stress_x(:,:) = Cd_ia     *rhoair(:,:)*wspeed(:,:)*p_as%u(:,:)! over ice
      atmos_fluxes%stress_y(:,:) = Cd_ia     *rhoair(:,:)*wspeed(:,:)*p_as%v(:,:)! over ice
    ELSE
      ! use wind stress provided by OMIP data
      atmos_fluxes%stress_xw(:,:) = p_as%topBoundCond_windStress_u(:,:)
      atmos_fluxes%stress_yw(:,:) = p_as%topBoundCond_windStress_v(:,:)

      atmos_fluxes%stress_x(:,:) = p_as%topBoundCond_windStress_u(:,:) ! over ice
      atmos_fluxes%stress_y(:,:) = p_as%topBoundCond_windStress_v(:,:) ! over ice
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print level (1-5          , fix)
    CALL dbg_print('omipBudOce:tafoK'              , tafoK                 , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:tafo'               , p_as%tafo             , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:ftdew'              , p_as%ftdew            , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:ftdewC'             , ftdewC                , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:pao'                , p_as%pao              , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:fa'                 , fa                    , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:fw'                 , fw                    , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:esta'               , esta                  , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:estw'               , estw                  , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:sphumida'           , sphumida              , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:sphumidw'           , sphumidw              , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:rhoair'             , rhoair                , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:dragl'              , dragl                 , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:drags'              , drags                 , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:fu10'               , p_as%fu10             , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:fu10lim'            , fu10lim               , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:stress_xw'          , atmos_fluxes%stress_xw, str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:stress_yw'          , atmos_fluxes%stress_yw, str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:p_as%windStr-u',p_as%topBoundCond_windStress_u,str_module,idt_src, in_subset=p_patch%cells%owned)
    idt_src=3  ! output print level (1-5          , fix)
    CALL dbg_print('omipBudOce:Tsurf ocean'        , Tsurf                 , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:atmflx%SWnetw'      , atmos_fluxes%SWnetw   , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:atmflx%LWnetw'      , atmos_fluxes%LWnetw   , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:atmflx%sensw'       , atmos_fluxes%sensw    , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:atmflx%latw'        , atmos_fluxes%latw     , str_module, idt_src, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE calc_omip_budgets_oce

!**********************************************************************
!---------------------------- WIND STRESS -----------------------------
!**********************************************************************
  !-------------------------------------------------------------------------
  !
  !> Sets the surface stress the ocean sees because of the presence of sea ice.
  !! Calculated as the average of atm-ocean wind stress (atmos_fluxes%stress_xw)
  !! and ice-ocean stress (calculated in this routine, not stored).
  !! Wind-stress is either calculated in bulk-formula or from atmosphere via coupling.
  !!
  !! Note: this ice-ocean stress is calculated again in mo_ice_fem_evp on the FEM grid.
  !! ToDo: store ice-ocean stress, so that I-O and O-I stresses are guaranteed to match.
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-06-05)
  !! Rewritten by Vladimir Lapin, MPI-M (2017-02).

  ! difference of ice and ocean velocities determines ocean stress below sea ice
  ! resulting stress on ocean surface is stored in atmos_fluxes%topBoundCond_windStress_u

  SUBROUTINE update_ocean_surface_stress(p_patch_3D, p_ice, p_os, atmos_fluxes, p_oce_sfc)

    TYPE(t_patch_3d),TARGET,  INTENT(IN)    :: p_patch_3D
    TYPE(t_sea_ice),          INTENT(IN)    :: p_ice
    TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os
    TYPE (t_atmos_fluxes),    INTENT(IN)    :: atmos_fluxes
    TYPE (t_ocean_surface),   INTENT(INOUT)   :: p_oce_sfc

 !  INPUT variables:
 !  p_ice%u(:,:),               : zonal ice velocity on centre                  [m/s]
 !  p_ice%u(:,:),               : meridional ice velocity on centre             [m/s]
 !  p_os%p_diag%u(:,1,:)        : zonal ocean velocity on centre                [m/s]
 !  p_os%p_diag%v(:,1,:)        : meridional ocean velocity on centre           [m/s]
 !  atmos_fluxes%stress_xw(:,:) : zonal stress, water
 !  atmos_fluxes%stress_yw(:,:) : meridional stress, water
 !  p_ice%concSum(:,:)          : total ice concentration within a grid cell
 !
 !  OUTPUT variables:
 !  p_oce_sfc%TopBC_WindStress_u: zonal stress that the ocean receives
 !  p_oce_sfc%TopBC_WindStress_v: meridional stress that the ocean receives
 !  p_oce_sfc%TopBC_WindStress_cc: its cartesian components

 !  Local variables
    REAL(wp) :: drag_coeff
!    REAL(wp) :: C_io               ! Drag coefficient for ice-ocean stress      [m/s]
    REAL(wp) :: delu                ! Zonal I-O velocity mismatch                [m/s]
    REAL(wp) :: delv                ! Meridional I-O velocity mismatch           [m/s]
    REAL(wp) :: delabs              ! I-O mismatch abs value                     [m/s]
    ! Ranges
    TYPE(t_patch), POINTER        :: p_patch
    TYPE(t_subset_range), POINTER :: all_cells
    ! Indexing
    INTEGER  :: i_startidx_c, i_endidx_c, jc, jb

!--------------------------------------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
!--------------------------------------------------------------------------------------------------

    IF ( i_sea_ice > 0 ) THEN ! sea ice is on

      drag_coeff = rho_ref*Cd_io
      ! TODO: The ice-ocean drag coefficient should depend on the depth of the upper most ocean
      ! velocity point: Cd_io = ( kappa/log(z/z0) )**2, with z0 ~= 0.4 cm

      ! for runs without ice dynamics, set ocean-ice stress to zero (no deceleration below sea ice)
      IF (stress_ice_zero) drag_coeff = 0.0_wp

!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, jc, delu, delv, delabs) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          delu = p_ice%u(jc,jb) - p_os%p_diag%u(jc,1,jb)
          delv = p_ice%v(jc,jb) - p_os%p_diag%v(jc,1,jb)
          delabs = SQRT( delu**2 + delv**2 )

          p_oce_sfc%TopBC_WindStress_u(jc,jb) = atmos_fluxes%stress_xw(jc,jb)*( 1._wp - p_ice%concSum(jc,jb) )   &
            &                                       + drag_coeff*delabs*delu * p_ice%concSum(jc,jb)
          p_oce_sfc%TopBC_WindStress_v(jc,jb) = atmos_fluxes%stress_yw(jc,jb)*( 1._wp - p_ice%concSum(jc,jb) )   &
            &                                       + drag_coeff*delabs*delv * p_ice%concSum(jc,jb)
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ELSE   !  sea ice is off

      ! apply wind stress directly
      p_oce_sfc%TopBC_WindStress_u(:,:) = atmos_fluxes%stress_xw(:,:)
      p_oce_sfc%TopBC_WindStress_v(:,:) = atmos_fluxes%stress_yw(:,:)

    ENDIF

!--------------------------------------------------------------------------------------------------

    ! After final updating of zonal and merdional components cartesian coordinates are calculated
!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
          CALL gvec2cvec(  p_oce_sfc%TopBC_WindStress_u(jc,jb),&
                         & p_oce_sfc%TopBC_WindStress_v(jc,jb),&
                         & p_patch%cells%center(jc,jb)%lon,&
                         & p_patch%cells%center(jc,jb)%lat,&
                         & p_oce_sfc%TopBC_WindStress_cc(jc,jb)%x(1),&
                         & p_oce_sfc%TopBC_WindStress_cc(jc,jb)%x(2),&
                         & p_oce_sfc%TopBC_WindStress_cc(jc,jb)%x(3))
        ELSE
          p_oce_sfc%TopBC_WindStress_u(jc,jb)         = 0.0_wp
          p_oce_sfc%TopBC_WindStress_v(jc,jb)         = 0.0_wp
          p_oce_sfc%TopBC_WindStress_cc(jc,jb)%x      = 0.0_wp
        ENDIF
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('sfc_flx: windStress_u',p_oce_sfc%TopBC_WindStress_u, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('sfc_flx: windStress_v',p_oce_sfc%TopBC_WindStress_v, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('sfc_flx: windStress_cc1',p_oce_sfc%TopBC_WindStress_cc%x(1), &
      &             str_module,3, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE update_ocean_surface_stress

!**********************************************************************
!---------------------------- OTHER -----------------------------
!**********************************************************************

  !-------------------------------------------------------------------------
  !>
  !! Balance sea level to zero over global ocean
  !!
  !! Balance sea level to zero over global ocean
  !! This routine uses parts of mo_ocean_diagnostics
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


END MODULE mo_ocean_bulk_forcing
