!>
!! @brief Contains subroutines for initializing the ECHAM physics
!! package in ICOHAM.
!!
!! @author Hui Wan, MPI-M 
!!
!! @par Revision History
!! First version by Hui Wan, 2010-07-20
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_echam_phy_init

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish, message, message_text

  USE mo_sync,                 ONLY: sync_c, sync_patch_array

  USE mo_netcdf_read,          ONLY: netcdf_read_oncells_2D

  ! model configuration
  USE mo_dynamics_config,      ONLY: nnow 
  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: nlev, iqv, iqt, ntracer
  USE mo_vertical_coord_table, ONLY: vct
  USE mo_echam_phy_config,     ONLY: phy_config => echam_phy_config, &
                                   & configure_echam_phy
  USE mo_echam_conv_config,    ONLY: configure_echam_convection

  ! test cases
  USE mo_ha_testcases,         ONLY: ape_sst_case
  USE mo_ape_params,           ONLY: ape_sst
  USE mo_physical_constants,   ONLY: tmelt, Tf, albi, albedoW

!   USE mo_math_utilities,      ONLY: sphere_cell_mean_char_length

  ! radiation
  USE mo_radiation_config,     ONLY: ssi, tsi
  USE mo_srtm_config,          ONLY: setup_srtm, ssi_amip
  USE mo_lrtm_setup,           ONLY: lrtm_setup
  USE mo_newcld_optics,        ONLY: setup_newcld_optics

  ! vertical diffusion
  USE mo_echam_vdiff_params,   ONLY: init_vdiff_params, z0m_min
  USE mo_vdiff_solver,         ONLY: init_vdiff_solver

  ! cumulus convection
  USE mo_convect_tables,       ONLY: init_convect_tables

  ! stratiform clouds and cloud cover
  USE mo_echam_cloud_params,   ONLY: init_cloud_tables, sucloud, cvarmin

  ! air-sea-land interface
  USE mo_ext_data_state,       ONLY: ext_data
  USE mo_icoham_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd, &
                                   & init_sfc_indices

  ! domain and indices
  USE mo_model_domain,         ONLY: t_patch
  USE mo_loopindices,          ONLY: get_indices_c

  ! atmospheric state
  USE mo_icoham_dyn_types,     ONLY: t_hydro_atm
  USE mo_eta_coord_diag,       ONLY: half_level_pressure, full_level_pressure
  USE mo_echam_phy_memory,     ONLY: construct_echam_phy_state,    &
                                   & prm_field, t_echam_phy_field, &
                                   & prm_tend,  t_echam_phy_tend,  &
                                   & mean_charlen
  ! for coupling
  USE mo_coupling_config,      ONLY: is_coupled_run
  USE mo_icon_cpl_exchg,       ONLY: ICON_cpl_get_init, ICON_cpl_put_init
  USE mo_icon_cpl_def_field,   ONLY: ICON_cpl_get_nbr_fields, ICON_cpl_get_field_ids
  USE mo_timer,                ONLY: ltimer, timers_level, timer_start, timer_stop, &
    & timer_prep_echam_phy


  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: init_echam_phy, initcond_echam_phy
  PUBLIC  :: additional_restart_init

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !! Top-level routine for initialization of ECHAM6 physics.
  !! It calls a series of subroutines to initialize tunable parameters,
  !! lookup tables, and the physics state vectors "prm_field" and "prm_tend".
  !!
  !! @par Revision History
  !! Initial version by Hui Wan, MPI-M (2010-07)
  !! name change to init_echam_phy by Levi Silvers
  !!
  SUBROUTINE init_echam_phy( p_patch, ltestcase, ctest_name, &
                                nlev, vct_a, vct_b, ceta        )

    TYPE(t_patch),   INTENT(IN) :: p_patch(:)
    LOGICAL,         INTENT(IN) :: ltestcase
    CHARACTER(LEN=*),INTENT(IN) :: ctest_name
    INTEGER,         INTENT(IN) :: nlev
    REAL(wp),        INTENT(IN) :: vct_a(:), vct_b(:), ceta(:)

    INTEGER :: khydromet, ktrac
    INTEGER :: jg, ndomain


    IF (timers_level > 1) CALL timer_start(timer_prep_echam_phy)

    !-------------------------------------------------------------------
    ! Initialize parameters and lookup tables
    !-------------------------------------------------------------------
    ! Main switches (phy_config%lrad, phy_config%lcond, etc.)

    CALL configure_echam_phy (ltestcase, ctest_name)

    ! For radiation:

    IF (phy_config%lrad) THEN
      ssi(:) = ssi_amip(:)
      tsi    = SUM(ssi(:))
      CALL setup_srtm
      CALL lrtm_setup
      CALL setup_newcld_optics
    END IF

    ! For cumulus convection: 
    ! - assign value to echam_conv_config%nmctop;
    ! - allocate echam_conv_config%cevapcu(:) and assign values.

    IF (phy_config%lconv) THEN
      CALL configure_echam_convection(nlev, vct_a, vct_b, ceta)
    END IF

    ! For surface processes: 
    ! nsfc_type, iwtr, etc. are set in this subroutine. 
    ! See mo_icoham_sfc_indicies.f90 for further details.

    CALL init_sfc_indices( ltestcase, ctest_name )

    ! For turbulent mixing:
    ! Allocate memory for the tri-diagonal solver needed by the implicit
    ! time stepping scheme; Compute time-independent parameters.

    IF (phy_config%lvdiff) THEN
      ! Currently the tracer indices are sorted such that we count
      ! the water substances first, and then other species like 
      ! aerosols and their precursors. "ntracer" is the total number 
      ! of tracers (including water substances) handled in the model;
      ! "iqt" is the starting index for non-water species.
      ! Before more sophisticated meta-data structure becomes available, 
      ! it is assumed here that all tracers are subject to turbulent mixing.

      khydromet = iqt - 2        ! # of hydrometeors
      ktrac = ntracer - iqt + 1  ! # of non-water species 

      CALL init_vdiff_solver( khydromet, ktrac, nlev )
      CALL init_vdiff_params( nlev, nlev+1, nlev+1, vct )
    ENDIF

    ! Lookup tables for saturation vapour pressure

    IF (phy_config%lconv.OR.phy_config%lcond.OR.phy_config%lvdiff) &
    CALL init_convect_tables

    ! For large scale condensation:

    IF (phy_config%lcond) THEN
      CALL init_cloud_tables
      CALL sucloud( nlev, vct        &
!!$        &         , lmidatm=.FALSE.  &
        &         , lcouple=.FALSE.  &
        &         , lipcc=.FALSE.    &
!!$        &         , lham=.FALSE.     &
        &         )
    END IF

    !-------------------------------------------------------------------
    ! Allocate memory for the state vectors "prm_field" and "prm_tend"
    !-------------------------------------------------------------------
    CALL construct_echam_phy_state( ntracer, p_patch )


    ! general 
   !--------------------------------------------------------------
    !< characteristic gridlength needed by sso and sometimes by
    !! convection and turbulence
    !--------------------------------------------------------------
    ndomain = SIZE(p_patch)
 
    DO jg= 1,ndomain
      ! read it directly from the patch%geometry_info
      mean_charlen(jg) = p_patch(jg)%geometry_info%mean_characteristic_length
!       CALL sphere_cell_mean_char_length (p_patch(jg)%n_patch_cells_g, &
!         & mean_charlen(jg))
    ENDDO
    
    IF (timers_level > 1) CALL timer_stop(timer_prep_echam_phy)

  END SUBROUTINE init_echam_phy
  !-------------
  !>
  !! Loop over all grid levels and give proper values to some components
  !! of the state vectors "prm_field" and "prm_tend".
  !! This subroutine plays a role similar to "init_g3" in ECHAM6.
  !!
  !! @par Revision History
  !! Initial version by Hui Wan, MPI-M (2010-07)
  !!
  SUBROUTINE initcond_echam_phy( p_patch, p_hydro_state, ltestcase, ctest_name )

    TYPE(t_patch)    ,INTENT(IN) :: p_patch(:)
    TYPE(t_hydro_atm),INTENT(IN) :: p_hydro_state(:)
    LOGICAL,          INTENT(IN) :: ltestcase
    CHARACTER(LEN=*), INTENT(IN) :: ctest_name

    ! local variables and pointers

    INTEGER  :: ndomain, nblks_c, jg, jb, jbs, jc, jcs, jce, jk
    REAL(wp) :: zprat, zn1, zn2, zcdnc, zlat
    LOGICAL  :: lland, lglac

    TYPE(t_echam_phy_field),POINTER :: field => NULL()
    TYPE(t_echam_phy_tend) ,POINTER :: tend  => NULL()
    !----

    ! in case of coupling
    INTEGER               :: nbr_fields
    INTEGER               :: nbr_points
    INTEGER               :: nbr_hor_points

    INTEGER               :: field_shape(3)
    INTEGER, ALLOCATABLE  :: field_id(:)
    REAL(wp), ALLOCATABLE :: buffer(:,:)

    INTEGER               :: info, ierror !< return values form cpl_put/get calls
    INTEGER               :: return_status

    CHARACTER(len=*), PARAMETER :: slm_fn = 'bc_slm.nc'
    LOGICAL       :: lexist

    ! total number of domains/ grid levels

    ndomain = SIZE(prm_field)
    IF (ndomain.eq.0) CALL finish('init_phy_memory', &
       & 'ERROR: array prm_field has zero length')

    !-------------------------
    ! Loop over all domains
    !-------------------------
    DO jg = 1,ndomain

      field => prm_field(jg)
      tend  => prm_tend (jg)

      !----------------------------------------
      ! Loop over all blocks in domain jg
      !----------------------------------------
      nblks_c = p_patch(jg)%nblks_int_c
      jbs     = p_patch(jg)%cells%start_blk(2,1)

        ! For idealized test cases

      IF (ltestcase) THEN
        SELECT CASE (ctest_name)
        CASE('AMIP') !Note that there are only two surface types in this case
                    ! water and land, as long as no ice routines are implemented 
                    ! in the ECHAM-physics

!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
          INQUIRE (file=slm_fn, exist=lexist)
          IF (lexist) THEN
                return_status = netcdf_read_oncells_2D(slm_fn,'slm', field% lsmask, p_patch(jg)) 
          ELSE
            WRITE (message_text,*) 'Could not open file ',slm_fn
            CALL message('',message_text)
            CALL finish ('initcond_echam_phy:read slm ', 'run terminated.')
          ENDIF
          DO jb = jbs,nblks_c
            CALL get_indices_c( p_patch(jg), jb,jbs,nblks_c, jcs,jce, 2)
            field% tsfc_tile(jcs:jce,jb,iwtr) = tmelt
            field% tsfc_tile(jcs:jce,jb,ilnd) = tmelt
            field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
            field% seaice(jcs:jce,jb) = 0._wp   ! zero sea ice fraction
          END DO
!$OMP END PARALLEL DO

        CASE('APE') !Note that there is only one surface type in this case
                    !except ljsbach=.true. with land and ocean (two surface types)

          IF (phy_config%ljsbach) THEN
!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
            DO jb = jbs,nblks_c
              CALL get_indices_c( p_patch(jg), jb,jbs,nblks_c, jcs,jce, 2)
              DO jc = jcs,jce
                zlat = p_patch(jg)%cells%center(jc,jb)%lat
!                field% tsfc_tile(jc,jb,ilnd) = ext_data(jg)%atm%sst(jc,jb)   ! SST (preliminary)
!                field% tsfc_tile(jc,jb,iwtr) = ext_data(jg)%atm%sst(jc,jb)   ! SST
                field% tsfc_tile(jc,jb,ilnd) = ext_data(jg)%atm%sst_mon(jc,1,jb)   ! SST (preliminary)
                field% tsfc_tile(jc,jb,iwtr) = ext_data(jg)%atm%sst_mon(jc,1,jb)   ! SST
                IF (ext_data(jg)%atm%lsm_ctr_c(jc,jb) > 0) THEN ! this is land
                  field% tsfc(jc,jb) = field% tsfc_tile(jc,jb,ilnd)
                  field% lsmask(jc,jb) = 1._wp   ! a grid box is either completely land ...
                ELSE ! this is ocean
                  field% tsfc(jc,jb) = field% tsfc_tile(jc,jb,iwtr)
                  field% lsmask(jc,jb) = 0._wp   ! ... or a grid box is completely ocean
                END IF
              END DO
              field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
              field% seaice(jcs:jce,jb) = 0._wp   ! zero sea ice fraction
            END DO
!$OMP END PARALLEL DO
          ELSE
!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
            DO jb = jbs,nblks_c
              CALL get_indices_c( p_patch(jg), jb,jbs,nblks_c, jcs,jce, 2)
              DO jc = jcs,jce
                zlat = p_patch(jg)%cells%center(jc,jb)%lat
               !field% tsfc_tile(jc,iwtr,jb) = ape_sst(ape_sst_case,zlat)   ! SST
               !field% tsfc     (jc,     jb) = field% tsfc_tile(jc,iwtr,jb)
                field% tsfc_tile(jc,jb,iwtr) = ape_sst(ape_sst_case,zlat)   ! SST
                field% tsfc     (jc,     jb) = field% tsfc_tile(jc,jb,iwtr)
              END DO
              field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
              field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
              field% seaice(jcs:jce,jb) = 0._wp   ! zeor sea ice fraction
            END DO
!$OMP END PARALLEL DO
          END IF ! ljsbach

          IF ( is_coupled_run() ) CALL finish('ERROR: Use testcase APEc for a coupled run')

        CASE('APEi') 
          ! The same as APE, except we don't allow JSBACH and whenever SST reaches tmelt we put 1 m
          ! thick ice with concentration of 0.9 on top

!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
          DO jb = jbs,nblks_c
            CALL get_indices_c( p_patch(jg), jb,jbs,nblks_c, jcs,jce, 2)
            DO jc = jcs,jce
              zlat = p_patch(jg)%cells%center(jc,jb)%lat
              ! SST must reach Tf where there's ice. It may be better to modify ape_sst it self.
              field% tsfc_tile(jc,jb,iwtr) = ape_sst(ape_sst_case,zlat) + Tf
              ! Initialise the ice - Tsurf, T1 & T2 must be in degC
              field% tsfc_tile  (jc,jb,iice) = Tf + tmelt
              field% Tsurf      (jc,1, jb  ) = Tf
              field% T1         (jc,1, jb  ) = Tf
              field% T2         (jc,1, jb  ) = Tf
              field% hs         (jc,1, jb  ) = 0._wp
              IF ( field%tsfc_tile(jc,jb,iwtr) <= Tf + tmelt ) THEN
!                ! Set the ice surface temperature to the same value as the lowest model level above
!                ! surface. This is copied from the JWw and LDF cases.
!                field%tsfc_tile(jc,jb,iice) = &
!                  &     p_hydro_state(jg)%prog(nnow(jg))%temp(jc,nlev,jb)

                field%Tsurf (jc,1,jb) = field% tsfc_tile(jc,jb,iice) - tmelt
                field%conc  (jc,1,jb) = 0.9_wp
                field%hi    (jc,1,jb) = 1.0_wp
                field%seaice(jc,  jb) = field%conc(jc,1,jb)
              ELSE
                field%conc  (jc,1,jb) = 0._wp
                field%hi    (jc,1,jb) = 0._wp
                field%seaice(jc,  jb) = field%conc(jc,1,jb)
              ENDIF
              field% tsfc(jc,jb) = field%seaice(jc,jb)*field%tsfc_tile(jc,jb,iice) &
                &       + ( 1._wp - field%seaice(jc,jb) )*field%tsfc_tile(jc,jb,iwtr)
            END DO
            field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
            field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          END DO
!$OMP END PARALLEL DO
          field% albvisdir_ice(:,:,:) = albi ! albedo in the visible range for direct radiation
          field% albnirdir_ice(:,:,:) = albi ! albedo in the NIR range for direct radiation 
          field% albvisdif_ice(:,:,:) = albi ! albedo in the visible range for diffuse radiation
          field% albnirdif_ice(:,:,:) = albi ! albedo in the NIR range for diffuse radiation
          field% albvisdir_wtr(:,:) = albedoW ! albedo in the visible range for direct radiation
          field% albnirdir_wtr(:,:) = albedoW ! albedo in the NIR range for direct radiation 
          field% albvisdif_wtr(:,:) = albedoW ! ! albedo in the visible range for diffuse radiation
          field% albnirdif_wtr(:,:) = albedoW ! albedo in the NIR range for diffuse radiation

        CASE('APEc') 
          ! The same as APEi, except we initialize with no ice and don't modify the surface
          ! temperature. This is meant for a coupled run.

!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
          DO jb = jbs,nblks_c
            CALL get_indices_c( p_patch(jg), jb,jbs,nblks_c, jcs,jce, 2)
            DO jc = jcs,jce
              zlat = p_patch(jg)%cells%center(jc,jb)%lat
              field% tsfc_tile(jc,jb,iwtr) = ape_sst(ape_sst_case,zlat)
              ! Initialise the ice - Tsurf, T1 & T2 must be in degC
              field% tsfc_tile  (jc,jb,iice) = Tf + tmelt
              field% Tsurf      (jc,1, jb  ) = Tf
              field% T1         (jc,1, jb  ) = Tf
              field% T2         (jc,1, jb  ) = Tf
              field% hs         (jc,1, jb  ) = 0._wp
              field%conc  (jc,1,jb) = 0._wp
              field%hi    (jc,1,jb) = 0._wp
              field%seaice(jc,  jb) = field%conc(jc,1,jb)
              field% tsfc(jc,jb) = field%seaice(jc,jb)*field%tsfc_tile(jc,jb,iice) &
                &       + ( 1._wp - field%seaice(jc,jb) )*field%tsfc_tile(jc,jb,iwtr)
            END DO
            field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
            field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          END DO
!$OMP END PARALLEL DO
          field% albvisdir_ice(:,:,:) = albi ! albedo in the visible range for direct radiation
          field% albnirdir_ice(:,:,:) = albi ! albedo in the NIR range for direct radiation 
          field% albvisdif_ice(:,:,:) = albi ! albedo in the visible range for diffuse radiation
          field% albnirdif_ice(:,:,:) = albi ! albedo in the NIR range for diffuse radiation
          field% albvisdir_wtr(:,:) = albedoW ! albedo in the visible range for direct radiation
          field% albnirdir_wtr(:,:) = albedoW ! albedo in the NIR range for direct radiation 
          field% albvisdif_wtr(:,:) = albedoW ! ! albedo in the visible range for diffuse radiation
          field% albnirdif_wtr(:,:) = albedoW ! albedo in the NIR range for diffuse radiation

! This shouldn't be necessary!
          IF ( is_coupled_run() ) THEN

             ALLOCATE(buffer(nproma*nblks_c,4))

             nbr_hor_points = p_patch(jg)%n_patch_cells
             nbr_points     = nproma * p_patch(jg)%nblks_c

             !
             !  see drivers/mo_atmo_model.f90:
             !
             !   field_id(1) represents "TAUX"   wind stress component
             !   field_id(2) represents "TAUY"   wind stress component
             !   field_id(3) represents "SFWFLX" surface fresh water flux
             !   field_id(4) represents "SFTEMP" surface temperature
             !   field_id(5) represents "THFLX"  total heat flux
             !   field_id(6) represents "ICEATM" ice temperatures and melt potential
             !
             !   field_id(7) represents "SST"    sea surface temperature
             !   field_id(9) represents "OCEANU" u component of ocean surface current
             !   field_id(9) represents "OCEANV" v component of ocean surface current
             !   field_id(10)represents "ICEOCE" ice thickness, concentration and temperatures
             !
             !
             CALL ICON_cpl_get_nbr_fields ( nbr_fields )
             ALLOCATE(field_id(nbr_fields))
             CALL ICON_cpl_get_field_ids ( nbr_fields, field_id )
             !
             field_shape(1) = 1
             field_shape(2) = nbr_hor_points
             field_shape(3) = 1
             !
             ! Send fields away 
             ! ----------------
             !
             ! Is there really anything to send or can the ocean live without?
             !
             ! Send fields away
             ! ----------------
             !
             ! TAUX
             !
             buffer(:,1) = RESHAPE ( field%u_stress_tile(:,:,iwtr), (/ nbr_points /) )
             CALL ICON_cpl_put_init ( field_id(1), field_shape, &
                                      buffer(1:nbr_hor_points,1:1), ierror )
             !
             ! TAUY
             !
             buffer(:,1) = RESHAPE ( field%v_stress_tile(:,:,iwtr), (/ nbr_points /) )
             CALL ICON_cpl_put_init ( field_id(2), field_shape, &
                                      buffer(1:nbr_hor_points,1:1), ierror )
             !
             ! SFWFLX Note: the evap_tile should be properly updated and added
             !
             buffer(:,1) = RESHAPE ( field%rsfl(:,:), (/ nbr_points /) ) + &
                  &        RESHAPE ( field%rsfc(:,:), (/ nbr_points /) ) + &
                  &        RESHAPE ( field%ssfl(:,:), (/ nbr_points /) ) + &
                  &        RESHAPE ( field%ssfc(:,:), (/ nbr_points /) )
             buffer(:,2) = RESHAPE ( field%evap_tile(:,:,iwtr), (/ nbr_points /) )

             field_shape(3) = 2
             CALL ICON_cpl_put_init ( field_id(3), field_shape, &
                                      buffer(1:nbr_hor_points,1:2), ierror )
             !
             ! SFTEMP
             !
             buffer(:,1) =  RESHAPE ( field%temp(:,nlev,:), (/ nbr_points /) )
             field_shape(3) = 1
             CALL ICON_cpl_put_init ( field_id(4), field_shape, &
                                      buffer(1:nbr_hor_points,1:1), ierror )
             !
             ! THFLX, total heat flux
             !
             buffer(:,1) =  RESHAPE ( field%swflxsfc_tile(:,:,iwtr), (/ nbr_points /) ) !net shortwave flux for ocean
             buffer(:,2) =  RESHAPE ( field%lwflxsfc_tile(:,:,iwtr), (/ nbr_points /) ) + &
              &             RESHAPE ( field%shflx_tile(:,:,iwtr),    (/ nbr_points /) ) + &
              &             RESHAPE ( field%lhflx_tile(:,:,iwtr),    (/ nbr_points /) ) !net non-solar fluxes for ocean
             field_shape(3) = 2
             CALL ICON_cpl_put_init ( field_id(5), field_shape, &
                                      buffer(1:nbr_hor_points,1:2), ierror )
             !
             ! ICEATM, Ice state determined by atmosphere
             !
             buffer(:,1) =  RESHAPE ( field%Qtop(:,1,:), (/ nbr_points /) ) !Melt-potential for ice - top
             buffer(:,2) =  RESHAPE ( field%Qbot(:,1,:), (/ nbr_points /) ) !Melt-potential for ice - bottom
             buffer(:,3) =  RESHAPE ( field%T1  (:,1,:), (/ nbr_points /) ) !Temperature of upper ice layer
             buffer(:,4) =  RESHAPE ( field%T2  (:,1,:), (/ nbr_points /) ) !Temperature of lower ice layer
             field_shape(3) = 4
             CALL ICON_cpl_put_init ( field_id(6), field_shape, &
                                      buffer(1:nbr_hor_points,1:4), ierror )

             ! Receive fields, only assign values if something was received ( info > 0 )
             ! -------------------------------------------------------------------------
             !
             ! I guess that only the SST is really needed.
             !
             !
             ! SST
             !
             field_shape(3) = 1
             CALL ICON_cpl_get_init ( field_id(7), field_shape, &
                                      buffer(1:nbr_hor_points,1:1), info, ierror )

             IF ( info > 0 ) THEN
                buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
                field%tsfc_tile(:,:,iwtr) = RESHAPE (buffer(:,1), (/ nproma, nblks_c /) )
                field%tsfc     (:,:)      = field%tsfc_tile(:,:,iwtr)
                CALL sync_patch_array(sync_c, p_patch(jg), field%tsfc_tile(:,:,iwtr))
                CALL sync_patch_array(sync_c, p_patch(jg), field%tsfc     (:,:))
             ENDIF
             !
             ! OCEANU
             !
             CALL ICON_cpl_get_init ( field_id(8), field_shape, &
                                      buffer(1:nbr_hor_points,1:1), info, ierror )
             IF ( info > 0 ) THEN
                buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
                field%ocu(:,:) = RESHAPE (buffer(:,1), (/ nproma, nblks_c /) )
                CALL sync_patch_array(sync_c, p_patch(jg), field%ocu(:,:))
             ENDIF
             !
             ! OCEANV
             !
             CALL ICON_cpl_get_init ( field_id(9), field_shape, &
                                      buffer(1:nbr_hor_points,1:1), info, ierror )
             IF ( info > 0 ) THEN
                buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
                field%ocv(:,:) = RESHAPE (buffer(:,1), (/ nproma, nblks_c /) )
                CALL sync_patch_array(sync_c, p_patch(jg), field%ocv(:,:))
             ENDIF

             !
             ! ICEOCE
             !
             field_shape(3) = 4
             CALL ICON_cpl_get_init ( field_id(10), field_shape, &
                                      buffer(1:nbr_hor_points,1:4), info, ierror )
             IF ( info > 0 ) THEN
               buffer(nbr_hor_points+1:nbr_points,1:4) = 0.0_wp
               field%hi  (:,1,:) = RESHAPE (buffer(:,1), (/ nproma, nblks_c /) )
               field%conc(:,1,:) = RESHAPE (buffer(:,2), (/ nproma, nblks_c /) )
               field%T1  (:,1,:) = RESHAPE (buffer(:,3), (/ nproma, nblks_c /) )
               field%T2  (:,1,:) = RESHAPE (buffer(:,4), (/ nproma, nblks_c /) )
               field%seaice(:,:) = field%conc(:,1,:)
               CALL sync_patch_array(sync_c, p_patch(jg), field%hi  (:,1,:))
               CALL sync_patch_array(sync_c, p_patch(jg), field%conc(:,1,:))
               CALL sync_patch_array(sync_c, p_patch(jg), field%seaice(:,:))
               CALL sync_patch_array(sync_c, p_patch(jg), field%T1  (:,1,:))
               CALL sync_patch_array(sync_c, p_patch(jg), field%T2  (:,1,:))
             ENDIF

             DEALLOCATE(field_id)
             DEALLOCATE(buffer)

          ENDIF


        CASE('JWw-Moist','LDF-Moist')

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jcs,jce,jk,zlat,zprat,lland,lglac,zn1,zn2,zcdnc) ICON_OMP_DEFAULT_SCHEDULE
          DO jb = jbs,nblks_c
            CALL get_indices_c( p_patch(jg), jb,jbs,nblks_c, jcs,jce, 2)
           ! Set the surface temperature to the same value as the lowest model
            ! level above surface. For this test case, currently we assume
            ! there is no land or sea ice.

           !field% tsfc_tile(jcs:jce,iwtr,jb) = p_hydro_state(jg)%prog(nnow(jg))% &
           !                                  & temp(jcs:jce,nlev,jb)
           !field% tsfc     (jcs:jce,     jb) = field% tsfc_tile(jcs:jce,iwtr,jb)
            field% tsfc_tile(jcs:jce,jb,iwtr) = p_hydro_state(jg)%prog(nnow(jg))% &
                                              & temp(jcs:jce,nlev,jb)
            field% tsfc     (jcs:jce,     jb) = field% tsfc_tile(jcs:jce,jb,iwtr)

            field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
            field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
            field% seaice(jcs:jce,jb) = 0._wp   ! zero sea ice fraction
          END DO
!$OMP END DO  NOWAIT
!$OMP END PARALLEL
          
        END SELECT
      ENDIF ! ltestcase

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jcs,jce,jk,zlat,zprat,lland,lglac,zn1,zn2,zcdnc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,nblks_c
        CALL get_indices_c( p_patch(jg), jb,jbs,nblks_c, jcs,jce, 2)
        ! Compute pressure at half and full levels

        CALL half_level_pressure( p_hydro_state(jg)%prog(nnow(jg))%pres_sfc(:,jb), &! in
                                & nproma, jce,                          &! in
                                & field%presi_old(:,:,jb)               )! out

        CALL full_level_pressure( field%presi_old(:,:,jb), nproma, jce, &! in
                                & field%presm_old(:,:,jb)               )! out

        ! Initialize the flag lfland (.TRUE. if the fraction of land in 
        ! a grid box is larger than zero). In ECHAM a local array
        ! is initialized in each call of the subroutine "physc"

        DO jc = jcs,jce
          field%lfland(jc,jb) = field%lsmask(jc,jb).GT.0._wp
          field%lfglac(jc,jb) = field%glac  (jc,jb).GT.0._wp
          ! DWD NWP version
        ! field%lfland(jc,jb) = ext_data(jg)%atm%lsm_atm_c(jc,jb) > 0
        ! field%lfglac(jc,jb) = ext_data(jg)%atm%soiltyp  (jc,jb) == 1 ! soiltyp=ice
        ENDDO


        ! Initialize cloud droplet number concentration (acdnc)
        ! (In ECHAM6 this is done in subroutine "physc" using a 
        ! "IF (lstart) THEN" block.)

        DO jk = 1,nlev
          DO jc = jcs,jce
             zprat=(MIN(8._wp,80000._wp/field%presm_old(jc,jk,jb)))**2

             lland = field%lfland(jc,jb)
             lglac = lland.AND.field%glac(jc,jb).GT.0._wp
             IF (lland.AND.(.NOT.lglac)) THEN
               zn1= 50._wp
               zn2=220._wp
             ELSE
               zn1= 50._wp
               zn2= 80._wp
             ENDIF
             IF (field%presm_old(jc,jk,jb).LT.80000._wp) THEN
                zcdnc=1.e6_wp*(zn1+(zn2-zn1)*(EXP(1._wp-zprat)))
             ELSE
                zcdnc=zn2*1.e6_wp
             ENDIF
             field% acdnc(jc,jk,jb) = zcdnc
          END DO !jc
        END DO   !jk
      ENDDO      !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! Assign initial values for some components of the "field" and 
      ! "tend" state vectors.

!$OMP PARALLEL 
!$OMP WORKSHARE
      field% q(:,:,:,iqv)  = p_hydro_state(jg)%prog(nnow(jg))% tracer(:,:,:,iqv)
     !field% q(:,:,:,iqc)  = p_hydro_state(jg)%prog(nnow(jg))% tracer(:,:,:,iqc)
     !field% q(:,:,:,iqi)  = p_hydro_state(jg)%prog(nnow(jg))% tracer(:,:,:,iqi)
     !field% q(:,:,:,iqt:) = p_hydro_state(jg)%prog(nnow(jg))% tracer(:,:,:,iqt:)

      field% xvar  (:,:,:) = cvarmin*field% q(:,:,:,iqv)
      field% xskew (:,:,:) = 2._wp

      ! Other variabels (cf. subroutine init_g3 in ECHAM6)

      field% topmax(:,  :) = 99999._wp
      field% thvsig(:,  :) = 1.e-2_wp
      field% tke   (:,:,:) = 1.e-4_wp

      field% cosmu0    (:,  :) = 0._wp
      field% flxdwswtoa(:,  :) = 0._wp
      field% swflxsfc    (:,:) = 0._wp
      field% lwflxsfc    (:,:) = 0._wp
      field% dlwflxsfc_dT(:,:) = 0._wp
      field% swflxtoa    (:,:) = 0._wp
      field% lwflxtoa    (:,:) = 0._wp
      field% swflxsfc_avg(:,:) = 0._wp
      field% lwflxsfc_avg(:,:) = 0._wp
      field% swflxtoa_avg(:,:) = 0._wp
      field% lwflxtoa_avg(:,:) = 0._wp
      field% dlwflxsfc_dT_avg(:,:) = 0._wp
      field% aclc  (:,:,:) = 0._wp
      field% aclcac(:,:,:) = 0._wp
      field% aclcov(:,  :) = 0._wp
      field% qvi   (:,  :) = 0._wp
      field% xlvi  (:,  :) = 0._wp
      field% xivi  (:,  :) = 0._wp
      field% aprl  (:,  :) = 0._wp
      field% aprc  (:,  :) = 0._wp
      field% aprs  (:,  :) = 0._wp
      field% rsfl  (:,  :) = 0._wp
      field% ssfl  (:,  :) = 0._wp
      field% rsfc  (:,  :) = 0._wp
      field% ssfc  (:,  :) = 0._wp
      field% omega (:,:,:) = 0._wp

      field%totprec_avg(:,:) = 0._wp
      field%  evap_avg(:,  :) = 0._wp
      field% lhflx_avg(:,  :) = 0._wp
      field% shflx_avg(:,  :) = 0._wp
      field%dshflx_dT_tile    (:,:,:)= 0._wp
      field%dshflx_dT_avg_tile(:,:,:)= 0._wp

      field% u_stress_avg(:,  :) = 0._wp
      field% v_stress_avg(:,  :) = 0._wp

      field% rtype (:,  :) = 0._wp
      field% rintop(:,  :) = 0._wp

      field% albvisdir(:,  :) = 0.07_wp ! albedo in the visible range for direct radiation
                                             ! (set to the albedo of water for testing)
      field% albnirdir(:,  :) = 0.07_wp ! albedo in the NIR range for direct radiation 
                                             ! (set to the albedo of water for testing)
      field% albvisdif(:,  :) = 0.07_wp ! albedo in the visible range for diffuse radiation
                                             ! (set to the albedo of water for testing)
      field% albnirdif(:,  :) = 0.07_wp ! albedo in the NIR range for diffuse radiation
                                             ! (set to the albedo of water for testing)

      tend% x_dtr(:,:,:)   = 0._wp  !"xtec" in ECHAM
!$OMP END WORKSHARE

      IF (phy_config%ljsbach) THEN
!!$ TR for JSBACH testing (energy balance)
!$OMP WORKSHARE
      field% surface_temperature    (:,  :) = field% tsfc_tile(:,:,iwtr) !! initialize surface temperature == sst for testinf
      field% surface_temperature_old(:,  :) = field% tsfc_tile(:,:,iwtr) !! initialize surface temperature == sst for testing
      field% surface_temperature_rad(:,  :) = field% tsfc_tile(:,:,iwtr) !! initialize surface temperature == sst for testing
      field% surface_temperature_eff(:,  :) = field% tsfc_tile(:,:,iwtr) !! initialize surface temperature == sst for testing
      field% c_soil_temperature1    (:,  :) = 0._wp
      field% c_soil_temperature2    (:,  :) = 0._wp
      field% c_soil_temperature3    (:,  :) = 0._wp
      field% c_soil_temperature4    (:,  :) = 0._wp
      field% c_soil_temperature5    (:,  :) = 0._wp
      field% d_soil_temperature1    (:,  :) = 1._wp
      field% d_soil_temperature2    (:,  :) = 1._wp
      field% d_soil_temperature3    (:,  :) = 1._wp
      field% d_soil_temperature4    (:,  :) = 1._wp
      field% d_soil_temperature5    (:,  :) = 1._wp
      field% soil_temperature1      (:,  :) = field% tsfc_tile(:,:,iwtr) !! initialize soil temperature == sst for testindg
      field% soil_temperature2      (:,  :) = field% tsfc_tile(:,:,iwtr) !! initialize soil temperature == sst for testindg
      field% soil_temperature3      (:,  :) = field% tsfc_tile(:,:,iwtr) !! initialize soil temperature == sst for testindg
      field% soil_temperature4      (:,  :) = field% tsfc_tile(:,:,iwtr) !! initialize soil temperature == sst for testindg
      field% soil_temperature5      (:,  :) = field% tsfc_tile(:,:,iwtr) !! initialize soil temperature == sst for testindg
      field% heat_capacity          (:,  :) = 1._wp
      field% ground_heat_flux       (:,  :) = 0._wp
      field% swnet                  (:,  :) = 0._wp
      field% time_steps_soil        (:,  :) = 0._wp

!!$ TR for JSBACH testing (hydrology)
      field% moisture1              (:,  :) = 0.0109355_wp
      field% moisture2              (:,  :) = 0.0427327_wp
      field% moisture3              (:,  :) = 0.153602_wp
      field% moisture4              (:,  :) = 0.251354_wp
      field% moisture5              (:,  :) = 0._wp
      field% moisture_all           (:,  :) = 0.446431_wp
      field% csat                   (:,  :) = 1.0_wp
      field% cair                   (:,  :) = 1.0_wp
      field% csat_transpiration     (:,  :) = 0._wp
      field% sat_surface_specific_humidity (:,  :) = 0.01_wp
      field% skin_reservoir         (:,  :) = 0._wp
      field% snow_fract             (:,  :) = 0._wp
      field% snow                   (:,  :) = 0._wp
      field% snow_canopy            (:,  :) = 0._wp
      field% snow_melt              (:,  :) = 0._wp
      field% snow_acc               (:,  :) = 0._wp
      field% snow_melt_acc          (:,  :) = 0._wp
      field% glacier_runoff_acc     (:,  :) = 0._wp
      field% runoff_acc             (:,  :) = 0._wp
      field% drainage_acc           (:,  :) = 0._wp
!$OMP END WORKSHARE
      END IF ! ljsbach
!$OMP END PARALLEL

      IF (phy_config%lvdiff) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,nblks_c
        CALL get_indices_c( p_patch(jg), jb,jbs,nblks_c, jcs,jce, 2)
        field% coriol(jcs:jce,jb) = p_patch(jg)%cells%f_c(jcs:jce,jb)
      ENDDO 
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$OMP PARALLEL WORKSHARE
       !field% coriol(:,:)   = p_patch(jg)%cells%f_c(:,:)
        field% ustar (:,:)   = 1._wp
        field% kedisp(:,:)   = 0._wp
        field% tkem0 (:,:,:) = 1.e-4_wp
        field% tkem1 (:,:,:) = 1.e-4_wp
        field% thvvar(:,:,:) = 1.e-4_wp
        field% ocu   (:,:)   = 0._wp
        field% ocv   (:,:)   = 0._wp
        field% mixlen(:,:,:) = -999._wp
!$OMP END PARALLEL WORKSHARE
       !IF (iwtr<=nsfc_type) field% z0m_tile(:,iwtr,:) = 1e-3_wp !see init_surf in echam (or z0m_oce?)
       !IF (iice<=nsfc_type) field% z0m_tile(:,iice,:) = 1e-3_wp !see init_surf in echam (or z0m_ice?)
       !IF (ilnd<=nsfc_type) field% z0m_tile(:,ilnd,:) = z0m_min ! or maybe a larger value?
        IF (iwtr<=nsfc_type) field% z0m_tile(:,:,iwtr) = 1e-3_wp !see init_surf in echam (or z0m_oce?)
        IF (iice<=nsfc_type) field% z0m_tile(:,:,iice) = 1e-3_wp !see init_surf in echam (or z0m_ice?)
        IF (ilnd<=nsfc_type) field% z0m_tile(:,:,ilnd) = z0m_min ! or maybe a larger value?
      ENDIF

      ! Initialization of tendencies is necessary for doing I/O with
      ! the NAG compiler 

      tend% temp_radsw(:,:,:) = 0._wp
      tend% temp_radlw(:,:,:) = 0._wp

      tend% temp_cld(:,:,:)   = 0._wp
      tend%    q_cld(:,:,:,:) = 0._wp

      tend% temp_cnv(:,:,:)   = 0._wp
      tend%    q_cnv(:,:,:,:) = 0._wp
      tend%    u_cnv(:,:,:)   = 0._wp
      tend%    v_cnv(:,:,:)   = 0._wp

      tend% temp_vdf(:,:,:)   = 0._wp
      tend%    q_vdf(:,:,:,:) = 0._wp
      tend%    u_vdf(:,:,:)   = 0._wp
      tend%    v_vdf(:,:,:)   = 0._wp

      tend% temp_gwh(:,:,:)   = 0._wp
      tend%    u_gwh(:,:,:)   = 0._wp
      tend%    v_gwh(:,:,:)   = 0._wp

!!$      field% debug_2d_1(:,  :) = 0.0_wp
!!$      field% debug_2d_2(:,  :) = 0.0_wp
!!$      field% debug_2d_3(:,  :) = 0.0_wp
!!$      field% debug_2d_4(:,  :) = 0.0_wp
!!$      field% debug_2d_5(:,  :) = 0.0_wp
!!$      field% debug_2d_6(:,  :) = 0.0_wp
!!$      field% debug_2d_7(:,  :) = 0.0_wp
!!$      field% debug_2d_8(:,  :) = 0.0_wp
!!$
!!$      field% debug_3d_1(:,:,:) = 0.0_wp
!!$      field% debug_3d_2(:,:,:) = 0.0_wp
!!$      field% debug_3d_3(:,:,:) = 0.0_wp
!!$      field% debug_3d_4(:,:,:) = 0.0_wp
!!$      field% debug_3d_5(:,:,:) = 0.0_wp
!!$      field% debug_3d_6(:,:,:) = 0.0_wp
!!$      field% debug_3d_7(:,:,:) = 0.0_wp
!!$      field% debug_3d_8(:,:,:) = 0.0_wp

      NULLIFY( field,tend )
    ENDDO !domain loop

  END SUBROUTINE initcond_echam_phy
  !-------------
  !>
  !!
  SUBROUTINE additional_restart_init( p_patch, ltestcase, ctest_name )

    TYPE(t_patch),   INTENT(IN) :: p_patch(:)
    LOGICAL,         INTENT(IN) :: ltestcase
    CHARACTER(LEN=*),INTENT(IN) :: ctest_name

    INTEGER :: ndomain, nblks_c, jg, jb, jbs, jc, jcs, jce
    REAL(wp):: zlat

    TYPE(t_echam_phy_field),POINTER :: field => NULL()

    CHARACTER(LEN=*),PARAMETER :: routine = 'additional_restart_init'

    !----
    ! total number of domains/ grid levels

    ndomain = SIZE(prm_field)
    IF (ndomain.eq.0) CALL finish('init_phy_memory', &
       & 'ERROR: array prm_field has zero length')

    !-------------------------
    ! Loop over all domains
    !-------------------------
    DO jg = 1,ndomain

      field => prm_field(jg)

      nblks_c = p_patch(jg)%nblks_int_c
      jbs     = p_patch(jg)%cells%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,nblks_c
        CALL get_indices_c( p_patch(jg), jb,jbs,nblks_c, jcs,jce, 2)

        !---------------------------------------------------------------------
        ! Re-initialize SST, sea ice and glacier if necessary
        !---------------------------------------------------------------------
        IF (ltestcase) THEN

          SELECT CASE (ctest_name)
          CASE('APE')
          ! For an aqua-planet experiment, re-initialization is necessary if 
          ! the restart file in use was generated during a differently configured 
          ! experiment (e.g., an APE exp with a different SST setup, or 
          ! a real-world simulation such as AMIP, etc). 

            DO jc = jcs,jce
              zlat = p_patch(jg)%cells%center(jc,jb)%lat
             !field% tsfc_tile(jc,iwtr,jb) = ape_sst(ape_sst_case,zlat)   ! SST
             !field% tsfc     (jc,     jb) = field% tsfc_tile(jc,iwtr,jb)
              field% tsfc_tile(jc,jb,iwtr) = ape_sst(ape_sst_case,zlat)   ! SST
              field% tsfc     (jc,     jb) = field% tsfc_tile(jc,jb,iwtr)
            END DO
            field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
            field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
            field% seaice(jcs:jce,jb) = 0._wp   ! zeor sea ice fraction

          END SELECT

        ELSE 
          CALL finish(TRIM(routine),'ltestcase = .FALSE. '//                     &
                     & 'Implement re-initialization of SST, sea ice and glacier.')
        END IF

        !--------------------------------------------------------------------
        ! Initialize the flag lfland (.TRUE. if the fraction of land in 
        ! a grid box is larger than zero). In ECHAM a local array
        ! is initialized in each call of the subroutine "physc".
        ! Note that this initialization is needed for all resumed integrations
        ! regardless of the choice of "ltestcase" and "ctest_name", because
        ! logical variables can not yet be stored in restart files.
        !--------------------------------------------------------------------

        DO jc = jcs,jce
          field%lfland(jc,jb) = field%lsmask(jc,jb).GT.0._wp
          field%lfglac(jc,jb) = field%glac  (jc,jb).GT.0._wp
          ! DWD NWP version
        ! field%lfland(jc,jb) = ext_data(jg)%atm%lsm_atm_c(jc,jb) > 0
        ! field%lfglac(jc,jb) = ext_data(jg)%atm%soiltyp  (jc,jb) == 1 ! soiltyp=ice
        ENDDO !jc
      ENDDO   !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      !----------------------------------------
      ! Reset accumulated variables
      !----------------------------------------

      field% aclcac     (:,:,:) = 0._wp
      field% aclcov     (:,  :) = 0._wp
      field%totprec_avg(:,:)   = 0._wp

      field% qvi   (:,  :) = 0._wp
      field% xlvi  (:,  :) = 0._wp
      field% xivi  (:,  :) = 0._wp

      field% aprl  (:,  :) = 0._wp
      field% aprc  (:,  :) = 0._wp
      field% aprs  (:,  :) = 0._wp

      field%  evap_avg(:,:) = 0._wp
      field% lhflx_avg(:,:) = 0._wp
      field% shflx_avg(:,:) = 0._wp

      field% u_stress_avg(:,:) = 0._wp
      field% v_stress_avg(:,:) = 0._wp

     field% swflxsfc_avg(:,:) = 0._wp
     field% lwflxsfc_avg(:,:) = 0._wp
     field% swflxtoa_avg(:,:) = 0._wp
     field% lwflxtoa_avg(:,:) = 0._wp

      NULLIFY( field )
    ENDDO !jg

  END SUBROUTINE additional_restart_init

END MODULE mo_echam_phy_init

