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
MODULE mo_echam_phy_init

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish

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
                                   & prm_tend,  t_echam_phy_tend
  ! for coupling
  USE mo_master_control,       ONLY: is_coupled_run
  USE mo_icon_cpl_exchg,       ONLY: ICON_cpl_put, ICON_cpl_get
  USE mo_icon_cpl_def_field,   ONLY: ICON_cpl_get_nbr_fields, ICON_cpl_get_field_ids

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: prepare_echam_phy, initcond_echam_phy
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
  !!
  SUBROUTINE prepare_echam_phy( p_patch, ltestcase, ctest_name, &
                                nlev, vct_a, vct_b, ceta        )

    TYPE(t_patch),   INTENT(IN) :: p_patch(:)
    LOGICAL,         INTENT(IN) :: ltestcase
    CHARACTER(LEN=*),INTENT(IN) :: ctest_name
    INTEGER,         INTENT(IN) :: nlev
    REAL(wp),        INTENT(IN) :: vct_a(:), vct_b(:), ceta(:)

    INTEGER :: khydromet, ktrac

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

  END SUBROUTINE prepare_echam_phy
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
    INTEGER               :: field_shape(3)
    INTEGER, ALLOCATABLE  :: field_id(:)
    REAL(wp), ALLOCATABLE :: buffer(:,:)

    INTEGER               :: info, ierror !< return values form cpl_put/get calls

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
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jcs,jce,jk,zlat,zprat,lland,lglac,zn1,zn2,zcdnc)
      DO jb = jbs,nblks_c
        CALL get_indices_c( p_patch(jg), jb,jbs,nblks_c, jcs,jce, 2)

        ! For idealized test cases

        IF (ltestcase) THEN
          SELECT CASE (ctest_name)
          CASE('APE') !Note that there is only one surface type in this case

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

            IF ( is_coupled_run() ) THEN

               ALLOCATE(buffer(nproma*nblks_c,1))
               !
               !  see drivers/mo_atmo_model.f90:
               !
               !   field_id(1) represents "TAUX"   wind stress component
               !   field_id(2) represents "TAUY"   wind stress component
               !   field_id(3) represents "SFWFLX" surface fresh water flux
               !   field_id(4) represents "SHFLX"  sensible heat flux
               !   field_id(5) represents "LHFLX"  latent heat flux
               !
               !   field_id(6) represents "SST"    sea surface temperature
               !   field_id(7) represents "OCEANU" u component of ocean surface current
               !   field_id(8) represents "OCEANV" v component of ocean surface current
               !
               CALL ICON_cpl_get_nbr_fields ( nbr_fields )
               ALLOCATE(field_id(nbr_fields))
               CALL ICON_cpl_get_field_ids ( nbr_fields, field_id )
               !
               field_shape(1) = 1
               field_shape(2) = p_patch(jg)%n_patch_cells
               field_shape(3) = 1
               !
               ! Send fields away 
               ! ----------------
               !
               ! Is there really anything to send or can the ocean live without?
               !
               ! Receive fields, only assign values if something was received ( info > 0 )
               ! -------------------------------------------------------------------------
               !
               ! I guess that only the SST is really needed.
               !
               CALL ICON_cpl_get ( field_id(6), field_shape, buffer, info, ierror )
               IF ( info > 0 ) &
               field%tsfc_tile(:,:,iwtr) = RESHAPE (buffer(:,1), (/ nproma, nblks_c /) )

               DEALLOCATE(field_id)
               DEALLOCATE(buffer)

            ENDIF

          CASE('JWw-Moist','LDF-Moist')
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
          END SELECT
        ENDIF ! ltestcase

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
!$OMP END DO
!$OMP END PARALLEL

      ! Assign initial values for some components of the "field" and 
      ! "tend" state vectors.

!$OMP PARALLEL WORKSHARE
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

      field%  evap_ac(:,  :) = 0._wp
      field% lhflx_ac(:,  :) = 0._wp
      field% shflx_ac(:,  :) = 0._wp

      field% u_stress_ac(:,  :) = 0._wp
      field% v_stress_ac(:,  :) = 0._wp

      field% rtype (:,  :) = 0._wp
      field% rintop(:,  :) = 0._wp

       tend% x_dtr(:,:,:) = 0._wp  !"xtec" in ECHAM
!$OMP END PARALLEL WORKSHARE

      IF (phy_config%lvdiff) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce)
      DO jb = jbs,nblks_c
        CALL get_indices_c( p_patch(jg), jb,jbs,nblks_c, jcs,jce, 2)
        field% coriol(jcs:jce,jb) = p_patch(jg)%cells%f_c(jcs:jce,jb)
      ENDDO 
!$OMP END DO
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
!$OMP DO PRIVATE(jb,jc,jcs,jce,zlat)
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
!$OMP END DO
!$OMP END PARALLEL

      !----------------------------------------
      ! Reset accumulated variables
      !----------------------------------------

      field% aclcac(:,:,:) = 0._wp
      field% aclcov(:,  :) = 0._wp

      field% qvi   (:,  :) = 0._wp
      field% xlvi  (:,  :) = 0._wp
      field% xivi  (:,  :) = 0._wp

      field% aprl  (:,  :) = 0._wp
      field% aprc  (:,  :) = 0._wp
      field% aprs  (:,  :) = 0._wp

      field%  evap_ac(:,  :) = 0._wp
      field% lhflx_ac(:,  :) = 0._wp
      field% shflx_ac(:,  :) = 0._wp

      field% u_stress_ac(:,  :) = 0._wp
      field% v_stress_ac(:,  :) = 0._wp

      NULLIFY( field )
    ENDDO !jg

  END SUBROUTINE additional_restart_init

END MODULE mo_echam_phy_init

