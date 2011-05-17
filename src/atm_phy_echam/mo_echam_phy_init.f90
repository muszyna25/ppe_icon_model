!>
!! @brief Contains subroutines for initializing the ECHAM physics
!! package in ICOHAM.
!!
!! @author Hui Wan, MPI-M 
!!
!! @par Revision History
!! First version by Hui Wan, 2010-07-20
!!!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
  USE mo_exception,            ONLY: finish, message

  ! model configuration
  USE mo_dynamics_nml,       ONLY: nnow
  USE mo_run_nml,            ONLY: nlev, nlevp1, nvclev, nproma, &
    &                              iqv, iqt, ntracer, ltestcase
  USE mo_vertical_coord_table,ONLY: vct
  USE mo_echam_phy_nml,      ONLY: lrad, lvdiff, lconv, lcond
  USE mo_vertical_coord_table, ONLY: vct
  ! test cases
  USE mo_hydro_testcases,    ONLY: ctest_name, ape_sst_case
  USE mo_ape_params,         ONLY: ape_sst

  ! radiation
  USE mo_radiation_nml,        ONLY: ssi, tsi
  USE mo_srtm_config,          ONLY: setup_srtm, ssi_amip
  USE mo_lrtm_setup,           ONLY: lrtm_setup
  USE mo_newcld_optics,        ONLY: setup_newcld_optics

  ! vertical diffusion
  USE mo_echam_vdiff_params,   ONLY: init_vdiff_params, z0m_min
  USE mo_vdiff_solver,         ONLY: init_vdiff_solver

  ! cumulus convection
  USE mo_convect_tables,       ONLY: init_convect_tables
  USE mo_echam_conv_params,    ONLY: cuparam

  ! stratiform clouds and cloud cover
  USE mo_echam_cloud_params,   ONLY: init_cloud_tables, sucloud, cvarmin

  ! air-sea-land interface
  USE mo_icoham_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd, &
    &                                init_sfc_indices

  ! domain and indices
  USE mo_model_domain,         ONLY: t_patch
  USE mo_loopindices,          ONLY: get_indices_c

  ! atmospheric state
  USE mo_icoham_dyn_types,     ONLY: t_hydro_atm
  USE mo_eta_coord_diag,       ONLY: half_level_pressure, full_level_pressure
  USE mo_echam_phy_memory,     ONLY: construct_echam_phy_state,    &
    &                                prm_field, t_echam_phy_field, &
    &                                prm_tend,  t_echam_phy_tend


  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: init_echam_phy

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !! The top-level routine for the initialization of ECHAM6 physics.
  !! It calls a series of subroutines to initialize tunable parameters,
  !! lookup tables, and the physics state vectors "prm_field" and "prm_tend".
  !!
  !! @par Revision History
  !! Initial version by Hui Wan, MPI-M (2010-07)
  !!
  SUBROUTINE init_echam_phy( p_patch, p_hydro_state )

    TYPE(t_patch)      ,INTENT(IN) :: p_patch(:)
    TYPE(t_hydro_atm),INTENT(IN) :: p_hydro_state(:)
    INTEGER :: khydromet, ktrac

    !-------------------------------------------------------------------
    ! Initialize parameters and lookup tables
    !-------------------------------------------------------------------
    ! For radiation

    IF (lrad) THEN
      ssi(:) = ssi_amip(:)
      tsi    = SUM(ssi(:))
      CALL setup_srtm
      CALL lrtm_setup
      CALL setup_newcld_optics
    END IF

    ! For surface processes: 
    ! nsfc_type, iwtr, etc. are set in this subroutine. 
    ! See mo_icoham_sfc_indicies.f90 for further details.

    CALL init_sfc_indices( ltestcase, ctest_name )

    ! For turbulent mixing:
    ! Allocate memory for the tri-diagonal solver needed by the implicit
    ! time stepping scheme; Compute time-independent parameters.

    IF (lvdiff) THEN
      ! Currently the tracer indices are sorted such that we count
      ! the water substances first, and then other species like 
      ! aerosols and their precursors. "ntracer" is the total number 
      ! of tracers (including water substances) handled in the model;
      ! "iqt" is the starting index for non-water species.
      ! Before more sophisticated meta-data structure becomes available, 
      ! it is assumed here that all tracers are subject to turbulent mixing.

      khydromet = iqt - 2        ! # of hydrometeors
      ktrac = ntracer - iqt + 1  ! # of non-water species 

     !CALL init_vdiff_solver( khydromet, ktrac, nproma, nlev, nsfc_type )
      CALL init_vdiff_solver( khydromet, ktrac, nlev )
      CALL init_vdiff_params( nlev, nlevp1, nvclev, vct )
    ENDIF

    ! Lookup tables for saturation vapour pressure

    IF (lconv.OR.lcond.OR.lvdiff)  CALL init_convect_tables

    ! For cumulus convection

    IF (lconv) CALL cuparam

    ! For large scale condensation

    IF (lcond) THEN
      CALL init_cloud_tables
      CALL sucloud( nlev, vct        &
!!$        &         , lmidatm=.FALSE.  &
        &         , lcouple=.FALSE.  &
        &         , lipcc=.FALSE.    &
!!$        &         , lham=.FALSE.     &
        &         )
    END IF

    !-------------------------------------------------------------------
    ! Allocate memory for the state vectors "prm_field" and "prm_tend";
    ! Give intial values to some of their components.
    !-------------------------------------------------------------------
    CALL construct_echam_phy_state( p_patch )
    CALL init_phy_memory( p_patch, p_hydro_state )

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
  SUBROUTINE init_phy_memory( p_patch, p_hydro_state )

    TYPE(t_patch)      ,INTENT(IN) :: p_patch(:)
    TYPE(t_hydro_atm),INTENT(IN) :: p_hydro_state(:)

    ! local variables and pointers

    INTEGER  :: ndomain, nblks_c, jg, jb, jbs, jc, jcs, jce, jk
    REAL(wp) :: zprat, zn1, zn2, zcdnc, zlat
    LOGICAL  :: lland, lglac
    INTEGER  :: nexp

    TYPE(t_echam_phy_field),POINTER :: field => NULL()
    TYPE(t_echam_phy_tend) ,POINTER :: tend  => NULL()
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
      tend  => prm_tend (jg)

      !----------------------------------------
      ! Loop over all blocks in domain jg
      !----------------------------------------
      nblks_c = p_patch(jg)%nblks_int_c
      jbs     = p_patch(jg)%cells%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jcs,jce,jk)
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
             nexp=2
             zprat=(MIN(8._wp,80000._wp/field%presm_old(jc,jk,jb)))**nexp

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

      field% rtype (:,  :) = 0._wp
      field% rintop(:,  :) = 0._wp

       tend% x_dtr(:,:,:) = 0._wp  !"xtec" in ECHAM
!$OMP END PARALLEL WORKSHARE

      IF (lvdiff) THEN
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

  END SUBROUTINE init_phy_memory
  !-------------

END MODULE mo_echam_phy_init

