!>
!! @brief Interface between ICOHAM dynamics+transport and ECHAM physics
!!
!! @author Kristina Froehlich (DWD)
!! @author Marco Giorgetta (MPI-M)
!! @author Hui Wan (MPI-M)
!!
!! @par Revision History
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

MODULE mo_interface_icoham_echam

  USE mo_kind                  ,ONLY: wp
  USE mo_exception             ,ONLY: finish !, message, message_text, print_value

  USE mo_impl_constants        ,ONLY: min_rlcell_int
  USE mo_impl_constants_grf    ,ONLY: grf_bdywidth_e, grf_bdywidth_c

  USE mo_coupling_config       ,ONLY: is_coupled_run
  USE mo_parallel_config       ,ONLY: nproma
  USE mo_run_config            ,ONLY: nlev, ntracer, ltimer
  USE mo_echam_phy_config      ,ONLY: echam_phy_config

  USE mo_model_domain          ,ONLY: t_patch
  USE mo_intp_data_strc        ,ONLY: t_int_state
  USE mo_intp_rbf              ,ONLY: rbf_vec_interpol_cell

  USE mo_loopindices           ,ONLY: get_indices_c, get_indices_e
  USE mo_sync                  ,ONLY: sync_c, sync_e, sync_patch_array, sync_patch_array_mult

  USE mo_icoham_dyn_types      ,ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_eta_coord_diag        ,ONLY: half_level_pressure, full_level_pressure

  USE mo_datetime              ,ONLY: t_datetime
  USE mo_echam_phy_memory      ,ONLY: prm_field, prm_tend
  USE mo_echam_phy_bcs         ,ONLY: echam_phy_bcs_global
  USE mo_echam_phy_main        ,ONLY: echam_phy_main
  USE mo_interface_echam_ocean ,ONLY: interface_echam_ocean

  USE mo_timer                 ,ONLY: timer_start, timer_stop,        &
    &                                 timer_dyn2phy, timer_phy2dyn,   &
    &                                 timer_echam_phy, timer_coupling

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: interface_icoham_echam

  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_interface_icoham_echam'

CONTAINS
  !>
  !! SUBROUTINE echam_physics -- the Interface between ICON dynamics and
  !! ECHAM physics
  !!
  !! This subroutine is called in the time loop of the ICOHAM model.
  !! It takes the following as input:
  !! <ol>
  !! <li> prognostic and diagnostic variables of the dynamical core;
  !! <li> tendency of the prognostic varibles induced by adiabatic dynamics;
  !! <li> time step;
  !! <li> information about the dynamics grid;
  !! <li> interplation coefficients.
  !! </ol>
  !!
  !! The output includes tendencies of the prognostic variables caused by
  !! the parameterisations.
  !!
  !! Note that each call of this subroutine deals with a single grid level
  !! rather than the entire grid tree.

  SUBROUTINE interface_icoham_echam( pdtime, psteplen ,& !in
    &                                datetime         ,& !in
    &                                patch            ,& !in
    &                                pt_int_state     ,& !in
    &                                dyn_prog_old     ,& !in
    &                                dyn_diag_old     ,& !in
    &                                dyn_prog_new     ,& !in
    &                                dyn_tend         )  !inout

    !
    !> Arguments:
    !
    REAL(wp)              , INTENT(in)            :: pdtime          !< time step
    REAL(wp)              , INTENT(in)            :: psteplen        !< 2*time step in case of leapfrog
    TYPE(t_datetime)      , INTENT(in)            :: datetime

    TYPE(t_patch)         , INTENT(in)   , TARGET :: patch           !< grid/patch info
    TYPE(t_int_state)     , INTENT(in)   , TARGET :: pt_int_state    !< interpolation state

    TYPE(t_hydro_atm_prog), INTENT(inout)         :: dyn_prog_old
    TYPE(t_hydro_atm_diag), INTENT(in)            :: dyn_diag_old
    TYPE(t_hydro_atm_prog), INTENT(in)            :: dyn_prog_new

    TYPE(t_hydro_atm_prog), INTENT(inout)         :: dyn_tend

    ! Local array bounds

    INTEGER  :: i_nchdom             !< number of child patches
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: rl_start, rl_end
    INTEGER  :: jg                   !< grid index
    INTEGER  :: jcs, jce             !< start and end cell indices
    INTEGER  :: jes, jee, je         !< start and end edge indices
    INTEGER  :: jk                   !< level in column index
    INTEGER  :: jb, jbs, jbe         !< row in block index, start and end indices
    INTEGER  :: jcn,jbn              !< jc and jb of neighbor cells sharing an edge je

    ! Local variables

    REAL(wp) :: zvn1, zvn2
    REAL(wp), POINTER :: zdudt(:,:,:), zdvdt(:,:,:)

    LOGICAL  :: any_uv_tend
    LOGICAL  :: ltrig_rad

    INTEGER  :: return_status

    ! Local parameters

    CHARACTER(*), PARAMETER :: method_name = "interface_icoham_echam"

    !-------------------------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_dyn2phy)

    ! Inquire current grid level and the total number of grid cells
    i_nchdom = MAX(1,patch%n_childdom)
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    jg    = patch%id

    !-------------------------------------------------------------------------
    ! Dynamics to physics: remap dynamics variables to physics grid
    !-------------------------------------------------------------------------
    ! Currently this includes
    !  - reconstructing of u- and v-wind and their tendencies at cell centers
    !  - copying scalar fields from the dynamics state and from the tendency
    !    state to the physics state and physics tendencies, repsectively.
    !  - computing pressure values at the "new" time step
    ! Once a physics grid of different resolution is intruduced, conservative
    ! re-mapping will be called here.


    ! LL The physics runs only on the owned cells
    !  but the following rbf_vec_interpol_cell may use the halos(?)
    !
    ! - prm_field(jg)%u
    ! - prm_field(jg)%v
    CALL sync_patch_array( SYNC_E, patch, dyn_prog_old%vn )

    CALL rbf_vec_interpol_cell( dyn_prog_old%vn       ,&! in
      &                         patch                 ,&! in
      &                         pt_int_state          ,&! in
      &                         prm_field(jg)%u       ,&! out
      &                         prm_field(jg)%v       ,&! out
      &                         opt_rlstart=rl_start  ,&! in
      &                         opt_rlend=rl_end      ) ! in

!$OMP PARALLEL WORKSHARE

    ! Fill the physics state variables, which are used by echam:
    !
    prm_field(jg)%      geom(:,:,:)   = dyn_diag_old%     geo_mc(:,:,:)
    !
    prm_field(jg)%       vor(:,:,:)   = dyn_diag_old% rel_vort_c(:,:,:)
    !
    prm_field(jg)%      temp(:,:,:)   = dyn_prog_old%       temp(:,:,:)
    prm_field(jg)%        tv(:,:,:)   = dyn_diag_old%      tempv(:,:,:)
    !
    prm_field(jg)% presm_old(:,:,:)   = dyn_diag_old%    pres_mc(:,:,:)
    !
    prm_field(jg)%      qtrc(:,:,:,:) = dyn_prog_old%     tracer(:,:,:,:)
    !
    ! cloud water+ice
    prm_field(jg)%        qx(:,:,:)   = dyn_diag_old%         qx(:,:,:)
    !
    ! vertical velocity in p-system
    prm_field(jg)%     omega(:,:,:)   = dyn_diag_old%   wpres_mc(:,:,:)
    !
    prm_field(jg)%      geoi(:,:,:)   = dyn_diag_old%     geo_ic(:,:,:)
    !
    prm_field(jg)% presi_old(:,:,:)   = dyn_diag_old%    pres_ic(:,:,:)

!$OMP END PARALLEL WORKSHARE

    !---------------------------------
    ! Additional diagnostic variables

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c( patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      ! Pressure at time step "new" (i.e., n+1)

      CALL half_level_pressure( dyn_prog_new%pres_sfc(:,jb),     nproma, jce, &! in
        &                       prm_field(jg)%presi_new(:,:,jb)               )! out

      CALL full_level_pressure( prm_field(jg)%presi_new(:,:,jb), nproma, jce, &! in
        &                       prm_field(jg)%presm_new(:,:,jb)               )! out
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !--------------------------------
    ! transfer tendencies


    ! LL The physics runs only on the owned cells
    !    but the following rbf_vec_interpol_cell may use the halos(?)
    CALL sync_patch_array( SYNC_E, patch, dyn_tend%vn )

    CALL rbf_vec_interpol_cell( dyn_tend%vn,          &! in
      &                         patch, pt_int_state,  &! in
      &                         prm_tend(jg)%u,       &! out
      &                         prm_tend(jg)%v,       &! out
      &   opt_rlstart=rl_start, opt_rlend=rl_end     ) ! in

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c( patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      prm_tend(jg)%    temp(jcs:jce,:,jb)   = dyn_tend%   temp(jcs:jce,:,jb)
      prm_tend(jg)%    qtrc(jcs:jce,:,jb,:) = dyn_tend% tracer(jcs:jce,:,jb,:)
      prm_tend(jg)%qtrc_dyn(jcs:jce,:,jb,:) = dyn_tend% tracer(jcs:jce,:,jb,:)
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !
    !=====================================================================================

    IF (ltimer) THEN
      CALL timer_stop (timer_dyn2phy)
      CALL timer_start(timer_echam_phy)
    END IF

    !=====================================================================================
    !
    ! (3) Prepare boundary conditions for ECHAM physics
    !
    CALL echam_phy_bcs_global( datetime     ,&! in
      &                        jg           ,&! in
      &                        patch        ,&! in
      &                        pdtime       ,&! in
      &                        ltrig_rad    ) ! out
    !
    !=====================================================================================

    !=====================================================================================
    !
    ! (4) Call echam physics and compute the total physics tendencies.
    !     This includes the atmospheric processes (proper ECHAM) and
    !     the land processes, which are vertically implicitly coupled
    !     to the parameterization of vertical turbulent fluxes.
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce),  ICON_OMP_GUIDED_SCHEDULE

    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      ! Like in ECHAM, the subroutine *echam_phy_main* has direct access to the memory
      ! buffers prm_field and prm_tend. In addition it can also directly access
      ! the grid/patch information on which the computations are performed.
      ! Thus the argument list contains only
      ! - jg: the grid index in the grid hierarchy
      ! - jb: the row index in the block
      ! - jcs and jce: start and end indices of columns in a row
      ! - nproma: the block length
      ! - a few other globally valid arguments

      CALL echam_phy_main( jg           ,&! in
        &                  jb           ,&! in
        &                  jcs          ,&! in
        &                  jce          ,&! in
        &                  nproma       ,&! in
        &                  datetime     ,&! in
        &                  pdtime       ,&! in
        &                  psteplen     ,&! in
        &                  ltrig_rad     )

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    !
    !=====================================================================================

    IF (ltimer) CALL timer_stop(timer_echam_phy)

    !=====================================================================================
    !
    ! (5) Couple to ocean surface if an ocean is present and this is a coupling time step.
    !
    !
    IF ( is_coupled_run() ) THEN
      IF (ltimer) CALL timer_start(timer_coupling)

      CALL interface_echam_ocean( jg, patch )

      IF (ltimer) CALL timer_stop(timer_coupling)
    END IF
    !
    !=====================================================================================

    IF (ltimer) CALL timer_start(timer_phy2dyn)

    !=====================================================================================
    !
    !     Copy  physics tandencies in temp. and tracers from the physics to the dynamics
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c( patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      dyn_tend%   temp(jcs:jce,:,jb)   = prm_tend(jg)% temp(jcs:jce,:,jb)
      dyn_tend% tracer(jcs:jce,:,jb,:) = prm_tend(jg)% qtrc(jcs:jce,:,jb,:)
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    !
    CALL sync_patch_array( SYNC_C, patch, dyn_tend%temp )
    CALL sync_patch_array_mult(SYNC_C, patch, ntracer, f4din=dyn_tend% tracer)
    !
    !=====================================================================================


    !=====================================================================================
    !
    ! (6) Convert physics tandencies in the wind components (u,v) to tendencies in
    !     normal wind vn.
    !
    !
    any_uv_tend = echam_phy_config%lconv     .OR. &
      &           echam_phy_config%lvdiff    .OR. &
      &           echam_phy_config%lgw_hines .OR. &
      &           echam_phy_config%lssodrag

    IF (any_uv_tend) THEN

      ALLOCATE(zdudt(nproma,nlev,patch%nblks_c), &
        &      zdvdt(nproma,nlev,patch%nblks_c), &
        &      stat=return_status)
      IF (return_status > 0) THEN
        CALL finish (method_name, 'ALLOCATE(zdudt,zdvdt)')
      END IF
      zdudt(:,:,:) = 0.0_wp
      zdvdt(:,:,:) = 0.0_wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
        zdudt(jcs:jce,:,jb) = prm_tend(jg)% u_phy(jcs:jce,:,jb)
        zdvdt(jcs:jce,:,jb) = prm_tend(jg)% v_phy(jcs:jce,:,jb)
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! Now derive the physics-induced normal wind tendency, and add it to the
      ! total tendency.
      CALL sync_patch_array_mult(SYNC_C, patch, 2, zdudt, zdvdt)

      jbs   = patch%edges%start_blk(grf_bdywidth_e+1,1)
      jbe   = patch%nblks_e

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,jes,jee,jcn,jbn,zvn1,zvn2) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,jbe
        CALL get_indices_e(patch, jb,jbs,jbe, jes,jee, grf_bdywidth_e+1)

        DO jk = 1,nlev

          DO je = jes,jee
            jcn  =   patch%edges%cell_idx(je,jb,1)
            jbn  =   patch%edges%cell_blk(je,jb,1)
            zvn1 =   zdudt(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,1)%v1 &
              &    + zdvdt(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,1)%v2

            jcn  =   patch%edges%cell_idx(je,jb,2)
            jbn  =   patch%edges%cell_blk(je,jb,2)
            zvn2 =   zdudt(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,2)%v1 &
              &    + zdvdt(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,2)%v2

            dyn_tend%vn(je,jk,jb)        =   dyn_tend%vn(je,jk,jb)              &
              &                            + pt_int_state%c_lin_e(je,1,jb)*zvn1 &
              &                            + pt_int_state%c_lin_e(je,2,jb)*zvn2

          END DO ! je
        END DO ! jk
      END DO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      DEALLOCATE(zdudt, zdvdt)

    END IF ! any_uv_tend
    !
    !=====================================================================================

    IF (ltimer) CALL timer_stop(timer_phy2dyn)

  END SUBROUTINE interface_icoham_echam
  !----------------------------------------------------------------------------

END MODULE mo_interface_icoham_echam
