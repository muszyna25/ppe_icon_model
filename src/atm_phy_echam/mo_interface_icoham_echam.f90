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
  USE mo_exception             ,ONLY: finish !, message
  USE mo_impl_constants        ,ONLY: min_rlcell_int
  USE mo_impl_constants_grf    ,ONLY: grf_bdywidth_e, grf_bdywidth_c

  USE mo_datetime              ,ONLY: t_datetime
  USE mo_model_domain          ,ONLY: t_patch
  USE mo_loopindices           ,ONLY: get_indices_c, get_indices_e
  USE mo_intp_data_strc        ,ONLY: t_int_state
  USE mo_intp_rbf              ,ONLY: rbf_vec_interpol_cell
  USE mo_sync                  ,ONLY: SYNC_C, SYNC_E, sync_patch_array, sync_patch_array_mult

  USE mo_timer                 ,ONLY: timer_start, timer_stop,        &
    &                                 timer_dyn2phy, timer_phy2dyn,   &
    &                                 timer_echam_phy, timer_coupling

  USE mo_parallel_config       ,ONLY: nproma, p_test_run
  USE mo_run_config            ,ONLY: nlev, ntracer, ltimer
  USE mo_echam_phy_config      ,ONLY: echam_phy_config
  USE mo_coupling_config       ,ONLY: is_coupled_run

  USE mo_icoham_dyn_types      ,ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_echam_phy_memory      ,ONLY: prm_field, prm_tend
  USE mo_eta_coord_diag        ,ONLY: half_level_pressure, full_level_pressure
  USE mo_echam_phy_bcs         ,ONLY: echam_phy_bcs_global
  USE mo_echam_phy_main        ,ONLY: echam_phy_main
  USE mo_interface_echam_ocean ,ONLY: interface_echam_ocean

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_icoham_echam

  CHARACTER(len=*), PARAMETER :: version = '$Id$'
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

  SUBROUTINE interface_icoham_echam( datetime,             &! in
    &                                pdtime, psteplen, jg, &! in
    &                                p_patch,              &! in
    &                                p_int_state,          &! in
    &                                dyn_prog_old,         &! in
    &                                dyn_diag_old,         &! in
    &                                dyn_prog_new,         &! in
    &                                dyn_tend             ) ! inout

    ! Arguments

    TYPE(t_datetime),      INTENT(IN)    :: datetime
    REAL(wp),              INTENT(IN)    :: pdtime        !< time step
    REAL(wp),              INTENT(IN)    :: psteplen      !< 2*time step in case of leapfrog
    INTEGER,               INTENT(IN)    :: jg            !< grid level/domain index
    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
    TYPE(t_int_state),TARGET,INTENT(IN)  :: p_int_state

    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: dyn_prog_old
    TYPE(t_hydro_atm_diag),INTENT(IN)    :: dyn_diag_old
    TYPE(t_hydro_atm_prog),INTENT(IN)    :: dyn_prog_new

    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: dyn_tend

    ! Local variables

    LOGICAL  :: any_uv_tend, ltrig_rad
    REAL(wp) :: ztime_radtran  !< time instance (in radian) at which radiative transfer is computed
    REAL(wp) :: zvn1, zvn2
    REAL(wp), POINTER :: zdudt(:,:,:), zdvdt(:,:,:)

    INTEGER :: jb,jbs   !< block index and its staring value
    INTEGER :: jcs,jce  !< start/end column index within each block
    INTEGER :: jc, jk   !< column index, vertical level index
    INTEGER :: jcn,jbn  !< column and block indices of a neighbour cell

    INTEGER:: i_nchdom  !< number of child patches
    INTEGER:: rl_start, rl_end, i_startblk, i_endblk

    INTEGER:: return_status
    CHARACTER(*), PARAMETER :: method_name = "interface_icoham_echam"

    !-------------------------------------------------------------------------
    IF (ltimer) CALL timer_start(timer_dyn2phy)

    ! Inquire current grid level and the total number of grid cells
    i_nchdom = MAX(1,p_patch%n_childdom)
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

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
    CALL sync_patch_array( SYNC_E, p_patch, dyn_prog_old%vn )

    CALL rbf_vec_interpol_cell( dyn_prog_old%vn,      &! in
      &                         p_patch, p_int_state, &! in
      &                         prm_field(jg)%u,      &! out
      &                         prm_field(jg)%v,      &! out
      &   opt_rlstart=rl_start, opt_rlend=rl_end     ) ! in    

!$OMP PARALLEL WORKSHARE
    prm_field(jg)%         q(:,:,:,:) = dyn_prog_old%     tracer(:,:,:,:)
    prm_field(jg)%      temp(:,:,:)   = dyn_prog_old%       temp(:,:,:)

    prm_field(jg)% presi_old(:,:,:)   = dyn_diag_old%    pres_ic(:,:,:)
    prm_field(jg)% presm_old(:,:,:)   = dyn_diag_old%    pres_mc(:,:,:)

    prm_field(jg)%        qx(:,:,:)   = dyn_diag_old%         qx(:,:,:)
    prm_field(jg)%        tv(:,:,:)   = dyn_diag_old%      tempv(:,:,:)
    prm_field(jg)%      geom(:,:,:)   = dyn_diag_old%     geo_mc(:,:,:)
    prm_field(jg)%      geoi(:,:,:)   = dyn_diag_old%     geo_ic(:,:,:)
    prm_field(jg)%     omega(:,:,:)   = dyn_diag_old%   wpres_mc(:,:,:)
    prm_field(jg)%       vor(:,:,:)   = dyn_diag_old% rel_vort_c(:,:,:)
!$OMP END PARALLEL WORKSHARE

    !---------------------------------
    ! Additional diagnostic variables

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c( p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

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
    CALL sync_patch_array( SYNC_E, p_patch, dyn_tend%vn )

    CALL rbf_vec_interpol_cell( dyn_tend%vn,          &! in
      &                         p_patch, p_int_state, &! in
      &                         prm_tend(jg)%u,       &! out
      &                         prm_tend(jg)%v,       &! out
      &   opt_rlstart=rl_start, opt_rlend=rl_end     ) ! in

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c( p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      prm_tend(jg)% temp(jcs:jce,:,jb)   = dyn_tend%   temp(jcs:jce,:,jb)
      prm_tend(jg)%    q(jcs:jce,:,jb,:) = dyn_tend% tracer(jcs:jce,:,jb,:)
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (ltimer)  THEN
      CALL timer_stop (timer_dyn2phy)
      CALL timer_start(timer_echam_phy)
    END IF

    !-------------------------------------------------------------------------
    ! Prepare some global parameters or parameter arrays
    !-------------------------------------------------------------------------
    !
    CALL echam_phy_bcs_global( datetime     ,&! in
      &                        jg           ,&! in
      &                        p_patch      ,&! in
      &                        pdtime       ,&! in
      &                        ltrig_rad    ,&! out
      &                        ztime_radtran) ! out

    !    WRITE(0,*)'radiation=',ltrig_rad, dt_rad
    !    WRITE(0,*)' vor PYHSC rad fluxes sw sfc',  MAXVAL(prm_field(jg)% swflxsfc_avg(:,:))
    !    WRITE(0,*)' vor PYHSC rad fluxes lw sfc', MINVAL(prm_field(jg)% lwflxsfc_avg(:,:))

    !---------------------------------------------------------------------------------
    ! For each block, call "echam_phy_main" to compute various parameterised processes
    !---------------------------------------------------------------------------------
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce),  ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk,i_endblk

      CALL get_indices_c(p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

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
        &                  ltrig_rad    ,&! in
        &                  ztime_radtran) ! in

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !    WRITE(0,*)' nach PYHSC rad fluxes sw sfc', MAXVAL( prm_field(jg)% swflxsfc_avg(:,:))
    !    WRITE(0,*)' nach PYHSC rad fluxes lw sfc', MINVAL( prm_field(jg)% lwflxsfc_avg(:,:))

    IF (ltimer)  THEN
      CALL timer_stop (timer_echam_phy)
    END IF

    !-------------------------------------------------------------------------
    ! If running in atm-oce coupled mode, exchange information 
    !-------------------------------------------------------------------------
    ! 
    IF ( is_coupled_run() ) THEN
      IF (ltimer) CALL timer_start(timer_coupling)

      CALL interface_echam_ocean( jg, p_patch )

      IF (ltimer) CALL timer_stop(timer_coupling)
    END IF

    !-------------------------------------------------------------------------
    ! Physics to dynamics: remap tendencies to the dynamics grid
    !-------------------------------------------------------------------------
    ! Currently this includes a simple copying of the temperature and tracer
    ! tendencies to the dynamics data structure, and convert the u- and v-wind
    ! tendencies to the normal wind tendency.
    ! Once a physics grid of different resolution is intruduced,
    ! conservative re-mapping will be called here.

    IF (ltimer)  THEN
      CALL timer_start(timer_phy2dyn)
    END IF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c( p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      dyn_tend%   temp(jcs:jce,:,jb)   = prm_tend(jg)% temp(jcs:jce,:,jb)
      dyn_tend% tracer(jcs:jce,:,jb,:) = prm_tend(jg)%    q(jcs:jce,:,jb,:)
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL sync_patch_array( SYNC_C, p_patch, dyn_tend%temp )
    CALL sync_patch_array_mult(SYNC_C, p_patch, ntracer, f4din=dyn_tend% tracer)

    any_uv_tend = echam_phy_config%lconv     .OR. &
      &           echam_phy_config%lvdiff    .OR. &
      &           echam_phy_config%lgw_hines .OR. &
      &           echam_phy_config%lssodrag

    IF (any_uv_tend) THEN

      ALLOCATE(zdudt(nproma,nlev,p_patch%nblks_c), &
        &      zdvdt(nproma,nlev,p_patch%nblks_c), &
        &      stat=return_status)
      IF (return_status > 0) THEN
        CALL finish (method_name, 'ALLOCATE(zdudt,zdvdt)')
      END IF
      IF (p_test_run) THEN
        zdudt(:,:,:) = 0.0_wp
        zdvdt(:,:,:) = 0.0_wp
      END IF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
        zdudt(jcs:jce,:,jb) =   prm_tend(jg)% u_phy(jcs:jce,:,jb)
        zdvdt(jcs:jce,:,jb) =   prm_tend(jg)% v_phy(jcs:jce,:,jb)
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! Now derive the physics-induced normal wind tendency, and add it to the
      ! total tendency.
      CALL sync_patch_array_mult(SYNC_C, p_patch, 2, zdudt, zdvdt)

      jbs   = p_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jcs,jce,jcn,jbn,zvn1,zvn2) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,p_patch%nblks_e
        CALL get_indices_e(p_patch, jb,jbs,p_patch%nblks_e, jcs,jce, grf_bdywidth_e+1)
        DO jk = 1,nlev
          DO jc = jcs,jce

            jcn  =   p_patch%edges%cell_idx(jc,jb,1)
            jbn  =   p_patch%edges%cell_blk(jc,jb,1)
            zvn1 =   zdudt(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(jc,jb,1)%v1 &
              &    + zdvdt(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(jc,jb,1)%v2

            jcn  =   p_patch%edges%cell_idx(jc,jb,2)
            jbn  =   p_patch%edges%cell_blk(jc,jb,2)
            zvn2 =   zdudt(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(jc,jb,2)%v1 &
              &    + zdvdt(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(jc,jb,2)%v2

            dyn_tend%vn(jc,jk,jb) =   dyn_tend%vn(jc,jk,jb)             &
              &                     + p_int_state%c_lin_e(jc,1,jb)*zvn1 &
              &                     + p_int_state%c_lin_e(jc,2,jb)*zvn2
          END DO ! jc
        END DO ! jk
      END DO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      DEALLOCATE(zdudt, zdvdt)

    END IF !any_uv_tend

    IF (ltimer) CALL timer_stop(timer_phy2dyn)

    !--------------------------
  END SUBROUTINE interface_icoham_echam

END MODULE mo_interface_icoham_echam
