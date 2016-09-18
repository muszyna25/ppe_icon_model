!>
!!
!!----------------------------------------------------------------------
!! This module provides the interface between dynamics+transport and
!! echam physics and controls the coupling procedure dependent on the
!! switch echam_phy_config%idcphycpl.
!!
!! @author Marco Giorgetta (MPI-M)
!!
!!----------------------------------------------------------------------
!!
!! @brief Interface between ICONAM dynamics+transport and ECHAM physics
!!
!! The coupling mechanism is controlled by echam_phy_config%idcphycpl:
!!
!! idcphycpl = 1: The whole physics is treated as "fast" phyiscs.
!!                The physics tendencies are computed from the
!!                provisional state reached after dynamics&transport
!!                and the full physics tendencies are then used to
!!                update and reach the final new state
!!
!! idcphycpl = 2: The whole physics is treated as "slow" phyiscs.
!!                The state after dynamics+transport is the final
!!                state for which the full phyiscs tandencies are
!!                computed. These will be used in the following
!!                timestep as forcing for the dynamics and for
!!                updating tracers after the transport.
!!
!!
!! @par Revision History
!!  first implementation by Marco Giorgetta, MPI-M (2014-03-27)
!!
!! @par Copyright
!! 2002-2014 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!      violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!      copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!      an according license agreement with DWD and MPI-M.
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

MODULE mo_interface_iconam_echam

  USE mo_kind                  ,ONLY: wp
  USE mo_exception             ,ONLY: finish

  USE mo_impl_constants        ,ONLY: min_rlcell_int, grf_bdywidth_c, grf_bdywidth_e

  USE mo_coupling_config       ,ONLY: is_coupled_run
  USE mo_parallel_config       ,ONLY: nproma
  USE mo_run_config            ,ONLY: nlev, ntracer, iqv, iqc, iqi
  USE mo_nonhydrostatic_config ,ONLY: lhdiff_rcf
  USE mo_diffusion_config      ,ONLY: diffusion_config
  USE mo_echam_phy_config      ,ONLY: echam_phy_config

  USE mo_model_domain          ,ONLY: t_patch
  USE mo_intp_data_strc        ,ONLY: t_int_state
  USE mo_intp_rbf              ,ONLY: rbf_vec_interpol_cell

  USE mo_loopindices           ,ONLY: get_indices_c, get_indices_e
  USE mo_sync                  ,ONLY: sync_c, sync_e, sync_patch_array, sync_patch_array_mult

  USE mo_nonhydro_types        ,ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nh_diagnose_pres_temp ,ONLY: diagnose_pres_temp
  USE mo_physical_constants    ,ONLY: rd, p0ref, rd_o_cpd, vtmpc1, grav

  USE mo_datetime              ,ONLY: t_datetime
  USE mo_echam_phy_memory      ,ONLY: prm_field, prm_tend
  USE mo_echam_phy_bcs         ,ONLY: echam_phy_bcs_global
  USE mo_echam_phy_main        ,ONLY: echam_phy_main
  USE mo_interface_echam_ocean ,ONLY: interface_echam_ocean
  
#ifndef __NO_JSBACH__
  USE mo_jsb_interface         ,ONLY: jsbach_start_timestep, jsbach_finish_timestep
#endif
  
  USE mo_timer                 ,ONLY: ltimer, timer_start, timer_stop,           &
    &                                 timer_dyn2phy, timer_d2p_prep, timer_d2p_sync, timer_d2p_couple, &
    &                                 timer_echam_bcs, timer_echam_phy, timer_coupling,                &
    &                                 timer_phy2dyn, timer_p2d_prep, timer_p2d_sync, timer_p2d_couple

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: interface_iconam_echam

  CHARACTER(len=*), PARAMETER :: module_name = 'mo_interface_iconam_echam'

CONTAINS
  !
  !-----------------------------------------------------------------------
  !
  !  This subroutine works as interface between dynamics+transport and
  !  echam physics.
  !
  !  Marco Giorgetta, MPI-M, 2014
  !
  SUBROUTINE interface_iconam_echam( dtadv_loc        ,& !in
    &                                datetime         ,& !in
    &                                patch            ,& !in
    &                                pt_int_state     ,& !in
    &                                p_metrics        ,& !in
    &                                pt_prog_new      ,& !inout
    &                                pt_prog_new_rcf  ,& !inout
    &                                pt_diag          )  !inout

    !
    !> Arguments:
    !
    REAL(wp)              , INTENT(in)            :: dtadv_loc       !< advective time step
    TYPE(t_datetime)      , INTENT(in)            :: datetime

    TYPE(t_patch)         , INTENT(in)   , TARGET :: patch           !< grid/patch info
    TYPE(t_int_state)     , INTENT(in)   , TARGET :: pt_int_state    !< interpolation state
    TYPE(t_nh_metrics)    , INTENT(in)            :: p_metrics

    TYPE(t_nh_diag)       , INTENT(inout), TARGET :: pt_diag         !< diagnostic variables
    TYPE(t_nh_prog)       , INTENT(inout), TARGET :: pt_prog_new     !< progn. vars after dynamics  for wind, temp. rho, ...
    TYPE(t_nh_prog)       , INTENT(inout), TARGET :: pt_prog_new_rcf !< progn. vars after advection for tracers

    ! Local array bounds

    INTEGER  :: i_nchdom             !< number of child patches
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: rl_start, rl_end
    INTEGER  :: jg                   !< grid index
    INTEGER  :: jc, jcs, jce         !< cell in row index, start and end indices
    INTEGER  :: je, jes, jee         !< edge in row index, start and end indices
    INTEGER  :: jk                   !< level in column index
    INTEGER  :: jb, jbs, jbe         !< row in block index, start and end indices
    INTEGER  :: jcn,jbn              !< jc and jb of neighbor cells sharing an edge je
    INTEGER  :: jt                   !< tracer index

    ! Local variables

    REAL(wp) :: z_exner              !< to save provisional new exner
    REAL(wp) :: z_qsum               !< summand of virtual increment
!!$    REAL(wp) :: z_ddt_qsum           !< summand of virtual increment

    REAL(wp) :: zvn1, zvn2
    REAL(wp), POINTER :: zdudt(:,:,:), zdvdt(:,:,:)

    LOGICAL  :: ltrig_rad
    TYPE(t_datetime)   :: datetime_radtran !< date and time for radiative transfer calculation

    INTEGER  :: return_status

    ! Local parameters

    CHARACTER(*), PARAMETER :: method_name = "interface_iconam_echam"

    !-------------------------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_dyn2phy)

    ! Inquire current grid level and the total number of grid cells
    i_nchdom = MAX(1,patch%n_childdom)
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    jg    = patch%id

    !=====================================================================================
    !
    ! (1) Complete prognostic and diagnostic state as needed for the computation
    !     of the phyiscs tendendies.
    !
    !     idcphycpl = 1 : fast physics coupling, dynamics and physics update sequentially
    !     idcphycpl = 2 : slow physics coupling, dynamics uses physics forcing for updating
    !
    !
    ! Update prognostic variables
    !
    IF (ltimer) CALL timer_start(timer_d2p_prep)

    SELECT CASE (echam_phy_config%idcphycpl)
      !
    CASE (1) ! idcphycpl
      ! In this case all ECHAM physics is treated as "fast" physics.
      ! Nothing needs to be updated here.
      !
    CASE (2) ! idcphycpl
      ! In this case all ECHAM physics is treated as "slow" physics.
      ! The provisional "new" tracer state, resulting from the advection
      ! step, still needs to be updated with the physics tracer tendencies
      ! computed by the physics called at the end of the previous time step
      ! for the then final "new" state that is used in this time step as
      ! the "now" state for the dynamics and advection.
      ! (The physics tendencies for the other variables have been used
      ! already in the dynamical core.)
      !
!$OMP PARALLEL
!$OMP DO PRIVATE(jt,jb,jk,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
      DO jt = 1,ntracer
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
          DO jk = 1,nlev
            DO jc = jcs, jce
              !
              pt_prog_new_rcf%tracer(jc,jk,jb,jt) =   pt_prog_new_rcf%tracer(jc,jk,jb,jt)             &
                &                                   + prm_tend(jg)% qtrc_phy(jc,jk,jb,jt) * dtadv_loc
              !
            END DO
          END DO
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
      !
    END SELECT ! idcphycpl

    ! Diagnostics
    !
    ! - pt_diag%tempv
    ! - pt_diag%temp
    ! - pt_diag%pres_sfc   surface pressure filtered to remove sound waves, see diagnose_pres_temp
    ! - pt_diag%pres_ifc   hydrostatic pressure at layer interface
    ! - pt_diag%pres       hydrostatic pressure at layer midpoint = SQRT(upper pres_ifc * lower pres_ifc)
    ! - pt_diag%dpres_mc   pressure thickness of layer
    !
    CALL diagnose_pres_temp( p_metrics                ,&
      &                      pt_prog_new              ,&
      &                      pt_prog_new_rcf          ,&
      &                      pt_diag                  ,&
      &                      patch                    ,&
      &                      opt_calc_temp=.TRUE.     ,&
      &                      opt_calc_pres=.TRUE.     ,&
      &                      opt_rlend=min_rlcell_int )

    IF (ltimer) CALL timer_stop(timer_d2p_prep)
    !
    ! - pt_diag%u
    ! - pt_diag%v
    IF (ltimer) CALL timer_start(timer_d2p_sync)
    CALL sync_patch_array( SYNC_E, patch, pt_prog_new%vn )
    IF (ltimer) CALL timer_stop(timer_d2p_sync)
    !
    IF (ltimer) CALL timer_start(timer_d2p_prep)

    CALL rbf_vec_interpol_cell( pt_prog_new%vn       ,&! in
      &                         patch                ,&! in
      &                         pt_int_state         ,&! in
      &                         pt_diag%u            ,&! out
      &                         pt_diag%v            ,&! out
      &                         opt_rlstart=rl_start ,&! in
      &                         opt_rlend  =rl_end   ) ! in

    IF (ltimer) CALL timer_stop(timer_d2p_prep)
    !
    ! Now the new prognostic and diagnostic state variables (pt_prog_new, pt_prog_new_rcf,
    ! pt_diag) of the dynamical core are ready to be used in the phyiscs.
    !
    ! idcphycpl = 1: This is the provisional "new" state that still needs to be
    !                updated by the physics to obtain the final "new" state.
    !
    ! idcphycpl = 2: This is the final "new" state, for which the phyiscs
    !                forcing is computed that will be applied in the next
    !                dynamical step(s) (pt_diag%ddt_exner_phy) and advection
    !                step (pt_diag%ddt_vn_phy, see above).
    !
    !=====================================================================================


    !=====================================================================================
    !
    ! (2) Copy the new prognostic state and the related diagnostics from the
    !     dynamics state variables to the phyiscs state variables
    !

    ! Loop over cells
    IF (ltimer) CALL timer_start(timer_d2p_couple)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      DO jk = 1,nlev
        DO jc = jcs, jce

          ! Fill the physics state variables, which are used by echam:
          !
          prm_field(jg)%        zf(jc,jk,jb)     = p_metrics%        z_mc(jc,jk,jb)
          prm_field(jg)%        dz(jc,jk,jb)     = p_metrics% ddqz_z_full(jc,jk,jb)
          prm_field(jg)%      geom(jc,jk,jb)     = p_metrics%  geopot_agl(jc,jk,jb)
          !
          prm_field(jg)%         u(jc,jk,jb)     = pt_diag%    u(jc,jk,jb)
          prm_field(jg)%         v(jc,jk,jb)     = pt_diag%    v(jc,jk,jb)
          prm_field(jg)%       vor(jc,jk,jb)     = pt_diag%  vor(jc,jk,jb)
          !
          prm_field(jg)%      temp(jc,jk,jb)     = pt_diag% temp(jc,jk,jb)
          prm_field(jg)%        tv(jc,jk,jb)     = pt_diag% tempv(jc,jk,jb)
          !
          prm_field(jg)% presm_old(jc,jk,jb)     = pt_diag%  pres(jc,jk,jb)
          prm_field(jg)% presm_new(jc,jk,jb)     = pt_diag%  pres(jc,jk,jb)
          !
          ! Air mass
          prm_field(jg)%       mair(jc,jk,jb)    = pt_prog_new %         rho(jc,jk,jb) &
            &                                     *prm_field(jg)%         dz(jc,jk,jb)
          !
          ! H2O mass (vap+liq+ice)
          prm_field(jg)%      mh2o(jc,jk,jb)     = ( pt_prog_new_rcf% tracer(jc,jk,jb,iqv)  &
            &                                       +pt_prog_new_rcf% tracer(jc,jk,jb,iqc)  &
            &                                       +pt_prog_new_rcf% tracer(jc,jk,jb,iqi)) &
            &                                      *prm_field(jg)%      mair(jc,jk,jb)
          !
          ! Dry air mass
          prm_field(jg)%      mdry(jc,jk,jb)     = prm_field(jg)%       mair(jc,jk,jb) &
            &                                     -prm_field(jg)%       mh2o(jc,jk,jb)
          !
          ! cloud water+ice
          IF (echam_phy_config%ldrymoist) THEN
            prm_field(jg)%        qx(jc,jk,jb)     = ( pt_prog_new_rcf% tracer(jc,jk,jb,iqc)  &
              &                                       +pt_prog_new_rcf% tracer(jc,jk,jb,iqi)) &
              &                                      *prm_field(jg)%      mair(jc,jk,jb)      &
              &                                      /prm_field(jg)%      mdry(jc,jk,jb)
          ELSE
            prm_field(jg)%        qx(jc,jk,jb)     = ( pt_prog_new_rcf% tracer(jc,jk,jb,iqc)  &
              &                                       +pt_prog_new_rcf% tracer(jc,jk,jb,iqi))
          END IF
          !
          ! vertical velocity in p-system
          prm_field(jg)%     omega(jc,jk,jb)     = -0.5_wp                                               &
            &                                      * (pt_prog_new%w(jc,jk,jb)+pt_prog_new%w(jc,jk+1,jb)) &
            &                                      * pt_prog_new%rho(jc,jk,jb) * grav
          !
          ! Tendencies passed to the ECHAM physics for internal upating are set to 0
          ! because the state passed to physics is already updated with tendencies
          ! due to dynamics and transport.
          prm_tend(jg)%          u(jc,jk,jb)     = 0.0_wp
          prm_tend(jg)%          v(jc,jk,jb)     = 0.0_wp
          !
          prm_tend(jg)%       temp(jc,jk,jb)     = 0.0_wp
          !

        END DO
      END DO

      DO jk = 1,nlev+1
        DO jc = jcs, jce

          prm_field(jg)%            zh(jc,jk,jb) = p_metrics%          z_ifc(jc,jk,jb)
          prm_field(jg)%          geoi(jc,jk,jb) = p_metrics% geopot_agl_ifc(jc,jk,jb)
          !
          prm_field(jg)%     presi_old(jc,jk,jb) = pt_diag% pres_ifc(jc,jk,jb)
          prm_field(jg)%     presi_new(jc,jk,jb) = pt_diag% pres_ifc(jc,jk,jb)

        END DO
      END DO

    END DO ! jb
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP DO PRIVATE(jt,jb,jk,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jt = 1,ntracer
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
        DO jk = 1,nlev
          DO jc = jcs, jce

            ! Tracer mass
            prm_field(jg)%      mtrc(jc,jk,jb,jt)  = pt_prog_new_rcf%  tracer(jc,jk,jb,jt) &
               &                                    *prm_field(jg)%      mair(jc,jk,jb)
            !
            ! Tracer mass fraction
            IF (echam_phy_config%ldrymoist) THEN
              prm_field(jg)%      qtrc(jc,jk,jb,jt)  = pt_prog_new_rcf% tracer(jc,jk,jb,jt) &
                &                                     *prm_field(jg)%     mair(jc,jk,jb)    &
                &                                     /prm_field(jg)%     mdry(jc,jk,jb)
            ELSE
              prm_field(jg)%      qtrc(jc,jk,jb,jt)  = pt_prog_new_rcf% tracer(jc,jk,jb,jt)
            END IF
            !
            ! Tendencies passed to the ECHAM physics for internal upating are set to 0
            ! because the state passed to physics is already updated with tendencies
            ! due to dynamics and transport.
            prm_tend(jg)%       qtrc(jc,jk,jb,jt)  = 0.0_wp
            !
            ! Advective tendencies, already accounted for, but needed
            ! for diagnostic purposes in the convection scheme
            IF (echam_phy_config%ldrymoist) THEN
              prm_tend(jg)%   qtrc_dyn(jc,jk,jb,jt)  = pt_diag% ddt_tracer_adv(jc,jk,jb,jt) &
                &                                     *prm_field(jg)%     mair(jc,jk,jb)    &
                &                                     /prm_field(jg)%     mdry(jc,jk,jb)
            ELSE
              prm_tend(jg)%   qtrc_dyn(jc,jk,jb,jt)  = pt_diag% ddt_tracer_adv(jc,jk,jb,jt)
            END IF
            !
            
          END DO
        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    IF (ltimer) CALL timer_stop(timer_d2p_couple)

    !
    !=====================================================================================

    IF (ltimer)  CALL timer_stop (timer_dyn2phy)

    !=====================================================================================
    !
    ! (3) Prepare boundary conditions for ECHAM physics
    !
    IF (ltimer) CALL timer_start(timer_echam_bcs)

    CALL echam_phy_bcs_global( datetime     ,&! in
      &                        jg           ,&! in
      &                        patch        ,&! in
      &                        dtadv_loc    ,&! in
      &                        ltrig_rad    ,&! out
      &                        datetime_radtran) ! out

    IF (ltimer) CALL timer_stop(timer_echam_bcs)
    !
    !=====================================================================================

    !=====================================================================================
    !
    ! (4) Call echam physics and compute the total physics tendencies.
    !     This includes the atmospheric processes (proper ECHAM) and
    !     the land processes, which are vertically implicitly coupled
    !     to the parameterization of vertical turbulent fluxes.
    !
    IF (ltimer) CALL timer_start(timer_echam_phy)

#ifndef __NO_JSBACH__
    IF (echam_phy_config%ljsbach) THEN
      CALL jsbach_start_timestep(jg)
    END IF
#endif

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
        &                  dtadv_loc    ,&! in
        &                  dtadv_loc    ,&! in
        &                  ltrig_rad     )

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (ltimer) CALL timer_stop(timer_echam_phy)
    !
    !=====================================================================================

#ifndef __NO_JSBACH__
    IF (echam_phy_config%ljsbach) THEN
      CALL jsbach_finish_timestep(jg, dtadv_loc)
    END IF
#endif
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
    ! (6) Convert physics tendencies to dynamics tendencies
    !
    IF (ltimer) CALL timer_start(timer_p2d_prep)
    !
!!$    !     (a) (dT/dt|phy, dqv/dt|phy, dqc/dt|phy, dqi/dt|phy) --> dexner/dt|phy
!!$    !
!!$    ! Loop over cells
!!$!$OMP PARALLEL
!!$!$OMP DO PRIVATE(jb,jk,jc,jcs,jce,z_qsum,z_ddt_qsum) ICON_OMP_DEFAULT_SCHEDULE
!!$    DO jb = i_startblk,i_endblk
!!$      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
!!$      DO jk = 1,nlev
!!$        DO jc = jcs, jce
!!$          z_qsum     = pt_prog_new_rcf%tracer(jc,jk,jb,iqc) + pt_prog_new_rcf%tracer(jc,jk,jb,iqi)
!!$          z_ddt_qsum = prm_tend(jg)% qtrc_phy(jc,jk,jb,iqc) + prm_tend(jg)% qtrc_phy(jc,jk,jb,iqi)
!!$          !
!!$          pt_diag%ddt_exner_phy(jc,jk,jb) =                                               &
!!$            &  rd_o_cpd / pt_prog_new%theta_v(jc,jk,jb)                                   &
!!$            &  * (  prm_tend(jg)%temp_phy(jc,jk,jb)                                       &
!!$            &     * (1._wp + vtmpc1*pt_prog_new_rcf%tracer(jc,jk,jb,iqv) - z_qsum )       &
!!$            &     + pt_diag%temp(jc,jk,jb)                                                &
!!$            &     * (        vtmpc1*prm_tend(jg)% qtrc_phy(jc,jk,jb,iqv) - z_ddt_qsum ) )
!!$          !
!!$          ! Additionally use this loop also to set the dynamical exner increment to zero.
!!$          ! (It is accumulated over one advective time step in solve_nh)
!!$          pt_diag%exner_dyn_incr(jc,jk,jb) = 0._wp
!!$        END DO
!!$      END DO
!!$    END DO !jb
!!$!$OMP END DO
!!$!$OMP END PARALLEL

    !     (b) (du/dt|phy, dv/dt|phy) --> dvn/dt|phy
    !
    ALLOCATE(zdudt(nproma,nlev,patch%nblks_c), &
      &      zdvdt(nproma,nlev,patch%nblks_c), &
      &      stat=return_status)
    IF (return_status > 0) THEN
      CALL finish (module_name//method_name, 'ALLOCATE(zdudt,zdvdt)')
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

    IF (ltimer) CALL timer_stop(timer_p2d_prep)

    ! Now derive the physics-induced normal wind tendency, and add it to the
    ! total tendency.
    IF (ltimer) CALL timer_start(timer_p2d_sync)
    CALL sync_patch_array_mult(SYNC_C, patch, 2, zdudt, zdvdt)
    IF (ltimer) CALL timer_stop(timer_p2d_sync)

    IF (ltimer) CALL timer_start(timer_p2d_prep)

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

          pt_diag%ddt_vn_phy(je,jk,jb) =   pt_int_state%c_lin_e(je,1,jb)*zvn1 &
            &                            + pt_int_state%c_lin_e(je,2,jb)*zvn2

        END DO ! je
      END DO ! jk

    END DO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DEALLOCATE(zdudt, zdvdt)

    IF (ltimer) CALL timer_stop(timer_p2d_prep)
    !
    !=====================================================================================


    !=====================================================================================
    !
    ! (7) Couple dynamics+transport and physics

    IF (ltimer) CALL timer_start(timer_p2d_couple)
    !
    SELECT CASE (echam_phy_config%idcphycpl)

    CASE (1) ! idcphycpl
      ! In this case all ECHAM physics is treated as "fast" physics:
      ! - The provisional "new" state is updated with the total phyiscs
      !   tendencies, providing the final "new" state
      ! - The physics forcing that is passed to the dynamical
      !   core must be set to zero

      ! Loop over edges
      jbs   = patch%edges%start_blk(grf_bdywidth_e+1,1)
      jbe   = patch%nblks_e

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,jes,jee) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,jbe
        CALL get_indices_e(patch, jb,jbs,jbe, jes,jee, grf_bdywidth_e+1)

        DO jk = 1, nlev
          DO je = jes, jee

            ! (1) Velocity
            !
            ! Update with the total phyiscs tendencies
            pt_prog_new%vn    (je,jk,jb) =   pt_prog_new%vn    (je,jk,jb)             &
              &                            + pt_diag%ddt_vn_phy(je,jk,jb) * dtadv_loc
            !
            ! Set physics forcing to zero so that it is not re-applied in the dynamical core
            pt_diag%ddt_vn_phy(je,jk,jb) = 0._wp

          END DO
        END DO

      END DO !jb
!$OMP END DO
!$OMP END PARALLEL

      ! Loop over cells
!$OMP PARALLEL
!$OMP DO PRIVATE(jt,jb,jk,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
      DO jt =1,ntracer    
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
          DO jc = jcs, jce
            prm_field(jg)% mtrcvi    (jc,jb,jt) = 0.0_wp
            prm_tend (jg)% mtrcvi_phy(jc,jb,jt) = 0.0_wp
          END DO
          DO jk = 1,nlev
            DO jc = jcs, jce

              ! (2.1) Tracer mixing ratio with respect to dry air
              !
              ! tracer mass tendency
              IF (echam_phy_config%ldrymoist) THEN
                prm_tend(jg)%   mtrc_phy(jc,jk,jb,jt)  = prm_tend(jg)%  qtrc_phy(jc,jk,jb,jt) &
                  &                                     *prm_field(jg)% mdry    (jc,jk,jb)
              ELSE
                prm_tend(jg)%   mtrc_phy(jc,jk,jb,jt)  = prm_tend(jg)%  qtrc_phy(jc,jk,jb,jt) &
                  &                                     *prm_field(jg)% mair    (jc,jk,jb)
              END IF
              !
              ! tracer path tendency
              prm_tend(jg)% mtrcvi_phy(jc,   jb,jt)  = prm_tend(jg)% mtrcvi_phy(jc,   jb,jt) &
                &                                     +prm_tend(jg)%   mtrc_phy(jc,jk,jb,jt)
              !
              ! new tracer mass
              prm_field(jg)%  mtrc    (jc,jk,jb,jt)  = prm_field(jg)% mtrc    (jc,jk,jb,jt) &
                &                                     +prm_tend(jg)%  mtrc_phy(jc,jk,jb,jt) &
                &                                     *dtadv_loc
              !
              ! new tracer path
              prm_field(jg)%  mtrcvi  (jc,   jb,jt)  = prm_field(jg)% mtrcvi  (jc,   jb,jt) &
                &                                     +prm_field(jg)% mtrc    (jc,jk,jb,jt)
              !
            END DO
          END DO
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL

      ! Loop over cells
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
        DO jc = jcs, jce
          prm_field(jg)% mh2ovi(jc,jb) = 0.0_wp ! initialize air path after physics
          prm_field(jg)% mairvi(jc,jb) = 0.0_wp ! initialize air path after physics
          prm_field(jg)% mdryvi(jc,jb) = 0.0_wp ! initialize air path after physics
        END DO
        DO jk = 1,nlev
          DO jc = jcs, jce

            ! new h2o mass
            prm_field(jg)% mh2o  (jc,jk,jb) = prm_field(jg)%      mtrc (jc,jk,jb,iqv) &
              &                              +prm_field(jg)%      mtrc (jc,jk,jb,iqc) &
              &                              +prm_field(jg)%      mtrc (jc,jk,jb,iqi)
            !
            ! new h2o path
            prm_field(jg)% mh2ovi(jc,   jb) = prm_field(jg)%      mh2ovi(jc,   jb) &
                &                            +prm_field(jg)%      mh2o  (jc,jk,jb)
            !
            IF (echam_phy_config%ldrymoist) THEN
              !
              ! new air mass
              prm_field(jg)% mair  (jc,jk,jb) = prm_field(jg)%      mdry (jc,jk,jb) &
                &                              +prm_field(jg)%      mh2o (jc,jk,jb)
              !
              ! new density
              pt_prog_new %     rho(jc,jk,jb) = prm_field(jg)%      mair  (jc,jk,jb) &
                &                              /prm_field(jg)%      dz    (jc,jk,jb)
              !
            ELSE
              !
              ! new dry air mass
              prm_field(jg)% mdry  (jc,jk,jb) = prm_field(jg)%      mair (jc,jk,jb) &
                &                              -prm_field(jg)%      mh2o (jc,jk,jb)
              !              
            END IF
            !
            ! new air path
            prm_field(jg)% mairvi(jc,   jb) = prm_field(jg)%      mairvi(jc,   jb) &
                &                            +prm_field(jg)%      mair  (jc,jk,jb)
            !
            ! new dry air path
            prm_field(jg)% mdryvi(jc,   jb) = prm_field(jg)%      mdryvi(jc,   jb) &
              &                              +prm_field(jg)%      mdry  (jc,jk,jb)
            !
          END DO
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL

      ! Loop over cells
!$OMP PARALLEL
!$OMP DO PRIVATE(jt,jb,jk,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
      DO jt =1,ntracer    
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
          DO jk = 1,nlev
            DO jc = jcs, jce
              !
              ! new tracer mass fraction with respect to dry air
              IF (echam_phy_config%ldrymoist) THEN
                prm_field(jg)%   qtrc   (jc,jk,jb,jt)  = prm_field(jg)%  mtrc(jc,jk,jb,jt) &
                  &                                     /prm_field(jg)%  mdry(jc,jk,jb)
              ELSE
                prm_field(jg)%   qtrc   (jc,jk,jb,jt)  = prm_field(jg)%  mtrc(jc,jk,jb,jt) &
                  &                                     /prm_field(jg)%  mair(jc,jk,jb)
              END IF
              !
              pt_prog_new_rcf% tracer (jc,jk,jb,jt)  = prm_field(jg)%  mtrc(jc,jk,jb,jt) &
                &                                     /prm_field(jg)%  mair(jc,jk,jb)
              !
            END DO
          END DO
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jcs,jce,z_qsum,z_exner) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

        DO jk = 1,nlev
          DO jc = jcs, jce

            ! (3) Exner function and virtual potential temperature
            !
            ! (a) Update T, then compute Temp_v, Exner and Theta_v
            !
            pt_diag        %temp  (jc,jk,jb    ) =   pt_diag%     temp    (jc,jk,jb)             &
              &                                    + prm_tend(jg)%temp_phy(jc,jk,jb) * dtadv_loc
            !
            z_qsum = pt_prog_new_rcf%tracer(jc,jk,jb,iqc) + pt_prog_new_rcf%tracer(jc,jk,jb,iqi)
            !
            pt_diag%tempv(jc,jk,jb) =   pt_diag%temp(jc,jk,jb)                                             &
              &                       * ( 1._wp +  vtmpc1 * pt_prog_new_rcf%tracer(jc,jk,jb,iqv) - z_qsum)
            !
            ! Save provisional "new" exner from the slow-physics-forced dynamics
            z_exner = pt_prog_new%exner(jc,jk,jb)
            !
            ! Compute final new exner
            pt_prog_new%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd/p0ref * pt_prog_new%rho(jc,jk,jb) * pt_diag%tempv(jc,jk,jb)))
            !
            ! Add exner change from fast phyiscs to exner_old (why?)
            pt_diag%exner_old(jc,jk,jb) = pt_diag%exner_old(jc,jk,jb) + pt_prog_new%exner(jc,jk,jb) - z_exner
            !
!!$            ! (b) Update Exner, then compute Temp_v
!!$            !
!!$            pt_prog_new%exner(jc,jk,jb) = pt_prog_new%exner(jc,jk,jb)                 &
!!$              &                         + pt_diag%ddt_exner_phy(jc,jk,jb) * dtadv_loc
!!$            pt_diag%exner_old(jc,jk,jb) = pt_diag%exner_old(jc,jk,jb)                 &
!!$              &                         + pt_diag%ddt_exner_phy(jc,jk,jb) * dtadv_loc
!!$            !
!!$            pt_diag%tempv(jc,jk,jb) = EXP(LOG(pt_prog_new%exner(jc,jk,jb)/rd_o_cpd)) &
!!$              &                     / (pt_prog_new%rho(jc,jk,jb)*rd/p0ref)
!!$            !
            !
            ! (a) and (b) Compute Theta_v
            !
            pt_prog_new%theta_v(jc,jk,  jb) = pt_diag%tempv(jc,jk,jb) / pt_prog_new%exner(jc,jk,jb)
            !
            ! Set physics forcing to zero so that it is not re-applied in the dynamical core
            pt_diag%ddt_exner_phy(jc,jk,jb) = 0._wp
            !
            ! Additionally use this loop also to set the dynamical exner increment to zero.
            ! (It is accumulated over one advective time step in solve_nh)
            pt_diag%exner_dyn_incr(jc,jk,jb) = 0._wp
            !
        END DO
      END DO

    END DO !jb
!$OMP END DO
!$OMP END PARALLEL

    CASE (2) ! idcphycpl
      ! In this case all ECHAM physics is treated as "slow" physics:
      ! - The full physics forcing has been computed for the final "new" state,
      !   which is the "now" state of the next time step, on which the forcing
      !   shall be applied.
      ! - Hence the full phyiscs forcing is passed on and nothing needs to be
      !   done here.

    END SELECT ! idcphycpl

    IF (ltimer) CALL timer_stop(timer_p2d_couple)
    !
    !=====================================================================================

    !=====================================================================================
    !
    ! Finally do some synchronization for the next dynamics and transport time step(s)
    !
    IF (ltimer) CALL timer_start(timer_p2d_sync)

    CALL sync_patch_array_mult( SYNC_E, patch, 1, pt_prog_new%vn )

    IF      (lhdiff_rcf .AND. diffusion_config(jg)%lhdiff_w) THEN
      CALL sync_patch_array_mult( SYNC_C                       ,&
        &                         patch                        ,&
        &                         ntracer+5                    ,&
        &                         pt_prog_new%w                ,&
        &                         pt_prog_new%rho              ,&
        &                         pt_diag%exner_old            ,&
        &                         pt_prog_new%exner            ,&
        &                         pt_prog_new%theta_v          ,&
        &                         f4din=pt_prog_new_rcf%tracer )
    ELSE IF (lhdiff_rcf) THEN
      CALL sync_patch_array_mult( SYNC_C                       ,&
        &                         patch                        ,&
        &                         ntracer+4                    ,&
        &                         pt_prog_new%rho              ,&
        &                         pt_diag%exner_old            ,&
        &                         pt_prog_new%exner            ,&
        &                         pt_prog_new%theta_v          ,&
        &                         f4din=pt_prog_new_rcf%tracer )
    ELSE
      CALL sync_patch_array_mult( SYNC_C                       ,&
        &                         patch                        ,&
        &                         ntracer+3                    ,&
        &                         pt_prog_new%rho              ,&
        &                         pt_prog_new%exner            ,&
        &                         pt_prog_new%theta_v          ,&
        &                         f4din=pt_prog_new_rcf%tracer )
    ENDIF

    IF (ltimer) CALL timer_stop(timer_p2d_sync)
    !
    !=====================================================================================

    IF (ltimer) CALL timer_stop(timer_phy2dyn)

    !=====================================================================================
    !
    ! Now the final new state (pt_prog_new/pt_prog_new_rcf) and
    ! the slow-physics forcing based on this new state are ready.
    ! The latter is zero if echam_phy_config%idcphycpl=1.
    !
    !=====================================================================================


  END SUBROUTINE interface_iconam_echam
  !----------------------------------------------------------------------------

END MODULE mo_interface_iconam_echam
