!>
!!
!!----------------------------------------------------------------------
!! This module provides the interface between dynamics+transport and
!! echam physics and controls the coupling procedure dependent on the
!! switch mpi_phy_config%idcphycpl.
!!
!! @author Marco Giorgetta (MPI-M)
!!
!!----------------------------------------------------------------------
!!
!! @brief Interface between ICONAM dynamics+transport and ECHAM physics
!!
!! The coupling mechanism is controlled by mpi_phy_config%idcphycpl:
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

  USE mo_coupling_config       ,ONLY: is_coupled_run
  USE mo_parallel_config       ,ONLY: nproma
  USE mo_master_config         ,ONLY: isRestart
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, &
    &                                  calculate_time_interpolation_weights
  USE mo_run_config            ,ONLY: nlev, ntracer, iqv, iqc, iqi, io3
  USE mo_nonhydrostatic_config ,ONLY: lhdiff_rcf
  USE mo_diffusion_config      ,ONLY: diffusion_config
  USE mo_mpi_phy_config        ,ONLY: mpi_phy_config, mpi_phy_tc, dt_zero

  USE mo_model_domain          ,ONLY: t_patch
  USE mo_intp_data_strc        ,ONLY: t_int_state
  USE mo_intp_rbf              ,ONLY: rbf_vec_interpol_cell
  USE mo_loopindices           ,ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants        ,ONLY: min_rlcell_int, grf_bdywidth_c, grf_bdywidth_e
  USE mo_sync                  ,ONLY: sync_c, sync_e, sync_patch_array, sync_patch_array_mult

  USE mo_nonhydro_types        ,ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nh_diagnose_pres_temp ,ONLY: diagnose_pres_temp
  USE mo_physical_constants    ,ONLY: rd, p0ref, rd_o_cpd, vtmpc1, grav, &
    &                                 amd, amo3       
  USE mtime                    ,ONLY: datetime , newDatetime , deallocateDatetime     ,&
    &                                 timedelta, newTimedelta, deallocateTimedelta    ,&
    &                                 max_timedelta_str_len  , getPTStringFromSeconds ,&
    &                                 OPERATOR(+), OPERATOR(>)

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

  USE mo_lcariolle_types         ,ONLY: l_cariolle_initialized_o3, t_avi, t_time_interpolation
  USE mo_linked_list,             ONLY: t_var_list
  USE mo_ext_data_state,          ONLY: ext_data
  USE mo_art_reaction_interface,  ONLY: art_reaction_interface
  USE mo_run_config,              ONLY: lart

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
  SUBROUTINE interface_iconam_echam( dt_loc          ,& !in
    &                                datetime_new    ,& !in
    &                                patch           ,& !in
    &                                pt_int_state    ,& !in
    &                                p_metrics       ,& !in
    &                                pt_prog_old     ,& !in
    &                                pt_prog_old_rcf ,& !in
    &                                pt_prog_new     ,& !inout
    &                                pt_prog_new_rcf ,& !inout
    &                                pt_diag         ,& !inout
    &                                p_prog_list)

    !
    !> Arguments:
    !
    REAL(wp)              , INTENT(in)            :: dt_loc          !< advective time step
    TYPE(datetime)        , POINTER               :: datetime_new    !< date and time at the end of this time step

    TYPE(t_patch)         , INTENT(in)   , TARGET :: patch           !< grid/patch info
    TYPE(t_int_state)     , INTENT(in)   , TARGET :: pt_int_state    !< interpolation state
    TYPE(t_nh_metrics)    , INTENT(in)            :: p_metrics

    TYPE(t_nh_diag)       , INTENT(inout), TARGET :: pt_diag         !< diagnostic variables
    TYPE(t_nh_prog)       , INTENT(inout), TARGET :: pt_prog_old     !< progn. vars before dynamics  for wind, temp. rho, ...
    TYPE(t_nh_prog)       , INTENT(inout), TARGET :: pt_prog_old_rcf !< progn. vars before advection for tracers
    TYPE(t_nh_prog)       , INTENT(inout), TARGET :: pt_prog_new     !< progn. vars after dynamics  for wind, temp. rho, ...
    TYPE(t_nh_prog)       , INTENT(inout), TARGET :: pt_prog_new_rcf !< progn. vars after advection for tracers

!ICON_ART
    TYPE(t_var_list), OPTIONAL, INTENT(in)        :: p_prog_list     !< current prognostic state list

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

    CHARACTER(len=max_timedelta_str_len) :: neg_dt_loc_string !< negative time delta as string
    TYPE(timedelta), POINTER             :: neg_dt_loc_mtime  !< negative time delta as mtime variable
    TYPE(datetime) , POINTER             :: datetime_old      !< date and time at the beginning of this time step

    REAL(wp) :: z_exner              !< to save provisional new exner
    REAL(wp) :: z_qsum               !< summand of virtual increment
!!$    REAL(wp) :: z_ddt_qsum           !< summand of virtual increment

    REAL(wp) :: zvn1, zvn2
    REAL(wp), POINTER :: zdudt(:,:,:), zdvdt(:,:,:)

    TYPE(t_avi) :: avi

    INTEGER  :: return_status

    ! Local parameters

    CHARACTER(*), PARAMETER :: method_name = "interface_iconam_echam"

    ! Temporary variables for Cariolle scheme (ozone)
    REAL(wp)    :: vmr_o3(nproma,nlev)
    TYPE(t_time_interpolation) :: time_interpolation
    EXTERNAL       lcariolle_lat_intp_li, lcariolle_pres_intp_li
    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights
    !-------------------------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_dyn2phy)

    IF (ltimer) CALL timer_start(timer_d2p_prep)

    ! Inquire current grid level and the total number of grid cells
    i_nchdom = MAX(1,patch%n_childdom)
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    jg    = patch%id

    ! The date and time needed for the radiation computation in the phyiscs is
    ! the date and time of the initial data for this step.
    ! As 'datetime_new' contains already the date and time of the end of this
    ! time step, we compute here the old datetime 'datetime_old':
    !
    CALL getPTStringFromSeconds(-dt_loc, neg_dt_loc_string)
    neg_dt_loc_mtime  => newTimedelta(neg_dt_loc_string)
    datetime_old      => newDatetime(datetime_new)
    datetime_old      =  datetime_new + neg_dt_loc_mtime
    CALL deallocateTimedelta(neg_dt_loc_mtime)

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
    SELECT CASE (mpi_phy_config(jg)%idcphycpl)
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
                &                                   + prm_tend(jg)% qtrc_phy(jc,jk,jb,jt) * dt_loc
              !
            END DO
          END DO
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
      !
    END SELECT ! idcphycpl
    !
    ! Now the new prognostic state variables (pt_prog_new and pt_prog_new_rcf) of
    ! the dynamics and transport are complete and ready to be used in the phyiscs.
    !
    ! idcphycpl = 1: The "new" state is provisional, updated only by dynamics,
    !                diffusion and tracer transport.
    !                In the following the phyiscs forcing is computed for this new
    !                provisional state and the provisional new state is updated
    !                to obtain the final "new" state X(t+dt).
    !
    ! idcphycpl = 2: The "new" state is the final state X(t+dt).
    !                In the following the phyiscs forcing is computed for this new
    !                final state so that it is available in the next time step.
    !
    !=====================================================================================


    !=====================================================================================
    !
    ! (2) Diagnostics
    !
    ! - pt_diag%tempv
    ! - pt_diag%temp
    ! - pt_diag%pres_sfc   surface pressure filtered to remove sound waves, see diagnose_pres_temp
    ! - pt_diag%pres_ifc   hydrostatic pressure at layer interface
    ! - pt_diag%pres       hydrostatic pressure at layer midpoint = SQRT(upper pres_ifc * lower pres_ifc)
    ! - pt_diag%dpres_mc   pressure thickness of layer
    !
    ! For the old state:
    !
    CALL diagnose_pres_temp( p_metrics                ,&
      &                      pt_prog_old              ,&
      &                      pt_prog_old_rcf          ,&
      &                      pt_diag                  ,&
      &                      patch                    ,&
      &                      opt_calc_temp=.TRUE.     ,&
      &                      opt_calc_pres=.TRUE.     ,&
      &                      opt_rlend=min_rlcell_int )

    IF (ltimer) CALL timer_stop(timer_d2p_prep)

    ! - pt_diag%u
    ! - pt_diag%v
    IF (ltimer) CALL timer_start(timer_d2p_sync)
    CALL sync_patch_array( SYNC_E, patch, pt_prog_old%vn )
    IF (ltimer) CALL timer_stop(timer_d2p_sync)

    IF (ltimer) CALL timer_start(timer_d2p_prep)

    CALL rbf_vec_interpol_cell( pt_prog_old%vn       ,&! in
      &                         patch                ,&! in
      &                         pt_int_state         ,&! in
      &                         pt_diag%u            ,&! out
      &                         pt_diag%v            ,&! out
      &                         opt_rlstart=rl_start ,&! in
      &                         opt_rlend  =rl_end   ) ! in

    !
    ! Store old diagnosed fields provisionally for diagnosing later the trends:
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      DO jk = 1,nlev
        DO jc = jcs, jce
          prm_tend(jg)% ua(jc,jk,jb) = pt_diag% u   (jc,jk,jb)
          prm_tend(jg)% va(jc,jk,jb) = pt_diag% v   (jc,jk,jb)
          prm_tend(jg)% ta(jc,jk,jb) = pt_diag% temp(jc,jk,jb)
        END DO
      END DO
    END DO ! jb
!$OMP END DO
!$OMP END PARALLEL
    !
    !
    ! For the new state:
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
    ! Now the old and new prognostic and diagnostic state variables of the dynamical core
    ! are ready to be used in the phyiscs.
    !
    !=====================================================================================


    !=====================================================================================
    !
    ! (3) Copy the new prognostic state and the related diagnostics from the
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

          ! Fill the time dependent physics state variables, which are used by echam:
          !
          prm_field(jg)%        ua(jc,jk,jb)     = pt_diag%    u(jc,jk,jb)
          prm_field(jg)%        va(jc,jk,jb)     = pt_diag%    v(jc,jk,jb)
          prm_field(jg)%       vor(jc,jk,jb)     = pt_diag%  vor(jc,jk,jb)
          !
          prm_field(jg)%        ta(jc,jk,jb)     = pt_diag% temp(jc,jk,jb)
          prm_field(jg)%        tv(jc,jk,jb)     = pt_diag% tempv(jc,jk,jb)
          !
          prm_field(jg)% presm_old(jc,jk,jb)     = pt_diag%  pres(jc,jk,jb)
          prm_field(jg)% presm_new(jc,jk,jb)     = pt_diag%  pres(jc,jk,jb)
          !
          ! density
          prm_field(jg)%       rho(jc,jk,jb)     = pt_prog_new %         rho(jc,jk,jb)
          !
          ! air mass
          prm_field(jg)%      mair(jc,jk,jb)     = pt_prog_new %         rho(jc,jk,jb) &
            &                                     *prm_field(jg)%         dz(jc,jk,jb)
          !
          ! H2O mass (vap+liq+ice)
          prm_field(jg)%      mh2o(jc,jk,jb)     = ( pt_prog_new_rcf% tracer(jc,jk,jb,iqv)  &
            &                                       +pt_prog_new_rcf% tracer(jc,jk,jb,iqc)  &
            &                                       +pt_prog_new_rcf% tracer(jc,jk,jb,iqi)) &
            &                                      *prm_field(jg)%      mair(jc,jk,jb)
          !
          ! dry air mass
          prm_field(jg)%      mdry(jc,jk,jb)     = prm_field(jg)%       mair(jc,jk,jb) &
            &                                     -prm_field(jg)%       mh2o(jc,jk,jb)
          !
          ! cloud water+ice
          IF (mpi_phy_config(jg)%ldrymoist) THEN
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
          ! Diagnose tendencies from the new and old diagnostic states and the local time step:
          !
          SELECT CASE (mpi_phy_config(jg)%idcphycpl)
            !
          CASE (1) ! idcphycpl
             ! In this case the new state is provisional and updated only by dynamics:
             ! The dynamical tendency is diagnosed as the difference of the new and old(="now") state.
             ! The old state for ua, va, and ta is provisionally stored in the tendency variables:
             prm_tend(jg)% ua_dyn(jc,jk,jb)     = (pt_diag%u   (jc,jk,jb)-prm_tend(jg)%ua(jc,jk,jb))/dt_loc
             prm_tend(jg)% va_dyn(jc,jk,jb)     = (pt_diag%v   (jc,jk,jb)-prm_tend(jg)%va(jc,jk,jb))/dt_loc
             prm_tend(jg)% ta_dyn(jc,jk,jb)     = (pt_diag%temp(jc,jk,jb)-prm_tend(jg)%ta(jc,jk,jb))/dt_loc
             !
             ! Initialize the total tendencies, to be computed later:
             prm_tend(jg)% ua    (jc,jk,jb)     = 0.0_wp
             prm_tend(jg)% va    (jc,jk,jb)     = 0.0_wp
             prm_tend(jg)% ta    (jc,jk,jb)     = 0.0_wp
             !
             ! Now reset the physics tendencies before entering the physics:
             prm_tend(jg)% ua_phy(jc,jk,jb)     = 0.0_wp
             prm_tend(jg)% va_phy(jc,jk,jb)     = 0.0_wp
             prm_tend(jg)% ta_phy(jc,jk,jb)     = 0.0_wp
             !
          CASE(2) ! idcphycpl
             ! In this case the new state is final:
             ! The total tendency is diagnosed as the difference of the new and old(="now") state.
             ! The old state for ua, va, and ta is provisionally stored in the tendency variables:
             prm_tend(jg)% ua    (jc,jk,jb)     = (pt_diag%u   (jc,jk,jb)-prm_tend(jg)%ua(jc,jk,jb))/dt_loc
             prm_tend(jg)% va    (jc,jk,jb)     = (pt_diag%v   (jc,jk,jb)-prm_tend(jg)%va(jc,jk,jb))/dt_loc
             prm_tend(jg)% ta    (jc,jk,jb)     = (pt_diag%temp(jc,jk,jb)-prm_tend(jg)%ta(jc,jk,jb))/dt_loc
             !
             ! The dynamic tendency can be diagnosed from the total and the old physics tendencies:
             prm_tend(jg)% ua_dyn(jc,jk,jb)     = prm_tend(jg)% ua(jc,jk,jb)-prm_tend(jg)% ua_phy(jc,jk,jb)
             prm_tend(jg)% va_dyn(jc,jk,jb)     = prm_tend(jg)% va(jc,jk,jb)-prm_tend(jg)% va_phy(jc,jk,jb)
             prm_tend(jg)% ta_dyn(jc,jk,jb)     = prm_tend(jg)% ta(jc,jk,jb)-prm_tend(jg)% ta_phy(jc,jk,jb)
             !
             ! Now reset the physics tendencies before entering the physics:
             prm_tend(jg)% ua_phy(jc,jk,jb)     = 0.0_wp
             prm_tend(jg)% va_phy(jc,jk,jb)     = 0.0_wp
             prm_tend(jg)% ta_phy(jc,jk,jb)     = 0.0_wp
             !
          END SELECT

        END DO
      END DO

      DO jk = 1,nlev+1
        DO jc = jcs, jce

          prm_field(jg)%     presi_old(jc,jk,jb) = pt_diag% pres_ifc(jc,jk,jb)
          prm_field(jg)%     presi_new(jc,jk,jb) = pt_diag% pres_ifc(jc,jk,jb)

        END DO
      END DO

    END DO ! jb
!$OMP END DO
!$OMP END PARALLEL

!   Initialize ozone mass mixing ratios for Cariolle scheme here. 
!   An approximative initialization 
!   that considers the atmosphere as being dry is enough.
    IF (mpi_phy_tc(jg)%dt_car > dt_zero) THEN
      IF (.NOT.isRestart().AND. .NOT. l_cariolle_initialized_o3) THEN
        ALLOCATE(avi%cell_center_lat(nproma))
        avi%ldown=.TRUE.
        current_time_interpolation_weights = calculate_time_interpolation_weights(datetime_old)
        time_interpolation%imonth1=current_time_interpolation_weights%month1_index
        time_interpolation%imonth2=current_time_interpolation_weights%month2_index
        time_interpolation%weight1=current_time_interpolation_weights%weight1
        time_interpolation%weight2=current_time_interpolation_weights%weight2
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
          avi%pres                     => prm_field(jg)%presm_old(:,:,jb)
          avi%cell_center_lat(jcs:jce) =  prm_field(jg)%clat(jcs:jce,jb)
          CALL lcariolle_init_o3(                                              &
           & jcs,                   jce,                nproma,                &
           & nlev,                  time_interpolation, lcariolle_lat_intp_li, &
           & lcariolle_pres_intp_li,avi,                vmr_o3                 )
          pt_prog_new_rcf% tracer(jcs:jce,:,jb,io3)=vmr_o3(jcs:jce,:)*amo3/amd
        END DO
        l_cariolle_initialized_o3 = .TRUE.
        DEALLOCATE(avi%cell_center_lat)
      END IF
    END IF

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
            IF (mpi_phy_config(jg)%ldrymoist) THEN
              prm_field(jg)%      qtrc(jc,jk,jb,jt)  = pt_prog_new_rcf% tracer(jc,jk,jb,jt) &
                &                                     *prm_field(jg)%     mair(jc,jk,jb)    &
                &                                     /prm_field(jg)%     mdry(jc,jk,jb)
            ELSE
              prm_field(jg)%      qtrc(jc,jk,jb,jt)  = pt_prog_new_rcf% tracer(jc,jk,jb,jt)
            END IF

            SELECT CASE (mpi_phy_config(jg)%idcphycpl)
               !
               ! Diagnose tendencies from the new and old states and the local time step:
               !
            CASE (1) ! idcphycpl
               ! In this case the new state is provisional and updated only by dynamics:
               ! The dynamical tendency is diagnosed as the difference of the new and old(="now") state.
               IF (mpi_phy_config(jg)%ldrymoist) THEN
                 prm_tend(jg)%   qtrc_dyn(jc,jk,jb,jt)  = ( pt_prog_new_rcf% tracer(jc,jk,jb,jt)          &
                      &                                    -pt_prog_old_rcf% tracer(jc,jk,jb,jt) )/dt_loc &
                      &                                   *prm_field(jg)%      mair(jc,jk,jb)             &
                      &                                   /prm_field(jg)%      mdry(jc,jk,jb)
                 !
                 !-FOR TESTING-------------------------------------------------------------------!
                 prm_tend(jg)%   qtrc_dyn(jc,jk,jb,jt)  = pt_diag% ddt_tracer_adv(jc,jk,jb,jt) & !
                   &                                     *prm_field(jg)%     mair(jc,jk,jb)    & !
                   &                                     /prm_field(jg)%     mdry(jc,jk,jb)      !
                 !-------------------------------------------------------------------------------!
                 !
              ELSE
                 prm_tend(jg)%   qtrc_dyn(jc,jk,jb,jt)  = ( pt_prog_new_rcf% tracer(jc,jk,jb,jt)          &
                      &                                    -pt_prog_old_rcf% tracer(jc,jk,jb,jt) )/dt_loc
                 !
                 !-FOR TESTING-------------------------------------------------------------------!
                 prm_tend(jg)%   qtrc_dyn(jc,jk,jb,jt)  = pt_diag% ddt_tracer_adv(jc,jk,jb,jt)   !
                 !-------------------------------------------------------------------------------!
                 !
               END IF
               !
               ! Initialize the total tendencies, to be computed later:
               prm_tend(jg)%     qtrc    (jc,jk,jb,jt)  = 0.0_wp
               !
               ! Now reset the physics tendencies before entering the physics
               prm_tend(jg)%     qtrc_phy(jc,jk,jb,jt)  = 0.0_wp
               !
            CASE (2) ! idcphycpl
               ! In this case the new state is final:
               ! The total tendency is diagnosed as the difference of the new and old(="now") state.
               IF (mpi_phy_config(jg)%ldrymoist) THEN
                  prm_tend(jg)%   qtrc    (jc,jk,jb,jt)  = ( pt_prog_new_rcf% tracer(jc,jk,jb,jt)          &
                       &                                    -pt_prog_old_rcf% tracer(jc,jk,jb,jt) )/dt_loc &
                       &                                   *prm_field(jg)%     mair(jc,jk,jb)              &
                       &                                   /prm_field(jg)%     mdry(jc,jk,jb)
               ELSE
                  prm_tend(jg)%   qtrc    (jc,jk,jb,jt)  = ( pt_prog_new_rcf% tracer(jc,jk,jb,jt)          &
                       &                                    -pt_prog_new_rcf% tracer(jc,jk,jb,jt) )/dt_loc
               END IF
               !
               ! And the dynamic tendency can be diagnosed from the total and the old physics tendencies:
               prm_tend(jg)%      qtrc_dyn(jc,jk,jb,jt)  =  prm_tend(jg)% qtrc    (jc,jk,jb,jt)  &
                    &                                      -prm_tend(jg)% qtrc_phy(jc,jk,jb,jt)
               !
               ! Now reset the physics tendencies before entering the physics
               prm_tend(jg)%      qtrc_phy(jc,jk,jb,jt)  = 0.0_wp
               !
            END SELECT
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

    CALL echam_phy_bcs_global( datetime_old ,&! in
      &                        jg           ,&! in
      &                        patch        ,&! in
      &                        dt_loc       ) ! out

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
#ifndef __NO_JSBACH__
    IF (mpi_phy_config(jg)%ljsb) THEN
      CALL jsbach_start_timestep(jg)
    END IF
#endif

    IF (ltimer) CALL timer_start(timer_echam_phy)

    ! Like in ECHAM, the subroutine *echam_phy_main* has direct access to the memory
    ! buffers prm_field and prm_tend. 

    CALL echam_phy_main( patch,           &! in
      &                  rl_start, rl_end,&! in  
      &                  datetime_old    ,&! in
      &                  dt_loc          ) ! in

    IF (ltimer) CALL timer_stop(timer_echam_phy)

    CALL deallocateDatetime(datetime_old)
    !
    !=====================================================================================

#ifndef __NO_JSBACH__
    IF (mpi_phy_config(jg)%ljsb) THEN
      CALL jsbach_finish_timestep(jg, dt_loc)
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
      zdudt(jcs:jce,:,jb) = prm_tend(jg)% ua_phy(jcs:jce,:,jb)
      zdvdt(jcs:jce,:,jb) = prm_tend(jg)% va_phy(jcs:jce,:,jb)
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (ltimer) CALL timer_stop(timer_p2d_prep)

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
    SELECT CASE (mpi_phy_config(jg)%idcphycpl)

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
              &                            + pt_diag%ddt_vn_phy(je,jk,jb) * dt_loc
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

              ! Diagnose the total tendencies
              prm_tend(jg)%qtrc(jc,jk,jb,jt) =   prm_tend(jg)%qtrc_dyn(jc,jk,jb,jt)  &
                &                              + prm_tend(jg)%qtrc_phy(jc,jk,jb,jt)

              ! (2.1) Tracer mixing ratio with respect to dry air
              !
              ! tracer mass tendency
              IF (mpi_phy_config(jg)%ldrymoist) THEN
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
                &                                     *dt_loc
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
          ! initialize vertical integrals
          prm_field(jg)% mh2ovi(jc,jb) = 0.0_wp
          prm_field(jg)% mairvi(jc,jb) = 0.0_wp
          prm_field(jg)% mdryvi(jc,jb) = 0.0_wp
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
            IF (mpi_phy_config(jg)%ldrymoist) THEN
              !
              ! new air mass
              prm_field(jg)% mair  (jc,jk,jb) = prm_field(jg)%      mdry (jc,jk,jb) &
                &                              +prm_field(jg)%      mh2o (jc,jk,jb)
              !
              ! new density
              prm_field(jg)%    rho(jc,jk,jb) = prm_field(jg)%      mair  (jc,jk,jb) &
                &                              /prm_field(jg)%      dz    (jc,jk,jb)
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
              IF (mpi_phy_config(jg)%ldrymoist) THEN
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

            ! Diagnose the total tendencies
            prm_tend(jg)%ua(jc,jk,jb) = prm_tend(jg)%ua_dyn(jc,jk,jb) + prm_tend(jg)%ua_phy(jc,jk,jb)
            prm_tend(jg)%va(jc,jk,jb) = prm_tend(jg)%va_dyn(jc,jk,jb) + prm_tend(jg)%va_phy(jc,jk,jb)
            prm_tend(jg)%ta(jc,jk,jb) = prm_tend(jg)%ta_dyn(jc,jk,jb) + prm_tend(jg)%ta_phy(jc,jk,jb)
            !
            ! (3) Exner function and virtual potential temperature
            !
            ! (a) Update T, then compute Temp_v, Exner and Theta_v
            !
            pt_diag        %temp  (jc,jk,jb    ) =   pt_diag%     temp    (jc,jk,jb)             &
              &                                    + prm_tend(jg)%ta_phy  (jc,jk,jb) * dt_loc
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
            ! Add exner change from fast phyiscs to exner_pr in order to avoid unphysical sound wave generation
            pt_diag%exner_pr(jc,jk,jb) = pt_diag%exner_pr(jc,jk,jb) + pt_prog_new%exner(jc,jk,jb) - z_exner
            !
!!$            ! (b) Update Exner, then compute Temp_v
!!$            !
!!$            pt_prog_new%exner(jc,jk,jb) = pt_prog_new%exner(jc,jk,jb)                 &
!!$              &                         + pt_diag%ddt_exner_phy(jc,jk,jb) * dt_loc
!!$            pt_diag%exner_old(jc,jk,jb) = pt_diag%exner_old(jc,jk,jb)                 &
!!$              &                         + pt_diag%ddt_exner_phy(jc,jk,jb) * dt_loc
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
        &                         pt_diag%exner_pr             ,&
        &                         pt_prog_new%exner            ,&
        &                         pt_prog_new%theta_v          ,&
        &                         f4din=pt_prog_new_rcf%tracer )
    ELSE IF (lhdiff_rcf) THEN
      CALL sync_patch_array_mult( SYNC_C                       ,&
        &                         patch                        ,&
        &                         ntracer+4                    ,&
        &                         pt_prog_new%rho              ,&
        &                         pt_diag%exner_pr             ,&
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
    ! The latter is zero if mpi_phy_config%idcphycpl=1.
    !
    !=====================================================================================

    IF (mpi_phy_tc(jg)%dt_art > dt_zero) THEN
       CALL art_reaction_interface(ext_data(jg),           & !> in
            &                      patch,                  & !> in
            &                      datetime_new,           & !> in
            &                      dt_loc,                 & !> in
            &                      p_prog_list,            & !> in
            &                      pt_prog_new,            &
            &                      p_metrics,              & !> in
            &                      pt_diag,                & !> inout
            &                      pt_prog_new_rcf%tracer)
    ENDIF

  END SUBROUTINE interface_iconam_echam
  !----------------------------------------------------------------------------

END MODULE mo_interface_iconam_echam
