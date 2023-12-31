!>
!!
!!------------------------------------------------------------------------------
!! This module provides the interface between dynamics+transport and physics.
!!
!! @author Marco Giorgetta (MPI-M)
!!
!!------------------------------------------------------------------------------
!!
!! @brief Interface between ICONAM dynamics+transport and ECHAM physics
!!
!! The whole physics is treated as "fast" phyiscs.
!! The physics tendencies are computed from the provisional state reached
!! after dynamics+transport and the full physics tendencies are then used
!! to update and reach the final new state
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

  USE mo_kind                  ,ONLY: wp, vp
  USE mo_exception             ,ONLY: warning, finish, print_value

#ifdef YAC_coupling
  USE mo_coupling_config       ,ONLY: is_coupled_run
#endif
  USE mo_parallel_config       ,ONLY: nproma
  USE mo_advection_config      ,ONLY: advection_config
  USE mo_run_config            ,ONLY: nlev, ntracer, iqv, iqc, iqi
  USE mo_nonhydrostatic_config ,ONLY: lhdiff_rcf
  USE mo_diffusion_config      ,ONLY: diffusion_config
  USE mo_echam_phy_config      ,ONLY: echam_phy_config

  USE mo_model_domain          ,ONLY: t_patch
  USE mo_intp_data_strc        ,ONLY: t_int_state
  USE mo_intp_rbf              ,ONLY: rbf_vec_interpol_cell
  USE mo_loopindices           ,ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants        ,ONLY: min_rlcell_int, grf_bdywidth_c, &
    &                                 min_rledge_int, grf_bdywidth_e
  USE mo_sync                  ,ONLY: sync_c, sync_e, sync_patch_array, sync_patch_array_mult

  USE mo_nonhydro_types        ,ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nh_diagnose_pres_temp ,ONLY: diagnose_pres_temp
  USE mo_math_constants        ,ONLY: rad2deg
  USE mo_physical_constants    ,ONLY: rd, p0ref, rd_o_cpd, vtmpc1, grav
  USE mtime                    ,ONLY: datetime , newDatetime , deallocateDatetime     ,&
    &                                 timedelta, newTimedelta, deallocateTimedelta    ,&
    &                                 max_timedelta_str_len  , getPTStringFromSeconds ,&
    &                                 OPERATOR(+), OPERATOR(>)

  USE mo_echam_phy_memory      ,ONLY: t_echam_phy_field, prm_field, &
    &                                 t_echam_phy_tend , prm_tend
  USE mo_echam_phy_bcs         ,ONLY: echam_phy_bcs
  USE mo_echam_phy_main        ,ONLY: echam_phy_main
#if defined(YAC_coupling) && !defined(__NO_ECHAM__)
  USE mo_echam_coupling        ,ONLY: interface_echam_ocean
#endif
  
#ifndef __NO_JSBACH__
  USE mo_jsb_interface         ,ONLY: jsbach_start_timestep, jsbach_finish_timestep
#endif
  
  USE mo_timer                 ,ONLY: ltimer, timer_start, timer_stop,                                 &
    &                                 timer_dyn2phy, timer_d2p_prep, timer_d2p_sync, timer_d2p_couple, &
    &                                 timer_echam_bcs, timer_echam_phy, timer_coupling,                &
    &                                 timer_phy2dyn, timer_p2d_prep, timer_p2d_sync, timer_p2d_couple
  USE mo_run_config,            ONLY: lart
#if defined( _OPENACC )
  USE mo_var_list_gpu          ,ONLY: gpu_update_var_list
#endif

  USE mo_upatmo_config         ,ONLY: upatmo_config
  USE mo_upatmo_impl_const,     ONLY: idamtr

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
  SUBROUTINE interface_iconam_echam( dt_loc          & !in
    &                               ,datetime_new    & !in
    &                               ,patch           & !in
    &                               ,pt_int_state    & !in
    &                               ,p_metrics       & !in
    &                               ,pt_prog_new     & !inout
    &                               ,pt_prog_new_rcf & !inout
    &                               ,pt_diag         ) !inout

    !
    !> Arguments:
    !
    REAL(wp)              , INTENT(in)            :: dt_loc          !< advective time step
    TYPE(datetime)        , POINTER               :: datetime_new    !< date and time at the end of this time step

    TYPE(t_patch)         , INTENT(inout), TARGET :: patch           !< grid/patch info
    TYPE(t_int_state)     , INTENT(in)   , TARGET :: pt_int_state    !< interpolation state
    TYPE(t_nh_metrics)    , INTENT(in)            :: p_metrics

    TYPE(t_nh_diag)       , INTENT(inout), TARGET :: pt_diag         !< diagnostic variables
    TYPE(t_nh_prog)       , INTENT(inout), TARGET :: pt_prog_new     !< progn. vars after dynamics  for wind, temp. rho, ...
    TYPE(t_nh_prog)       , INTENT(inout), TARGET :: pt_prog_new_rcf !< progn. vars after advection for tracers

    ! Local array bounds

    INTEGER  :: ncd                  !< number of child patches

    INTEGER  :: rls_c, rle_c
    INTEGER  :: jbs_c, jbe_c         !< start and end indices for rows of cells
    INTEGER  :: jcs,jce              !< cell start and end indices

    INTEGER  :: rls_e, rle_e
    INTEGER  :: jbs_e, jbe_e         !< start and end indices for rows of edges
    INTEGER  :: jes,jee              !< edge start and end indices

    INTEGER  :: jcn,jbn              !< jc and jb of neighbor cells sharing an edge je

    INTEGER  :: jg                   !< grid   index
    INTEGER  :: jb                   !< block  index
    INTEGER  :: jc                   !< cell   index
    INTEGER  :: je                   !< edge   index
    INTEGER  :: jk                   !< level  index
    INTEGER  :: jt                   !< tracer index

    ! Local variables

    CHARACTER(len=max_timedelta_str_len) :: neg_dt_loc_string !< negative time delta as string
    TYPE(timedelta)         , POINTER    :: neg_dt_loc_mtime  !< negative time delta as mtime variable
    TYPE(datetime)          , POINTER    :: datetime_old      !< date and time at the beginning of this time step

    TYPE(t_echam_phy_field) , POINTER    :: field
    TYPE(t_echam_phy_tend)  , POINTER    :: tend
 
    REAL(wp) :: z_exner              !< to save provisional new exner
    REAL(wp) :: z_qsum               !< summand of virtual increment
!!$    REAL(wp) :: z_ddt_qsum           !< summand of virtual increment

    REAL(wp) :: zvn1, zvn2
    REAL(wp), POINTER :: zdudt(:,:,:), zdvdt(:,:,:)

    INTEGER  :: return_status

    ! (For deep-atmosphere modification)
    REAL(wp) :: deepatmo_vol(patch%nlev)

    ! Local parameters

    CHARACTER(*), PARAMETER :: method_name = "interface_iconam_echam"

    INTEGER :: jt_end

    !-------------------------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_dyn2phy)

    IF (ltimer) CALL timer_start(timer_d2p_prep)

    ! Inquire current grid level and the total number of grid cells
    ncd = MAX(1,patch%n_childdom)

    ! cells
    rls_c = grf_bdywidth_c+1
    rle_c = min_rlcell_int
    jbs_c = patch%cells%start_blk(rls_c,  1)
    jbe_c = patch%cells%  end_blk(rle_c,ncd)

    ! edges
    rls_e = grf_bdywidth_e+1
    rle_e = min_rledge_int
    jbs_e = patch%edges%start_blk(rls_e,  1)
    jbe_e = patch%edges%  end_blk(rle_e,ncd)

    jg    = patch%id

    ! associate pointers
    field => prm_field(jg)
    tend  => prm_tend (jg)

    ! The date and time needed for the radiation computation in the physics is
    ! the date and time of the initial data for this step.
    ! As 'datetime_new' contains already the date and time of the end of this
    ! time step, we compute here the old datetime 'datetime_old':
    !
    CALL getPTStringFromSeconds(-dt_loc, neg_dt_loc_string)
    neg_dt_loc_mtime  => newTimedelta(neg_dt_loc_string)
    datetime_old      => newDatetime(datetime_new)
    datetime_old      =  datetime_new + neg_dt_loc_mtime
    CALL deallocateTimedelta(neg_dt_loc_mtime)

    !$ACC DATA PRESENT( pt_prog_new%vn, pt_prog_new%w, pt_prog_new%rho,                         &
    !$ACC               pt_prog_new%exner, pt_prog_new%theta_v,                                 &
    !$ACC               pt_prog_new_rcf%tracer,                                                 &
    !$ACC               pt_diag%u, pt_diag%v, pt_diag%temp, pt_diag%tempv,                      &
    !$ACC               pt_diag%ddt_tracer_adv,                                                 &
    !$ACC               pt_diag%ddt_vn_phy, pt_diag%exner_pr, pt_diag%ddt_exner_phy,            &
    !$ACC               pt_diag%exner_dyn_incr,                                                 &
    !$ACC               pt_int_state%c_lin_e,  advection_config(jg)%trHydroMass%list,           &
    !$ACC               patch%edges%cell_idx, patch%edges%primal_normal_cell,                   &
    !$ACC               field%pfull,                                                            &
    !$ACC               field%rho, field%mair, field%dz, field%mh2o,                            &
    !$ACC               field%mdry, field%mref, field%xref, field%wa, field%omega,              &
    !$ACC               field%clon, field%clat, field%mtrc, field%qtrc,                         &
    !$ACC               field%mtrcvi, field%mh2ovi, field%mairvi, field%mdryvi, field%mrefvi,   &
    !$ACC               tend%ua_phy, tend%va_phy, tend%ta_phy, tend%qtrc, tend%qtrc_dyn,        &
    !$ACC               tend%qtrc_phy, tend%mtrc_phy, tend%mtrcvi_phy, p_metrics%deepatmo_t1mc )&
    !$ACC       CREATE( deepatmo_vol )                                                          &
    !$ACC       COPYIN( echam_phy_config(jg:jg) )

    jt_end = ntracer

    ! Preparation for deep-atmosphere modifications:
    ! - We could do without deepatmo_vol and access p_metrics%deepatmo_t1mc directly. 
    !   However, we introduced deepatmo_vol as a further safety barrier.
    ! - Independent of l_shallowatmo the modifications are computationally expensive!
    !   The computation of field%mair and field/pt_prog_new%rho include 
    !   an additional multiplication or division, respectively.
    ! - The computational overhead could be avoided by the implementation 
    !   of an additional metric 3d-array in field, which contains the values 
    !   of the product field%dz(jc,jk,jb) * p_metrics%deepatmo_t1mc(jk,idamtr%t1mc%vol). 
    !   However, at the cost of the considerable memory consumption of an additional 3d-array.
    IF (upatmo_config(jg)%echam_phy%l_shallowatmo) THEN
      ! no cell volume modification
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR
      DO jk = 1, patch%nlev
        deepatmo_vol(jk) = 1._wp
      END DO
      !$ACC END PARALLEL
    ELSE
      ! cell volume modification factors from 'p_metrics'
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR
      DO jk = 1, patch%nlev
        deepatmo_vol(jk) = p_metrics%deepatmo_t1mc(jk,idamtr%t1mc%vol)
      END DO
      !$ACC END PARALLEL
    END IF

    !=====================================================================================
    !
    ! The "new" state is provisional, updated only by dynamics,
    ! diffusion and tracer transport.
    ! In the following the phyiscs forcing is computed for this new
    ! provisional state and the provisional new state is updated
    ! to obtain the final "new" state X(t+dt).
    !
    !=====================================================================================

    !=====================================================================================
    !
    ! (1) Handling of negative tracer mass fractions resulting from dynamics
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jt,jb,jk,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs_c,jbe_c
      !
      CALL get_indices_c(patch, jb,jbs_c,jbe_c, jcs,jce, rls_c,rle_c)
      IF (jcs>jce) CYCLE
      !
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO jt = 1,ntracer
        DO jk = 1,nlev
          DO jc = jcs, jce
            !
            IF (echam_phy_config(jg)%iqneg_d2p /= 0) THEN
                IF (pt_prog_new_rcf% tracer(jc,jk,jb,jt) < 0.0_wp) THEN
#ifndef _OPENACC
                  IF (echam_phy_config(jg)%iqneg_d2p == 1 .OR. echam_phy_config(jg)%iqneg_d2p == 3) THEN
                     CALL print_value('d2p:grid   index jg',jg)
                     CALL print_value('d2p:tracer index jt',jt)
                     CALL print_value('d2p:level  index jk',jk)
                     CALL print_value('d2p:pressure   [Pa]',field% pfull(jc,jk,jb))
                     CALL print_value('d2p:longitude [deg]',field% clon(jc,jb)*rad2deg)
                     CALL print_value('d2p:latitude  [deg]',field% clat(jc,jb)*rad2deg)
                     CALL print_value('d2p:pt_prog_new_rcf%tracer',pt_prog_new_rcf% tracer(jc,jk,jb,jt))
                  END IF
#endif
                  IF (echam_phy_config(jg)%iqneg_d2p == 2 .OR. echam_phy_config(jg)%iqneg_d2p == 3) THEN
                     pt_prog_new_rcf% tracer(jc,jk,jb,jt) = 0.0_wp
                  END IF
               END IF
            END IF
            !
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
      !
    END DO
!$OMP END DO
!$OMP END PARALLEL
    !
    !=====================================================================================

    !=====================================================================================
    !
    ! (2) Diagnostics
    !
    ! - pt_diag%tempv    = field%tv
    ! - pt_diag%temp     = field%ta
    ! - pt_diag%pres_sfc                surface pressure filtered to remove sound waves, see diagnose_pres_temp
    ! - pt_diag%pres_ifc = field%phalf  hydrostatic pressure at layer interface
    ! - pt_diag%pres     = field%pfull  hydrostatic pressure at layer midpoint = SQRT(upper pres_ifc * lower pres_ifc)
    ! - pt_diag%dpres_mc                pressure thickness of layer
    !
    CALL diagnose_pres_temp( p_metrics                ,&
      &                      pt_prog_new              ,&
      &                      pt_prog_new_rcf          ,&
      &                      pt_diag                  ,&
      &                      patch                    ,&
      &                      opt_calc_temp=.TRUE.     ,&
      &                      opt_calc_pres=.TRUE.     ,&
      &                      opt_rlend=min_rlcell_int ,& 
      &                      opt_lconstgrav=upatmo_config(jg)%echam_phy%l_constgrav )

    IF (ltimer) CALL timer_stop(timer_d2p_prep)

    ! - pt_diag%u
    ! - pt_diag%v
    IF (ltimer) CALL timer_start(timer_d2p_sync)
    CALL sync_patch_array( SYNC_E, patch, pt_prog_new%vn )
    IF (ltimer) CALL timer_stop(timer_d2p_sync)
    !
    IF (ltimer) CALL timer_start(timer_d2p_prep)

    CALL rbf_vec_interpol_cell( pt_prog_new%vn   ,&! in
      &                         patch            ,&! in
      &                         pt_int_state     ,&! in
      &                         pt_diag%u        ,&! out
      &                         pt_diag%v        ,&! out
      &                         opt_rlstart=rls_c  ,&! in
      &                         opt_rlend  =rle_c  ) ! in

    IF (ltimer) CALL timer_stop(timer_d2p_prep)
    
    !
    ! Now the new prognostic and diagnostic state variables of the dynamical core
    ! are ready to be used in the phyiscs.
    !
    !=====================================================================================


    !=====================================================================================
    !
    ! (3) Copy the new prognostic state and the related diagnostics from the
    !     dynamics state variables to the physics state variables
    !

    ! Loop over cells
    IF (ltimer) CALL timer_start(timer_d2p_couple)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs_c,jbe_c
      !
      CALL get_indices_c(patch, jb,jbs_c,jbe_c, jcs,jce, rls_c,rle_c)
      IF (jcs>jce) CYCLE
      !
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1,nlev
        DO jc = jcs, jce
          !
          ! Fill the time dependent physics state variables, which are used by echam:
          !
          ! density
          field%       rho(jc,jk,jb) = pt_prog_new% rho(jc,jk,jb)
          !
          ! air mass
          field%      mair(jc,jk,jb) = pt_prog_new% rho(jc,jk,jb) &
            &                         *field%        dz(jc,jk,jb) &
            ! (deep-atmosphere modification factor for cell volume)
            &                         *deepatmo_vol(jk)
          !
          ! H2O mass (vap+liq+ice)
          field%      mh2o(jc,jk,jb) = ( pt_prog_new_rcf% tracer(jc,jk,jb,iqv)  &
            &                           +pt_prog_new_rcf% tracer(jc,jk,jb,iqc)  &
            &                           +pt_prog_new_rcf% tracer(jc,jk,jb,iqi)) &
            &                         *field% mair(jc,jk,jb)
          !
          ! dry air mass
          field%      mdry(jc,jk,jb) = field% mair(jc,jk,jb) &
            &                         -field% mh2o(jc,jk,jb)
          !
          ! cloud water+ice
          IF (echam_phy_config(jg)%ldrymoist) THEN
            field%    mref(jc,jk,jb) = field% mdry(jc,jk,jb)
            field%    xref(jc,jk,jb) = field% mair(jc,jk,jb) &
              &                       /field% mdry(jc,jk,jb)
          ELSE
            field%    mref(jc,jk,jb) = field% mair(jc,jk,jb)
            field%    xref(jc,jk,jb) = 1._wp
          END IF
          !
          ! vertical velocity in p-system
          ! (deep-atmosphere modification of 'grav' is assumed to be negligible here)
          field%     omega(jc,jk,jb) = -0.5_wp                                                  &
            &                         * (pt_prog_new% w(jc,jk,jb) + pt_prog_new% w(jc,jk+1,jb)) &
            &                         *  pt_prog_new% rho(jc,jk,jb) * grav
          !
          ! Now reset the physics tendencies before entering the physics:
          tend% ua_phy(jc,jk,jb)     = 0.0_wp
          tend% va_phy(jc,jk,jb)     = 0.0_wp
          tend% ta_phy(jc,jk,jb)     = 0.0_wp
          !
        END DO
      END DO
      !$ACC END PARALLEL
      !
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1,nlev+1
        DO jc = jcs, jce
          !
          field%     wa(jc,jk,jb) = pt_prog_new% w(jc,jk,jb)          !
        END DO
      END DO
      !$ACC END PARALLEL
      !
    END DO ! jb
!$OMP END DO
!$OMP END PARALLEL

    CALL sync_patch_array( SYNC_C, patch, field%wa )

!$OMP PARALLEL
!$OMP DO PRIVATE(jt,jb,jk,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs_c,jbe_c
      !
      CALL get_indices_c(patch, jb,jbs_c,jbe_c, jcs,jce, rls_c,rle_c)
      IF (jcs>jce) CYCLE
      !
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO jt = 1,ntracer
        DO jk = 1,nlev
          DO jc = jcs, jce
            !
            ! Tracer mass
            !
            field%      mtrc(jc,jk,jb,jt)  = pt_prog_new_rcf% tracer(jc,jk,jb,jt) &
                 &                            *field%           mair  (jc,jk,jb)
            !
            ! Tracer mass fraction
            field%      qtrc(jc,jk,jb,jt)  = pt_prog_new_rcf% tracer(jc,jk,jb,jt) &
              &                             *field% xref(jc,jk,jb)

            ! Tracer transport tendency
            tend% qtrc_dyn(jc,jk,jb,jt)  = pt_diag% ddt_tracer_adv(jc,jk,jb,jt) &
              &                           *field% xref(jc,jk,jb)
            !
            ! Initialize the total tendencies, to be computed later:
            tend% qtrc    (jc,jk,jb,jt)  = 0.0_wp
            !
            ! Now reset the physics tendencies before entering the physics
            tend% qtrc_phy(jc,jk,jb,jt)  = 0.0_wp
            !
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
      !
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

    CALL echam_phy_bcs( patch        ,&! in
      &                 datetime_old ,&! in
      &                 dt_loc       ) ! out

    IF (ltimer) CALL timer_stop(timer_echam_bcs)
    !
    !=====================================================================================
    !
    ! (4) Call echam physics and compute the total physics tendencies.
    !     This includes the atmospheric processes (proper ECHAM) and
    !     the land processes, which are vertically implicitly coupled
    !     to the parameterization of vertical turbulent fluxes.
    !
#ifndef __NO_JSBACH__
    IF (echam_phy_config(jg)%ljsb) THEN
      CALL jsbach_start_timestep(jg, datetime_old, dt_loc)
    END IF
#endif

    IF (ltimer) CALL timer_start(timer_echam_phy)

    ! Like in ECHAM, the subroutine *echam_phy_main* has direct access to the memory
    ! buffers prm_field and prm_tend. 

    CALL echam_phy_main( patch            & ! in
      &                 ,datetime_old     & ! in
      &                 ,dt_loc           ) ! in

    IF (ltimer) CALL timer_stop(timer_echam_phy)

#ifndef __NO_JSBACH__
    IF (echam_phy_config(jg)%ljsb) THEN
      CALL jsbach_finish_timestep(jg, datetime_old, dt_loc)
    END IF
#endif

    CALL deallocateDatetime(datetime_old)
    !
    !=====================================================================================
    !
    ! (5) Couple to ocean surface if an ocean is present and this is a coupling time step.
    !
    !
#ifdef YAC_coupling
    IF ( is_coupled_run() ) THEN
#if defined( _OPENACC )
      CALL warning('GPU:interface_echam_ocean','GPU host synchronization should be removed when port is done!')
      CALL gpu_update_var_list('prm_field_D', .false., jg)
      CALL gpu_update_var_list('prm_tend_D', .false., jg)
#endif

      IF (ltimer) CALL timer_start(timer_coupling)

      CALL interface_echam_ocean( patch , pt_diag )

      IF (ltimer) CALL timer_stop(timer_coupling)

#if defined( _OPENACC )
      CALL warning('GPU:interface_echam_ocean','GPU device synchronization should be removed when port is done!')
      CALL gpu_update_var_list('prm_field_D', .true., jg)
      CALL gpu_update_var_list('prm_tend_D', .true., jg)
#endif
    END IF
#endif
    !
    !=====================================================================================

    IF (ltimer) CALL timer_start(timer_phy2dyn)

    !=====================================================================================
    !
    ! (6) Convert physics tendencies to dynamics tendencies
    !
    IF (ltimer) CALL timer_start(timer_p2d_prep)

    !     (a) diagnose again temperature, which is provisionally updated in phyiscs,
    !         from the "new" state after dynamics, so that the temperature field
    !         can be used for updating the model state
    !
    CALL diagnose_pres_temp( p_metrics                ,&
      &                      pt_prog_new              ,&
      &                      pt_prog_new_rcf          ,&
      &                      pt_diag                  ,&
      &                      patch                    ,&
      &                      opt_calc_temp=.TRUE.     ,&
      &                      opt_calc_pres=.FALSE.    ,&
      &                      opt_rlend=min_rlcell_int ,& 
      &                      opt_lconstgrav=upatmo_config(jg)%echam_phy%l_constgrav )

    !     (b) (du/dt|phy, dv/dt|phy) --> dvn/dt|phy
    !
    ALLOCATE(zdudt(nproma,nlev,patch%nblks_c), &
      &      zdvdt(nproma,nlev,patch%nblks_c), &
      &      stat=return_status)
    IF (return_status > 0) THEN
      CALL finish (module_name//method_name, 'ALLOCATE(zdudt,zdvdt)')
    END IF
    !$ACC DATA CREATE( zdudt, zdvdt )

    DO jb = 1, patch%nblks_c
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = 1, nproma
          zdudt(jc,jk,jb) = 0.0_wp
          zdvdt(jc,jk,jb) = 0.0_wp
        END DO
      END DO
      !$ACC END PARALLEL
    END DO

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs_c,jbe_c
      !
      CALL get_indices_c(patch, jb,jbs_c,jbe_c, jcs,jce, rls_c,rle_c)
      IF (jcs>jce) CYCLE
      !
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = jcs, jce
          zdudt(jc,jk,jb) = tend% ua_phy(jc,jk,jb)
          zdvdt(jc,jk,jb) = tend% va_phy(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
      !
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1,nlev+1
        DO jc = jcs, jce
          !
           pt_prog_new% w(jc,jk,jb)=field%     wa(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
      !
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL sync_patch_array( SYNC_C, patch, pt_prog_new%w )

    IF (ltimer) CALL timer_stop(timer_p2d_prep)

    IF (ltimer) CALL timer_start(timer_p2d_sync)
    CALL sync_patch_array_mult(SYNC_C, patch, 2, zdudt, zdvdt)
    IF (ltimer) CALL timer_stop(timer_p2d_sync)

    IF (ltimer) CALL timer_start(timer_p2d_prep)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,jes,jee,jcn,jbn,zvn1,zvn2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs_e,jbe_e
      !
      CALL get_indices_e(patch, jb,jbs_e,jbe_e, jes,jee, rls_e,rle_e)
      IF (jes>jee) CYCLE
      !
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE( jcn, jbn, zvn1, zvn2 ) COLLAPSE(2)
      DO jk = 1,nlev
        DO je = jes,jee
          !
          jcn  =   patch%edges%cell_idx(je,jb,1)
          jbn  =   patch%edges%cell_blk(je,jb,1)
          zvn1 =   zdudt(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,1)%v1 &
            &    + zdvdt(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,1)%v2
          !
          jcn  =   patch%edges%cell_idx(je,jb,2)
          jbn  =   patch%edges%cell_blk(je,jb,2)
          zvn2 =   zdudt(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,2)%v1 &
            &    + zdvdt(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,2)%v2
          !
          pt_diag%ddt_vn_phy(je,jk,jb) =   REAL(  pt_int_state%c_lin_e(je,1,jb)*zvn1      &
            &                                   + pt_int_state%c_lin_e(je,2,jb)*zvn2, vp)
          !
        END DO ! je
      END DO ! jk
      !$ACC END PARALLEL
      !
    END DO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA
    DEALLOCATE(zdudt, zdvdt)

    IF (ltimer) CALL timer_stop(timer_p2d_prep)
    !
    !=====================================================================================


    !=====================================================================================
    !
    ! (7) Couple dynamics+transport and physics

    IF (ltimer) CALL timer_start(timer_p2d_couple)
    !
    ! Here the physics forcing is treated as "fast" physics:
    ! - The provisional "new" state is updated with the total phyiscs
    !   tendencies, providing the final "new" state
    ! - The physics forcing that is passed to the dynamical
    !   core must be set to zero

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,jes,jee) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs_e,jbe_e
      !
      CALL get_indices_e(patch, jb,jbs_e,jbe_e, jes,jee, rls_e,rle_e)
      IF (jes>jee) CYCLE
      !
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO je = jes, jee
          !
          ! (1) Velocity
          !
          ! Update with the total physics tendencies
          pt_prog_new%vn    (je,jk,jb) =   pt_prog_new%vn    (je,jk,jb)             &
            &                            + pt_diag%ddt_vn_phy(je,jk,jb) * dt_loc
          !
          ! Set physics forcing to zero so that it is not re-applied in the dynamical core
          pt_diag%ddt_vn_phy(je,jk,jb) = 0._wp
          !
        END DO
      END DO
      !$ACC END PARALLEL
      !
    END DO !jb
!$OMP END DO
!$OMP END PARALLEL

    IF (lart) jt_end = advection_config(jg)%nname

    ! Loop over cells
!$OMP PARALLEL
!$OMP DO PRIVATE(jt,jb,jk,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs_c,jbe_c
      !
      CALL get_indices_c(patch, jb,jbs_c,jbe_c, jcs,jce, rls_c,rle_c)
      IF (jcs>jce) CYCLE
      !
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jt = 1,jt_end
        DO jc = jcs, jce
          field% mtrcvi    (jc,jb,jt) = 0.0_wp
          tend%  mtrcvi_phy(jc,jb,jt) = 0.0_wp
        END DO
      END DO
      !$ACC END PARALLEL
      !
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP SEQ
      DO jk = 1,nlev
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jt = 1,jt_end
          DO jc = jcs, jce
            !
            ! Diagnose the total tendencies
            tend% qtrc      (jc,jk,jb,jt) = tend% qtrc_dyn(jc,jk,jb,jt)  &
              &                            +tend% qtrc_phy(jc,jk,jb,jt)
            !
            ! (2.1) Tracer mixing ratio with respect to reference air mass
            !
            ! tracer mass tendency
            tend% mtrc_phy  (jc,jk,jb,jt) = tend% qtrc_phy(jc,jk,jb,jt) &
              &                            *field% mref(jc,jk,jb)
            !
            ! tracer path tendency
            tend% mtrcvi_phy(jc,   jb,jt) = tend% mtrcvi_phy(jc,   jb,jt) &
              &                            +tend% mtrc_phy  (jc,jk,jb,jt)
            !
            ! If the physics tendency is /= 0 then change the tracer mass
            ! and optionally check and correct negative values.
            IF (tend% mtrc_phy  (jc,jk,jb,jt) /= 0.0_wp) THEN
              !
              ! new tracer mass
              field% mtrc     (jc,jk,jb,jt) = field% mtrc(jc,jk,jb,jt) &
                &                            +tend% mtrc_phy(jc,jk,jb,jt) &
                &                            *dt_loc
              !
              ! Handling of negative tracer mass coming from physics
              !   qtrc as well as other fields are derived from mtrc.
              !   Therefore check mtrc for negative values.
              !
              IF (echam_phy_config(jg)%iqneg_p2d /= 0) THEN
                IF (field% mtrc(jc,jk,jb,jt) < 0.0_wp) THEN
#ifndef _OPENACC
                  IF (echam_phy_config(jg)%iqneg_p2d == 1 .OR. echam_phy_config(jg)%iqneg_p2d == 3) THEN
                    CALL print_value('p2d:grid   index jg',jg)
                    CALL print_value('p2d:tracer index jt',jt)
                    CALL print_value('p2d:level  index jk',jk)
                    CALL print_value('p2d:pressure   [Pa]',field% pfull(jc,jk,jb))
                    CALL print_value('p2d:longitude [deg]',field% clon(jc,jb)*rad2deg)
                    CALL print_value('p2d:latitude  [deg]',field% clat(jc,jb)*rad2deg)
                    CALL print_value('p2d:field%mtrc     ',field% mtrc(jc,jk,jb,jt))
                  END IF
#endif
                  IF (echam_phy_config(jg)%iqneg_p2d == 2 .OR. echam_phy_config(jg)%iqneg_p2d == 3) THEN
                    field% mtrc(jc,jk,jb,jt) = 0.0_wp
                  END IF
                END IF
              END IF
                  !
            END IF
            !
            ! new tracer path
            field% mtrcvi   (jc,   jb,jt) = field% mtrcvi  (jc,   jb,jt) &
            &                            +field% mtrc    (jc,jk,jb,jt)
            !
          END DO  ! jc
        END DO    ! jt
      END DO      ! jk
      !$ACC END PARALLEL
      !
    END DO  ! jb
!$OMP END DO
!$OMP END PARALLEL

    ! Loop over cells
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs_c,jbe_c
      !
      CALL get_indices_c(patch, jb,jbs_c,jbe_c, jcs,jce, rls_c,rle_c)
      IF (jcs>jce) CYCLE
      !
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR
      DO jc = jcs, jce
        ! initialize vertical integrals
        field% mh2ovi(jc,jb) = 0.0_wp
        field% mairvi(jc,jb) = 0.0_wp
        field% mdryvi(jc,jb) = 0.0_wp
        field% mrefvi(jc,jb) = 0.0_wp
        field% cptgzvi(jc,jb) = 0.0_wp
      END DO
      !$ACC END PARALLEL
      !
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP SEQ
      DO jk = 1,nlev
        !$ACC LOOP GANG VECTOR
        DO jc = jcs, jce
          !
          ! new h2o mass
          field% mh2o  (jc,jk,jb) = field% mtrc (jc,jk,jb,iqv) &
            &                      +field% mtrc (jc,jk,jb,iqc) &
            &                      +field% mtrc (jc,jk,jb,iqi)
          !
          IF (echam_phy_config(jg)%ldrymoist) THEN
            !
            ! new air mass
            field%       mair(jc,jk,jb) = field% mref(jc,jk,jb) &
              &                          +field% mh2o(jc,jk,jb)
            !
            ! new density
            field%       rho (jc,jk,jb) = field% mair(jc,jk,jb) &
              &                          /field% dz  (jc,jk,jb) &
              ! (deep-atmosphere modification factor for cell volume)
              &                          /deepatmo_vol(jk)
            !
            ! new density
            pt_prog_new% rho (jc,jk,jb) = field% mair(jc,jk,jb) &
              &                          /field% dz  (jc,jk,jb) &
              ! (deep-atmosphere modification factor for cell volume)
              &                          /deepatmo_vol(jk)
            !
          ELSE
            !
            ! new dry air mass
            field%       mdry(jc,jk,jb) = field% mref(jc,jk,jb) &
              &                          -field% mh2o(jc,jk,jb)
            !              
          END IF
          !
          ! h2o path
          ! DA these guys prevent collapsing loops!!
          ! DA TODO: figure out how to best optimize that
          field% mh2ovi(jc,   jb) = field% mh2ovi(jc,   jb) &
              &                    +field% mh2o  (jc,jk,jb)
          !
          ! air path
          field% mairvi(jc,   jb) = field% mairvi(jc,   jb) &
              &                    +field% mair  (jc,jk,jb)
          !
          ! dry air path
          field% mdryvi(jc,   jb) = field% mdryvi(jc,   jb) &
            &                      +field% mdry  (jc,jk,jb)
          !
          ! reference air path
          field% mrefvi(jc,   jb) = field% mrefvi(jc,   jb) &
            &                      +field% mref  (jc,jk,jb)
          !
          ! reference air path
          field% cptgzvi(jc,  jb) = field% cptgzvi(jc,   jb) &
            &                      +(field%cptgz(jc,jk,jb)*field%rho(jc,jk,jb)*field%dz(jc,jk,jb))
          !
        END DO
      END DO
      !$ACC END PARALLEL
      !
    END DO
!$OMP END DO
!$OMP END PARALLEL

    ! Loop over cells
!$OMP PARALLEL
!$OMP DO PRIVATE(jt,jb,jk,jc,jcs,jce,z_qsum,z_exner) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs_c,jbe_c
      !
      CALL get_indices_c(patch, jb,jbs_c,jbe_c, jcs,jce, rls_c,rle_c)
      IF (jcs>jce) CYCLE
      !
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP SEQ
      DO jt =1,jt_end
        !$ACC LOOP GANG(static:1) VECTOR COLLAPSE(2)
        DO jk = 1,nlev
          DO jc = jcs, jce
            !
            ! new tracer mass fraction with respect to reference air mass
            field%           qtrc   (jc,jk,jb,jt)  = field%  mtrc(jc,jk,jb,jt) &
              &                                     /field%  mref(jc,jk,jb)
            !
            IF (echam_phy_config(jg)%ldrymoist) THEN
              ! in this case mtrc or mair or both may have changed, and
              ! therefore always re-compute the tracer variable
              pt_prog_new_rcf% tracer (jc,jk,jb,jt)  = field%  mtrc(jc,jk,jb,jt) &
                &                                     /field%  mair(jc,jk,jb)
              !
            ELSE IF (echam_phy_config(jg)%l2moment) Then
              ! this is a special case for the 2 moment scheme as field%mtrc
              ! is changed by clipping after calculating changes due to physical
              ! tendencies, so we have to re-compute all tracer variables.
              pt_prog_new_rcf% tracer (jc,jk,jb,jt)  = field% mtrc(jc,jk,jb,jt) &
                &                                     /field%  mair(jc,jk,jb)
            ELSE
              ! in this case mair is unchanged, and the tracer variable
              ! is re-computed only if the mtrc is changed due to a physical
              ! tendency so that changes of only numerical origin are avoided
              IF (tend% mtrc_phy  (jc,jk,jb,jt) /= 0.0_wp) THEN
                pt_prog_new_rcf% tracer (jc,jk,jb,jt)  = field%  mtrc(jc,jk,jb,jt) &
                   &                                     /field%  mair(jc,jk,jb)
              END IF
              !
            END IF
            !
          END DO
        END DO
      END DO

      IF (lart) THEN
        !$ACC LOOP SEQ
        DO jt = jt_end+1,ntracer
          !$ACC LOOP GANG(static:1) VECTOR COLLAPSE(2)
          DO jk = 1,nlev
            DO jc = jcs, jce
              pt_prog_new_rcf% tracer(jc,jk,jb,jt) = prm_field(jg)%qtrc(jc,jk,jb,jt)  +prm_tend(jg)%qtrc_phy(jc,jk,jb,jt)*dt_loc

            ENDDO
          ENDDO
        ENDDO       
      ENDIF

      !$ACC LOOP GANG(static:1) VECTOR COLLAPSE(2)
      DO jk = 1,nlev
        DO jc = jcs, jce
          !
          ! Re-compute new exner, tempv and theta_v only if physics tendencies
          ! of temperature and water traces are non-zero in order to avoid
          ! purely numerical changes.
          IF (      tend% ta_phy   (jc,jk,jb)     /= 0.0_wp &
            & .OR. tend% mtrc_phy (jc,jk,jb,iqv) /= 0.0_wp &
            & .OR. tend% mtrc_phy (jc,jk,jb,iqc) /= 0.0_wp &
            & .OR. tend% mtrc_phy (jc,jk,jb,iqi) /= 0.0_wp ) THEN
            !
            ! (3) Exner function and virtual potential temperature
            !
            ! (a) Update T, then compute Temp_v, Exner and Theta_v
            !
            pt_diag% temp (jc,jk,jb) =   pt_diag% temp  (jc,jk,jb)             &
              &                        + tend%    ta_phy(jc,jk,jb) * dt_loc
            !
            z_qsum = SUM(pt_prog_new_rcf%tracer(jc,jk,jb,advection_config(jg)%trHydroMass%list))
            !
            pt_diag% tempv(jc,jk,jb) =   pt_diag%temp(jc,jk,jb)                                            &
              &                       * ( 1._wp +  vtmpc1 * pt_prog_new_rcf% tracer(jc,jk,jb,iqv) - z_qsum)
            !
            ! Save provisional "new" exner from the slow-physics-forced dynamics
            z_exner = pt_prog_new% exner(jc,jk,jb)
            !
            ! Compute final new exner
            pt_prog_new% exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd/p0ref * pt_prog_new% rho(jc,jk,jb) * pt_diag% tempv(jc,jk,jb)))
            !
            ! Add exner change from fast physics to exner_pr in order to avoid unphysical sound wave generation
            pt_diag% exner_pr(jc,jk,jb)  = pt_diag% exner_pr(jc,jk,jb) + pt_prog_new% exner(jc,jk,jb) - z_exner
            !
!!$            ! (b) Update Exner, then compute Temp_v
!!$            !
!!$            pt_prog_new%exner(jc,jk,jb) = pt_prog_new% exner(jc,jk,jb)                 &
!!$              &                         + pt_diag% ddt_exner_phy(jc,jk,jb) * dt_loc
!!$            pt_diag%exner_old(jc,jk,jb) = pt_diag% exner_old(jc,jk,jb)                 &
!!$              &                         + pt_diag% ddt_exner_phy(jc,jk,jb) * dt_loc
!!$            !
!!$            pt_diag%tempv(jc,jk,jb) = EXP(LOG(pt_prog_new%exner(jc,jk,jb)/rd_o_cpd)) &
!!$              &                     / (pt_prog_new%rho(jc,jk,jb)*rd/p0ref)
!!$            !
            !
            ! (a) and (b) Compute Theta_v
            !
            pt_prog_new% theta_v(jc,jk,  jb) = pt_diag% tempv(jc,jk,jb) / pt_prog_new% exner(jc,jk,jb)
            !
          END IF
          !
          ! Set physics forcing to zero so that it is not re-applied in the dynamical core
          pt_diag% ddt_exner_phy(jc,jk,jb) = 0._wp
          !
          ! Additionally use this loop also to set the dynamical exner increment to zero.
          ! (It is accumulated over one advective time step in solve_nh)
          pt_diag% exner_dyn_incr(jc,jk,jb) = 0._wp
          !
        END DO
      END DO
      !$ACC END PARALLEL
      !
    END DO !jb
!$OMP END DO
!$OMP END PARALLEL

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

    !$ACC END DATA

    IF (ltimer) CALL timer_stop(timer_p2d_sync)
    !
    !=====================================================================================

    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) CALL timer_stop(timer_phy2dyn)

    !=====================================================================================
    !
    ! Now the final new state (pt_prog_new/pt_prog_new_rcf) is ready.
    !
    !=====================================================================================

  END SUBROUTINE interface_iconam_echam
  !----------------------------------------------------------------------------

END MODULE mo_interface_iconam_echam
