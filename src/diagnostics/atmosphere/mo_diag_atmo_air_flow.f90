!>
!! Atmospheric diagnostics around motion of air
!!
!! This module contains subroutines for the computation of: 
!! * Zonal component of relative vorticity
!! * Meridional component of relative vorticity
!!
!! @par Revision History
!! Initial revision by Sebastian Borchert, DWD (2020-05-29)
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!
!----------------------------
#include "omp_definitions.inc"
!----------------------------
!
MODULE mo_diag_atmo_air_flow

  USE mo_kind,                      ONLY: wp, vp
  USE mo_exception,                 ONLY: finish, message
  USE mo_impl_constants,            ONLY: MAX_CHAR_LENGTH, SUCCESS, &
    &                                     min_rledge_int,           &
    &                                     start_prog_cells, end_prog_cells
  USE mo_upatmo_impl_const,         ONLY: idamtr
  USE mo_model_domain,              ONLY: t_patch
  USE mo_intp_data_strc,            ONLY: t_int_state
  USE mo_nonhydro_types,            ONLY: t_nh_prog, t_nh_metrics
  USE mo_parallel_config,           ONLY: nproma
  USE mo_loopindices,               ONLY: get_indices_c, get_indices_e
  USE mo_icon_interpolation_scalar, ONLY: cells2edges_scalar, cells2verts_scalar
  USE mo_intp_rbf,                  ONLY: rbf_vec_interpol_edge
  USE mo_math_gradients,            ONLY: grad_fd_tang
  USE mo_sync,                      ONLY: SYNC_E, sync_patch_array
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_opt_diag_atmo_vor

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hor_comps_of_rel_vorticity

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_diag_atmo_air_flow'

  !-----------------------------------------------------
  !     Edge-normal component of relative vorticity
  !-----------------------------------------------------

  REAL(wp), ALLOCATABLE :: vor_n(:,:,:)
  INTEGER               :: counter = 0

CONTAINS

  !>
  !! Horizontal components of relative vorticity
  !!
  !! Compute zonal and meridional components of relative vorticity.
  !!
  SUBROUTINE hor_comps_of_rel_vorticity( p_patch,      & !inout
    &                                    p_int_state,  & !in
    &                                    p_metrics,    & !in
    &                                    p_prog,       & !in
    &                                    both_comps,   & !in
    &                                    opt_vor_u,    & !optinout
    &                                    opt_vor_v,    & !optinout
    &                                    opt_startlev, & !optin
    &                                    opt_endlev,   & !optin
    &                                    opt_timer,    & !optin
    &                                    opt_verbose   ) !optin

    ! in/out variables

    TYPE(t_patch),                INTENT(INOUT) :: p_patch          !< horizontal grid
    TYPE(t_int_state),            INTENT(IN)    :: p_int_state      !< for horizontal interpolation
    TYPE(t_nh_metrics),           INTENT(IN)    :: p_metrics        !< vertical grid
    TYPE(t_nh_prog),              INTENT(IN)    :: p_prog           !< prognostic variables
    LOGICAL,                      INTENT(IN)    :: both_comps       !< u- and v-component requested for output?
    REAL(wp),           OPTIONAL, INTENT(INOUT) :: opt_vor_u(:,:,:) !< zonal component of relative vorticity
                                                                    !< (nproma, nlev, nblks_c) [s-1]
    REAL(wp),           OPTIONAL, INTENT(INOUT) :: opt_vor_v(:,:,:) !< meridional component of relative vorticity
                                                                    !< (nproma, nlev, nblks_c) [s-1]
    INTEGER,            OPTIONAL, INTENT(IN)    :: opt_startlev     !< start index for vertical loops
    INTEGER,            OPTIONAL, INTENT(IN)    :: opt_endlev       !< end index for vertical loops
    LOGICAL,            OPTIONAL, INTENT(IN)    :: opt_timer        !< switch for timer monitoring
    LOGICAL,            OPTIONAL, INTENT(IN)    :: opt_verbose      !< switch for message output

    ! local variables

    INTEGER  :: nlev, nblks_e
    INTEGER  :: reset
    INTEGER  :: istat
    LOGICAL  :: present_vor_u, present_vor_v, allocated_vor_n, timer, verbose

    LOGICAL, PARAMETER :: dealloc_vor_n_reg = .TRUE.  ! deallocate vor_n regularly?

    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':hor_comps_of_rel_vorticity'

    !------------------------------------------------

    ! timer monitoring
    IF (PRESENT(opt_timer)) THEN
      timer = opt_timer
    ELSE
      timer = .FALSE.
    ENDIF

    ! transmit messages if required
    IF (PRESENT(opt_verbose)) THEN
      verbose = opt_verbose
    ELSE
      verbose = .FALSE.
    ENDIF

    ! switch on timer monitoring if desired
    IF(timer) CALL timer_start(timer_opt_diag_atmo_vor)

    IF (verbose) CALL message(TRIM(routine), &
      & 'Computation of horizontal vorticity component(s) started')

    present_vor_u   = PRESENT(opt_vor_u)
    present_vor_v   = PRESENT(opt_vor_v)
    allocated_vor_n = ALLOCATED(vor_n)
    
    IF (.NOT. (present_vor_u .OR. present_vor_v)) THEN
      CALL finish(TRIM(routine), 'At least either of opt_vor_u and opt_vor_v should be present')
    ELSEIF (.NOT. ANY((/0, 1, 2/) == counter)) THEN
      CALL finish(TRIM(routine), 'Invalid counter')
    ENDIF

    IF (.NOT. allocated_vor_n) THEN
    
      ! number of vertical grid levels
      nlev = p_patch%nlev
      
      ! number of blocks
      nblks_e = p_patch%nblks_e

      ALLOCATE(vor_n(nproma,nlev,nblks_e), STAT=istat)
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Allocation of vor_n failed')        

    ENDIF

    IF (counter == 0 .OR. .NOT. allocated_vor_n) THEN

      ! compute edge-normal component of relative vorticity

      IF (verbose) CALL message(TRIM(routine), &
        & 'Update edge-normal component of relative vorticity')

      CALL nor_comp_of_rel_vorticity( p_patch     = p_patch,     & !inout
        &                             p_int_state = p_int_state, & !in
        &                             p_metrics   = p_metrics,   & !in
        &                             p_prog      = p_prog,      & !in
        &                             omega_n     = vor_n        ) !inout
      
    ENDIF

    IF (present_vor_u) THEN

      ! interpolate zonal components from edge-normal components

      IF (verbose) CALL message(TRIM(routine), &
        & 'Interpolate zonal components of relative vorticity from normal components')

      CALL rbf_vec_interpol_cell_compwise( p_patch      = p_patch,                     & !inout
        &                                  blk          = p_int_state%rbf_vec_blk_c,   & !in
        &                                  idx          = p_int_state%rbf_vec_idx_c,   & !in
        &                                  coeff        = p_int_state%rbf_vec_coeff_c, & !in
        &                                  vec_n        = vor_n,                       & !in
        &                                  vec_geo      = opt_vor_u,                   & !inout
        &                                  component    = 1,                           & !in
        &                                  opt_startlev = opt_startlev,                & !optin
        &                                  opt_endlev   = opt_endlev                   ) !optin

      counter = counter + 1
      
    ENDIF

    IF (present_vor_v) THEN

      ! interpolate meridional components from edge-normal components

      IF (verbose) CALL message(TRIM(routine), &
        & 'Interpolate meridional components of relative vorticity from normal components')

      CALL rbf_vec_interpol_cell_compwise( p_patch      = p_patch,                     & !inout
        &                                  blk          = p_int_state%rbf_vec_blk_c,   & !in
        &                                  idx          = p_int_state%rbf_vec_idx_c,   & !in
        &                                  coeff        = p_int_state%rbf_vec_coeff_c, & !in
        &                                  vec_n        = vor_n,                       & !in
        &                                  vec_geo      = opt_vor_v,                   & !inout
        &                                  component    = 2,                           & !in
        &                                  opt_startlev = opt_startlev,                & !optin
        &                                  opt_endlev   = opt_endlev                   ) !optin

      counter = counter + 1
      
    ENDIF

    ! clean-up

    reset = MERGE(2, 1, both_comps)
    
    IF (counter == reset) THEN
      counter = 0
      IF (dealloc_vor_n_reg) THEN
        DEALLOCATE(vor_n, STAT=istat)
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Deallocation of vor_n failed')        
      ENDIF
    ENDIF

    IF (verbose) CALL message(TRIM(routine), &
      & 'Computation of horizontal vorticity component(s) finished')

    ! switch of timer monitoring
    IF(timer) CALL timer_stop(timer_opt_diag_atmo_vor)

  END SUBROUTINE hor_comps_of_rel_vorticity

  !=====================================================================================================

  !>
  !! Edge-normal component of relative vorticity
  !!
  !! @par Revision History
  !! This is an excerpt of:
  !! * /src/atm_phy_nwp/mo_opt_nwp_diagnostics: compute_field_pv
  !! by Tobias Selz, LMU
  !!
  SUBROUTINE nor_comp_of_rel_vorticity( p_patch,     & !inout
    &                                   p_int_state, & !in
    &                                   p_metrics,   & !in
    &                                   p_prog,      & !in
    &                                   omega_n      ) !inout

    ! in/out variables

    TYPE(t_patch),              INTENT(INOUT) :: p_patch
    TYPE(t_int_state),  TARGET, INTENT(IN)    :: p_int_state
    TYPE(t_nh_metrics),         INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),            INTENT(IN)    :: p_prog
    REAL(wp),                   INTENT(INOUT) :: omega_n(:,:,:) !< edge-normal component of relative vorticity
                                                                !< (nproma, nlev, nblks_e) [s-1]

    ! local variables

    REAL(wp) :: vt     (nproma, p_patch%nlev,   p_patch%nblks_e), &
      &         w_vh   (nproma, p_patch%nlev+1, p_patch%nblks_v), & 
      &         w_eh   (nproma, p_patch%nlev+1, p_patch%nblks_e), &
      &         ddtw_eh(nproma, p_patch%nlev+1, p_patch%nblks_e)
    REAL(wp) :: invvdfac(p_patch%nlev)

    REAL(vp), POINTER :: ddtz(:,:,:), gamma(:,:,:)

    INTEGER  :: ivd1(p_patch%nlev), &
      &         ivd2(p_patch%nlev)

    REAL(wp) :: invgamma
    INTEGER  :: nlev
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk 
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: jb, je, jk

    !------------------------------------------------

    !-------------------------------------------
    !              Description
    !-------------------------------------------

    ! Computation of the edge-normal component of Rot(v), where v denotes the 3d wind vector.
    ! Fortunately, the computation of this component is part of the computation of the potential vorticity. 
    ! So we can simply copy the respective code parts from:
    !
    ! * /src/atm_phy_nwp/mo_opt_nwp_diagnostics: compute_field_pv
    ! 
    ! Please see there for further details.
    ! See also: Tobias Selz, 2019, "Calculation of the potential voriticity on the ICON grid", 
    ! in Reports on ICON, Issue 02, DOI 10.5676/DWD_pub/nwv/icon_002
    !

    !-------------------------------------------
    !               Computation
    !-------------------------------------------
    
    nlev = p_patch%nlev
    
    ddtz  => p_metrics%ddxt_z_full
    
    gamma => p_metrics%ddqz_z_full_e

    ! get vt at edges (p_diag%vt is not up to date)
    CALL rbf_vec_interpol_edge(p_prog%vn, p_patch, p_int_state, vt)
    
    ! interpolate w to edges
    CALL cells2edges_scalar(p_prog%w, p_patch, p_int_state%c_lin_e, w_eh)
    
    ! interpolate w to vertices
    CALL cells2verts_scalar(p_prog%w, p_patch, p_int_state%cells_aw_verts, w_vh)
    
    ! calculate horizontal derivatives of w
    CALL grad_fd_tang (w_vh, p_patch, ddtw_eh)
    
    ! auxiliary indices and factor for vertical derivatives 
    ! of full level variables
    DO jk = 1, nlev
      IF (jk == 1) THEN
        ivd1(jk)     = 1
        ivd2(jk)     = 2
        invvdfac(jk) = 1._wp
      ELSEIF (jk == nlev) THEN
        ivd1(jk)     = nlev - 1
        ivd2(jk)     = nlev
        invvdfac(jk) = 1._wp
      ELSE
        ivd1(jk)     = jk - 1
        ivd2(jk)     = jk + 1
        invvdfac(jk) = 0.5_wp
      END IF
    ENDDO  !jk
    
    rl_start   = 3
    rl_end     = min_rledge_int - 1
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)
    
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx, invgamma) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk
      
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      
      DO jk = 1, nlev
        
        DO je = i_startidx, i_endidx
          invgamma          = 1._wp / REAL(gamma(je,jk,jb), wp)
          ! edge-normal component of relative vorticity
          omega_n(je,jk,jb) = ( invgamma * invvdfac(jk) * (vt(je,ivd1(jk),jb) - vt(je,ivd2(jk),jb))       &
            &               - p_metrics%deepatmo_t1mc(jk,idamtr%t1mc%gradh) *                             &
            &                 ( 0.5_wp * (ddtw_eh(je,jk,jb) + ddtw_eh(je,jk+1,jb))                        &
            &               + invgamma * REAL(ddtz(je,jk,jb), vp) * (w_eh(je,jk+1,jb) - w_eh(je,jk,jb)) ) &
            &               + p_metrics%deepatmo_t1mc(jk,idamtr%t1mc%invr) * vt(je,jk,jb)                 &
            &                 )
        ENDDO  !je
      ENDDO  !jk
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL
    
    CALL sync_patch_array(SYNC_E, p_patch, omega_n)

    ! clean-up
    NULLIFY(ddtz, gamma)
    
  END SUBROUTINE nor_comp_of_rel_vorticity

  !=====================================================================================================

  !>
  !! Interpolate zonal or meridional vector components at cell center
  !!
  !! Interpolate zonal or meridional vector components at cell center
  !! from edge-normal vector components.
  !!
  !! @par Revision History
  !! This is a modified copy of:
  !! * /src/shr_horizontal/mo_intp_rbf: rbf_vec_interpol_cell
  !! (Please note: we cannot place this subroutine there, 
  !! because mo_intp_rbf is already occupied by openACC directives)
  !!
  SUBROUTINE rbf_vec_interpol_cell_compwise( p_patch,      & !inout
    &                                        blk,          & !in
    &                                        idx,          & !in
    &                                        coeff,        & !in
    &                                        vec_n,        & !in
    &                                        vec_geo,      & !inout
    &                                        component,    & !in
    &                                        opt_startlev, & !optin
    &                                        opt_endlev    ) !optin

    ! in/out variables

    TYPE(t_patch),     INTENT(INOUT) :: p_patch        !< horizontal grid
    INTEGER,           INTENT(IN)    :: blk(:,:,:)     !< block array defining the stencil of surrounding edges 
                                                       !< for vector rbf interpolation at each cell center
                                                       !< (rbf_vec_dim_c, nproma, nblks_c) 
    INTEGER,           INTENT(IN)    :: idx(:,:,:)     !< index array defining the stencil of surrounding edges 
                                                       !< for vector rbf interpolation at each cell center
                                                       !< (rbf_vec_dim_c, nproma, nblks_c)
    REAL(wp),          INTENT(IN)    :: coeff(:,:,:,:) !< array containing the coefficients used for
                                                       !< vector rbf interpolation at each cell center
                                                       !< (rbf_vec_dim_c, 2, nproma, nblks_c)
    REAL(wp),          INTENT(IN)    :: vec_n(:,:,:)   !< edge-normal vector component
                                                       !< (nproma, nlev, nblks_e)
    REAL(wp),          INTENT(INOUT) :: vec_geo(:,:,:) !< zonal or meridional vector component at cell center
                                                       !< (nproma, nlev, nblks_c)
    INTEGER,           INTENT(IN)    :: component      !< indicator for geographic component contained in vec_geo:
                                                       !< * 1: zonal component
                                                       !< * 2: meridional component
    INTEGER, OPTIONAL, INTENT(IN)    :: opt_startlev   !< start level for vertical loop
    INTEGER, OPTIONAL, INTENT(IN)    :: opt_endlev     !< end level for vertical loop

    ! local variables

    INTEGER  :: nlev
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk 
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: i_startlev, i_endlev
    INTEGER  :: jb, jc, jk

    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':rbf_vec_interpol_cell_compwise'
    
    !------------------------------------------------

    IF (.NOT. ANY((/1, 2/) == component)) CALL finish(TRIM(routine), 'Invalid component') 

    nlev = p_patch%nlev

    IF (PRESENT(opt_startlev)) THEN
      i_startlev = MAX(1, opt_startlev)
    ELSE
      i_startlev = 1
    ENDIF
    IF (PRESENT(opt_endlev)) THEN
      i_endlev   = MIN(opt_endlev, nlev)
    ELSE
      i_endlev   = nlev
    ENDIF

    ! prognostic domain
    rl_start   = start_prog_cells
    rl_end     = end_prog_cells
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL 
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk
      
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = i_startlev, i_endlev
#else
!CDIR UNROLL=2
      DO jk = i_startlev, i_endlev
        DO jc = i_startidx, i_endidx
#endif
          vec_geo(jc,jk,jb) = &
            & coeff(1,component,jc,jb) * vec_n(idx(1,jc,jb),jk,blk(1,jc,jb)) + &
            & coeff(2,component,jc,jb) * vec_n(idx(2,jc,jb),jk,blk(2,jc,jb)) + &
            & coeff(3,component,jc,jb) * vec_n(idx(3,jc,jb),jk,blk(3,jc,jb)) + &
            & coeff(4,component,jc,jb) * vec_n(idx(4,jc,jb),jk,blk(4,jc,jb)) + &
            & coeff(5,component,jc,jb) * vec_n(idx(5,jc,jb),jk,blk(5,jc,jb)) + &
            & coeff(6,component,jc,jb) * vec_n(idx(6,jc,jb),jk,blk(6,jc,jb)) + &
            & coeff(7,component,jc,jb) * vec_n(idx(7,jc,jb),jk,blk(7,jc,jb)) + &
            & coeff(8,component,jc,jb) * vec_n(idx(8,jc,jb),jk,blk(8,jc,jb)) + &
            & coeff(9,component,jc,jb) * vec_n(idx(9,jc,jb),jk,blk(9,jc,jb))
        ENDDO  !jk/jc
      ENDDO  !jc/jk
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE rbf_vec_interpol_cell_compwise

END MODULE mo_diag_atmo_air_flow

