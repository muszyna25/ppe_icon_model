!>
!! Nudging.
!!
!! This module contains procedures related to the nudging 
!! of the atmospheric state simulated by ICON towards driving data. 
!! Covered nudging types:
!! - Global nudging
!! For the nudging types:
!! - Lateral boundary nudging
!! - Upper boundary nudging 
!! please see:
!! - src/atm_dyn_iconam/mo_nh_nest_utilities: limarea_bdy_nudging
!! The nudging we refer to in this module is unrelated 
!! to the nudging in the context of:
!! - assimilation (assimilation_nml)
!! - large-scale forcing (ls_forcing_nml)
!! 
!! @par Revision History
!! The content of this module follows closely the nudging 
!! for the limited-area mode in:
!! - src/atm_dyn_iconam/mo_nh_nest_utilities
!! Modification for global nudging by Sebastian Borchert, DWD (2018-09)
!!
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

MODULE mo_nudging

  USE, INTRINSIC :: iso_c_binding, ONLY: c_int32_t

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_impl_constants,        ONLY: MAX_CHAR_LENGTH, SUCCESS, &
    &                                 min_rlcell, min_rledge,   &
    &                                 min_rlcell_int
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_physical_constants,    ONLY: rd, cvd_o_rd, p0ref, vtmpc1, &
    &                                 rcpd, grav, rd_o_cpd, rgrav, &
    &                                 dtdz_standardatm
  USE mo_math_constants,        ONLY: one_third, dbl_eps
  USE mo_run_config,            ONLY: iqv, iqc
  USE mo_parallel_config,       ONLY: nproma
  USE mo_limarea_config,        ONLY: latbc_config
  USE mo_nudging_config,        ONLY: t_nudging_config,    &
    &                                 indg_type, indg_var, &
    &                                 ithermdyn_type
  USE mo_upatmo_config,         ONLY: t_upatmo_config, idamtr
  USE mo_time_config,           ONLY: t_time_config
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_state, t_nh_prog, &
    &                                 t_nh_diag, t_nh_metrics
  USE mo_initicon_types,        ONLY: t_init_state, t_pi_atm
  USE mo_async_latbc_types,     ONLY: t_latbc_data
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mtime,                    ONLY: datetime, timedelta,             &
    &                                 deallocateDatetime, newDatetime, &
    &                                 OPERATOR(>=), OPERATOR(*), OPERATOR(+)
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e
  USE mo_async_latbc_utils,     ONLY: update_lin_interpolation
  USE mo_sync_latbc,            ONLY: update_lin_interc
  USE mo_nh_diagnose_pres_temp, ONLY: diagnose_pres_temp, &
    &                                 diag_temp, diag_pres
  USE mo_util_string,           ONLY: int2string, real2string
  USE mo_math_divrot,           ONLY: div
  USE mo_mpi,                   ONLY: my_process_is_stdio,          &
    &                                 get_my_mpi_work_communicator, &
    &                                 process_mpi_stdio_id, p_bcast, p_sum
  USE mo_timer,                 ONLY: timer_start, timer_stop, &
    &                                 timer_global_nudging
  USE mo_io_units,              ONLY: find_next_free_unit
  USE mo_io_config,             ONLY: inextra_2d

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nudging_interface

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nudging'

CONTAINS !..................................................................

  !>
  !! Nudging interface (for global nudging).
  !!
  !! This subroutine is meant as some kind of wrapper 
  !! for the procedures, which make up one nudging cycle 
  !! (except for the procedures related to the I/O of the driving data).
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert, DWD (2018-09)
  !!
  SUBROUTINE nudging_interface( p_patch,          & !in
    &                           p_nh_state,       & !inout
    &                           p_latbc_data,     & !in
    &                           latbc,            & !in
    &                           p_int_state,      & !in
    &                           mtime_datetime,   & !in
    &                           sim_time,         & !in
    &                           time_config,      & !in
    &                           ndyn_substeps,    & !in
    &                           nnew,             & !in
    &                           nnew_rcf,         & !in
    &                           last_latbc_tlev,  & !in
    &                           read_latbc_tlev,  & !in
    &                           upatmo_config,    & !in
    &                           nudging_config    ) !inout

    ! In/out variables
    TYPE(t_patch),            TARGET,  INTENT(IN)    :: p_patch           !< Grid/patch info
    TYPE(t_nh_state),         TARGET,  INTENT(INOUT) :: p_nh_state        !< Prognostic and diagnostic variables etc.
    TYPE(t_init_state),       TARGET,  INTENT(IN)    :: p_latbc_data(:)   !< Two time-level boundary data               
    TYPE(t_latbc_data),       TARGET,  INTENT(IN)    :: latbc             !< Data structure for async latbc prefetching
    TYPE(t_int_state),                 INTENT(IN)    :: p_int_state       !< For operations on the horizontal grid geometry
    TYPE(datetime),           POINTER, INTENT(IN)    :: mtime_datetime    !< Date/time information
    REAL(wp),                          INTENT(IN)    :: sim_time          !< Elapsed simulation time in seconds
    TYPE(t_time_config),               INTENT(IN)    :: time_config       !< Important times of current simulation
    INTEGER,                           INTENT(IN)    :: ndyn_substeps     !< Number of dynamics' substeps
    INTEGER,                           INTENT(IN)    :: nnew, nnew_rcf    !< Time level indices
    INTEGER,                           INTENT(IN)    :: last_latbc_tlev   !< Time level index for p_latbc_data
    INTEGER,                           INTENT(IN)    :: read_latbc_tlev   !< Time level index for p_latbc_data
    TYPE(t_upatmo_config),             INTENT(IN)    :: upatmo_config     !< Upper-atmosphere switches
    TYPE(t_nudging_config),            INTENT(INOUT) :: nudging_config    !< Nudging switches

    ! Local variables
    TYPE(t_pi_atm), POINTER :: p_latbc_old, p_latbc_new
    REAL(wp) :: tsrat, wfac_old, wfac_new
    INTEGER  :: jg
    LOGICAL  :: l_thermdyn, l_hydrostatic, l_qv, l_message
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER :: &
      & routine = modname//":nudging_interface"
    
    !----------------------------------------

    ! Some notes may be in order: 
    !
    ! - Although it comes at the expense of computational efficiency   
    !   (due to the interposition of this subroutine 
    !   before the actual call of the nudging procedures), 
    !   we introduced this interface, to reduce the optical overload 
    !   of 'src/atm_dyn_iconam/mo_nh_stepping'
    !
    ! - In contrast to 'nudging_config', 'latbc_config' became no input, 
    !   because subroutines such as 'update_lin_interpolation' change 'latbc_config' 
    !   and access it via a USE-binding in their modules. 
    !   So there is no reason for us, to do differently here
    !
    ! - This subroutine is called in 'mo_nh_stepping: integrate_nh' 
    !   only if global nudging is switched on (NOT nudging in general!) 
    !   and if we are on the primary domain ('jg = 1')
    !
    ! - The content of this subroutine follows closely the infrastructure 
    !   around the call of 'limarea_bdy_nudging' in 'mo_nh_stepping: integrate_nh'
    !
    ! - However, there is a lot of nudging infrastructure in 'mo_nh_stepping', 
    !   which cannot be included in this interface: 
    !   * Initialization and finalization of the nudging driving data processing 
    !   * Triggering of the read-in and postprocessing of the driving data
    !   In most cases a modification of this infrastructure for our purposes 
    !   is relatively moderate (search for "l_global_nudging"). 
    !   The prefetching of boundary data in case of asynchronous read-in 
    !   (num_prefetch_proc > 0) is triggered, if 'latbc_config%itype_latbc > 0'. 
    !   This should also work for global nudging, since 'latbc_config%itype_latbc = 1' 
    !   (i.e. time-dependent driving data) is enforced 
    !   in 'src/namelists/mo_nudging_nml: check_nudging'.

    ! Domain index
    jg = p_patch%id

    ! Better check
    IF (jg /= 1) RETURN

    !...............................................................
    !                       Global nudging
    !...............................................................

    IF (nudging_config%ltype(indg_type%globn)) THEN 

      IF (nudging_config%ltimer) CALL timer_start(timer_global_nudging)

      l_message = nudging_config%lmessage
      IF(l_message) CALL message(TRIM(routine), 'Start global nudging.')

      !---------------------------------------------------------------
      !                        Preparation
      !---------------------------------------------------------------

      IF (nudging_config%lasync) THEN
        !
        ! Asynchronous read-in of driving data:
        !
        ! The following subroutine updates the weights: 
        ! * latbc_config%lc1 
        ! * latbc_config%lc2 
        ! for interpolation of the driving data in time 
        CALL update_lin_interpolation(latbc, mtime_datetime)
        ! Set pointer to past and future state of driving data,
        ! from which their current state is estimated by 
        ! linear interpolation in time
        p_latbc_old => latbc%latbc_data( latbc%prev_latbc_tlev() )%atm
        p_latbc_new => latbc%latbc_data( latbc%new_latbc_tlev    )%atm
      ELSE
        !
        ! Synchronous read-in of driving data:
        !
        CALL update_lin_interc(mtime_datetime)
        p_latbc_old => p_latbc_data( last_latbc_tlev )%atm
        p_latbc_new => p_latbc_data( read_latbc_tlev )%atm
      ENDIF

      ! Get weights for time interpolation of driving data
      wfac_old = latbc_config%lc1
      wfac_new = latbc_config%lc2
      
      !---------------------------------------------------------------
      !                 Update diagnostic variables
      !---------------------------------------------------------------
      
      ! If we nudge the hydrostatic thermodynamic variables, 
      ! the diagnostic variables pressure and (virtual) temperature are required 
      ! (for the levels 'nudging_config%ilev_start' to 'nudging_config%ilev_end')

      l_thermdyn    = nudging_config%lvar(indg_var%thermdyn) 
      l_hydrostatic = nudging_config%thermdyn_type == ithermdyn_type%hydrostatic
      l_qv          = nudging_config%lvar(indg_var%qv)
      
      IF (l_thermdyn .AND. l_hydrostatic) THEN
        ! ('diagnose_pres_temp' updates the prognostic, halo and lateral boundary cells)
        CALL diagnose_pres_temp( p_metrics         = p_nh_state%metrics,           & !in
          &                      pt_prog           = p_nh_state%prog( nnew     ),  & !in 
          &                      pt_prog_rcf       = p_nh_state%prog( nnew_rcf ),  & !in
          &                      pt_diag           = p_nh_state%diag,              & !out
          &                      pt_patch          = p_patch,                      & !in
          &                      opt_calc_temp     = .TRUE.,                       & !optin
          &                      opt_calc_pres     = .TRUE.,                       & !optin 
          &                      opt_calc_temp_ifc = .FALSE.,                      & !optin
          &                      opt_slev          = nudging_config%ilev_start,    & !optin
          &                      opt_lconstgrav    = upatmo_config%dyn%l_constgrav ) !optin        
      ENDIF  !Update press and temp(v)?

      !---------------------------------------------------------------
      !                        Apply nudging
      !---------------------------------------------------------------
      
      ! Dynamics-physics time step ratio
      tsrat = REAL(ndyn_substeps, wp)

      CALL global_nudging( p_patch        = p_patch,                     & !in
        &                  p_prog         = p_nh_state%prog( nnew     ), & !inout
        &                  p_prog_rcf     = p_nh_state%prog( nnew_rcf ), & !inout
        &                  p_metrics      = p_nh_state%metrics,          & !in
        &                  p_diag         = p_nh_state%diag,             & !in
        &                  p_latbc_old    = p_latbc_old,                 & !in
        &                  p_latbc_new    = p_latbc_new,                 & !in
        &                  tsrat          = tsrat,                       & !in
        &                  wfac_old       = wfac_old,                    & !in
        &                  wfac_new       = wfac_new,                    & !in
        &                  nudging_config = nudging_config               ) !in

      !---------------------------------------------------------------
      !              Analysis of nudging success/impact
      !---------------------------------------------------------------

      IF (nudging_config%ldiagnose .AND. (.NOT. nudging_config%diag_kit%lfileclosed)) THEN

        ! The analysis takes place every 'idiagnose * dtime'
        IF ( nudging_config%diag_kit%ncount == 0 .OR.                           &
          &  MOD(nudging_config%diag_kit%ncount, nudging_config%idiagnose) == 0 ) THEN

          ! The analysis makes use of 'p_diag%tempv'
          IF ((l_thermdyn .OR. l_qv)) THEN
            CALL diagnose_pres_temp( p_metrics         = p_nh_state%metrics,           & !in
              &                      pt_prog           = p_nh_state%prog( nnew     ),  & !in 
              &                      pt_prog_rcf       = p_nh_state%prog( nnew_rcf ),  & !in
              &                      pt_diag           = p_nh_state%diag,              & !out
              &                      pt_patch          = p_patch,                      & !in
              &                      opt_calc_temp     = .TRUE.,                       & !optin
              &                      opt_calc_pres     = .FALSE.,                      & !optin 
              &                      opt_calc_temp_ifc = .FALSE.,                      & !optin
              &                      opt_slev          = MAX(1, p_patch%nlev - 2),     & !optin
              &                      opt_lconstgrav    = upatmo_config%dyn%l_constgrav ) !optin
          ENDIF  !IF (l_thermdyn ...)

          CALL nudging_diagnostics( p_patch        = p_patch,                 & !in
            &                       p_prog         = p_nh_state%prog( nnew ), & !in
            &                       p_metrics      = p_nh_state%metrics,      & !in
            &                       p_diag         = p_nh_state%diag,         & !in
            &                       p_int_state    = p_int_state,             & !in
            &                       p_latbc_old    = p_latbc_old,             & !in
            &                       p_latbc_new    = p_latbc_new,             & !in
            &                       mtime_datetime = mtime_datetime,          & !in
            &                       sim_time       = sim_time,                & !in
            &                       time_config    = time_config,             & !in
            &                       wfac_old       = wfac_old,                & !in
            &                       wfac_new       = wfac_new,                & !in
            &                       l_message      = l_message,               & !in
            &                       nudging_config = nudging_config           ) !inout

        ENDIF  !IF (nudging_config%diag_kit%ncount == 0 .OR. ... )

        ! The calling counter has to update every call of 'nudging_interface' (i.e. every 'dtime')
        nudging_config%diag_kit%ncount = nudging_config%diag_kit%ncount + 1
        
      ENDIF  !IF (nudging_config%ldiagnose .AND. (.NOT. nudging_config%diag_kit%lfileclosed))

      !---------------------------------------------------------------
      !                         Clean-up
      !---------------------------------------------------------------

      p_latbc_old => NULL()
      p_latbc_new => NULL()

      IF(l_message) CALL message(TRIM(routine), 'End of global nudging.')
      
      IF (nudging_config%ltimer) CALL timer_stop(timer_global_nudging)

    ENDIF  !IF (nudging_config%ltype(indg_type%globn))
    
  END SUBROUTINE nudging_interface

!---------------------------------------------------------------------------

  !>
  !! This routine executes global nudging.
  !!
  !! @par Revision History
  !! This subroutine is a conceptual copy of:
  !! - 'src/atm_dyn_iconam/mo_nh_nest_utilities: limarea_bdy_nudging'
  !! Modifications for global nudging by Sebastian Borchert, DWD (2018-09)
  !!
  SUBROUTINE global_nudging( p_patch,       & !in
    &                        p_prog,        & !inout
    &                        p_prog_rcf,    & !inout
    &                        p_metrics,     & !in
    &                        p_diag,        & !in
    &                        p_latbc_old,   & !in
    &                        p_latbc_new,   & !in
    &                        tsrat,         & !in
    &                        wfac_old,      & !in
    &                        wfac_new,      & !in
    &                        nudging_config ) !in

    ! In/out variables
    TYPE(t_patch),          INTENT(IN)    :: p_patch
    TYPE(t_nh_prog),        INTENT(INOUT) :: p_prog, p_prog_rcf
    TYPE(t_nh_metrics),     INTENT(IN)    :: p_metrics
    TYPE(t_nh_diag),        INTENT(IN)    :: p_diag
    TYPE(t_pi_atm),         INTENT(IN)    :: p_latbc_old, p_latbc_new  !< Past and future state of driving data
    REAL(wp),               INTENT(IN)    :: tsrat                     !< Ratio between advective and dynamical time step
    REAL(wp),               INTENT(IN)    :: wfac_old, wfac_new        !< Weights for time interpolation of driving data
    TYPE(t_nudging_config), INTENT(IN)    :: nudging_config            !< Nudging switches

    ! Local variables
    REAL(wp), ALLOCATABLE :: qv_tend(:,:,:)
    REAL(wp), ALLOCATABLE :: nudge_coeff_thermdyn(:), &
      &                      nudge_coeff_vn(:),       &
      &                      nudge_coeff_qv(:)
    REAL(wp) :: rho_tend, thv_tend, vn_tend
    REAL(wp) :: pres_tend, temp_tend, tempv_tend, qv_update
    REAL(wp) :: nudge_coeff
    INTEGER  :: jg, jc, je, jb, jk
    INTEGER  :: istart, iend
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: rl_start, rl_end
    INTEGER  :: istat
    LOGICAL  :: l_thermdyn, l_vn, l_qv, l_hydrostatic, l_qv_tend
    REAL(wp), PARAMETER :: rd_o_cvd    = 1._wp / cvd_o_rd
    REAL(wp), PARAMETER :: rd_o_p0ref  = rd / p0ref
    REAL(wp), PARAMETER :: rrd         = 1._wp / rd
    REAL(wp), PARAMETER :: eps_qc      = 1.e-10_wp
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER :: &
      & routine = modname//":global_nudging"

    !----------------------------------------

    !........................................................................................
    ! Some notes:
    ! 
    ! - Global nudging is applied to the prognostic and the halo partitions of the domain, 
    !   so we assume that the prognostic variables 'p_prog' as well as the driving data 
    !   on 'p_latbc_old/new' have been synchronized before they enter this subroutine
    !
    ! - Some diagnostic variables ('p_diag%...') are required, 
    !   so they should have been updated before they enter this subroutine
    ! 
    ! - The following variables are potentially (directly) updated by nudging increments: 
    ! 
    !   Const. driving data         | Time-dependent driving data
    !   ----------------------------------------------------------
    !                 { * rho       | * rho            * press
    !   Thermodynamic {             |             or  
    !                 { * theta_v   | * theta_v        * temp
    !                               |
    !                   * vn        | * vn
    !                               |
    !                               | * qv
    !   ----------------------------------------------------------
    !
    ! - The nudging formula reads:
    !       X(t) = X'(t) + ndyn_substeps * nudge_coeff_X(z) * [ X0(t) - X'(t) ],     (0)
    !                                                         |_______________|
    !                                                                 |
    !                                                                = dX
    !   where X' denotes the value of the variable before the nudging step, 
    !   ndyn_substeps is the number of dynamics time steps per fast physics time step, 
    !   nudge_coeff_X is the variable-specific, height-dependent nudging coefficient, 
    !   and X0 is the value of the variable from the nudging driving data. 
    !   Note: the nudging update X - X' does not depend on the time step, 
    !   so it is not a typical state transition process, 
    !   but rather a state displacement or shift.
    !
    ! - In case of time-dependent driving data, there are two possibilities 
    !   for the provided thermodynamic nudging increments: 
    ! 
    !   Non-hydrostatic variables | Hydrostatic variables
    !   --------------------------------------------------
    !   * drho                    | * dpress
    !                             | 
    !   * dtheta_v                | * dtemp
    !   --------------------------------------------------
    !
    !   The increments drho and dtheta_v are multiplied with the nudging coefficient 
    !   and are directly added to the corresponding prognostic variables of ICON. 
    ! 
    !   The increments dpress and dtemp are transformed into drho and dtheta_v 
    !   using the linear map: 
    !                       drho  = A * dtemp_v + B * dpress,                        (1)
    !                    dtheta_v = C * dtemp_v + D * dpress,                        (2)    
    !   where the factors A, B, C and D are completely determined by the ICON-state 
    !   before the nudging. In order to find expressions for (1) and (2), 
    !   we proceed in the following way:
    !   First, given the definition of the virtual temperature 
    !   for moist (not cloudy) air: 
    !                      temp_v = temp * ( 1 + vtmpc1 * qv ),                      (3)
    !   
    !   we use its differential: 
    !          dtemp_v = dtemp * ( 1 + vtmpc1 * qv ) + temp * vtmpc1 * dqv,          (4)
    !   to compute the virtual temperature increment dtemp_v.
    !   Next, given the equation of state: 
    !                        rho = press / ( rd * temp_v ),                          (5)
    !   we use its differential: 
    !       drho = dpress / ( rd * temp_v ) - press * dtemp_v / ( rd * temp_v**2 )
    !            = [ dpress / rd - rho * dtemp_v ] / temp_v,                         (6)
    !   to compute the density increment drho. 
    !   And finally, given the definition of the virtual potential temperature: 
    !          theta_v = temp_v * ( p0ref / press )**(rd/cpd) = temp_v / exner,      (7)
    !   we use its differential: 
    !     dtheta_v = dtemp_v / exner - temp_v / exner**2 * ( dexner / dpress ) * dpress  
    !              = dtemp_v / exner 
    !              - ( rd / cpd ) * ( temp_v / exner ) * ( 1 / press ) * dpress  
    !              = [ dtemp_v - dpress / ( cpd * rho ) ] / exner,                   (8)
    !   where: 
    !        dexner / dpress = ( rd / cpd ) * ( press / p0ref )**(rd/cp-1) / p0ref
    !                        = ( rd / cpd ) * ( exner / press )                      (9)
    !   and (5) have been used.
    !
    !   After the nudging increments (multiplied with the nudging coefficient) 
    !   have been added to rho and theta_v, the Exner pressure is updated by
    !   the following equivalent to (5):
    !              exner = ( rd * rho * theta_v / p0ref )**(rd/cvd).                 (10)
    !
    ! - Please note that a nudging of the thermodynamic variables 
    !   based on the hydrostatic set of variables (pres and temp) 
    !   does not imply that the nudging increments themselves 
    !   are hydrostatically balanced!
    !........................................................................................

    ! Domain index
    jg = p_patch%id
    
    ! Although this query should have been likely covered before calling this subroutine, 
    ! we might feel better, if we make sure that:
    IF (jg == 1 .AND. nudging_config%ltype(indg_type%globn)) THEN

      ! Which variables should be nudged?
      ! Thermodynamic variables
      l_thermdyn = nudging_config%lvar(indg_var%thermdyn) 
      ! Horizontal wind
      l_vn       = nudging_config%lvar(indg_var%vn)  
      ! Water vapor 
      ! (includes 'ltransport=.true. implicitly)
      l_qv       = nudging_config%lvar(indg_var%qv)

      ! Thermodynamic variables: 
      ! should we nudge the hydrostatic pressure & temperature, 
      ! or the density and virtual potential temperature?
      l_hydrostatic = nudging_config%thermdyn_type == ithermdyn_type%hydrostatic

      ! Do we have to compute water vapor tendencies? 
      l_qv_tend = l_qv .OR. ( l_thermdyn .AND. l_hydrostatic )

      ! Start and end indices for vertical loops
      istart = nudging_config%ilev_start
      iend   = nudging_config%ilev_end

      ! Precompute the nudging coefficient, to reduce the number of multiplications
      ! (one reason, why we do not compute these profiles only once 
      ! and store them for reuse is that 'tsrat = ndyn_substeps' may change during runtime)
      IF (l_thermdyn) THEN 
        ALLOCATE(nudge_coeff_thermdyn(istart:iend), STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Allocation of nudge_coeff_thermdyn failed!')
      ENDIF
      IF (l_vn) THEN 
        ALLOCATE(nudge_coeff_vn(istart:iend), STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Allocation of nudge_coeff_vn failed!')
      ENDIF
      IF (l_qv) THEN 
        ALLOCATE(nudge_coeff_qv(istart:iend), STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Allocation of nudge_coeff_qv failed!')
      ENDIF
      DO jk = istart, iend
        nudge_coeff = tsrat * p_metrics%nudgecoeff_vert(jk)
        IF (l_thermdyn) nudge_coeff_thermdyn(jk) = nudging_config%max_nudge_coeff_thermdyn * nudge_coeff
        IF (l_vn)       nudge_coeff_vn(jk)       = nudging_config%max_nudge_coeff_vn       * nudge_coeff
        IF (l_qv)       nudge_coeff_qv(jk)       = nudging_config%max_nudge_coeff_qv       * nudge_coeff
      ENDDO  !jk

      ! Allocate the field for the water vapor increments
      IF (l_qv_tend) THEN 
        ALLOCATE(qv_tend(nproma,istart:iend,p_patch%nblks_c), STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Allocation of qv_tend failed!')
      ENDIF

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
 
      !...............................................................
      !               Nudging of cell-based variables
      !...............................................................
      
      ! The horizontal loop runs over the entire prognostic domain + halo cells, 
      ! in order to avoid a synchronization afterwards.
      ! (In case of the limited-area mode, the following setting should 
      ! guarantee that the lateral boundary interpolation zone and 
      ! its halo cells are excluded from nudging.)
      rl_start   = grf_bdywidth_c + 1
      rl_end     = min_rlcell
      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)
      
      !---------------------------------------------------------------
      !                  Get water vapor increments
      !---------------------------------------------------------------
      
      ! For efficiency reasons, we put the query for the nudging variables 
      ! outside the jb-loop, but have to pay for it with some code overhead. 
      ! 
      ! Water vapor increments are only required (and allocated!), 
      ! if we nudge water vapor and/or if we nudge the hydrostatic set 
      ! of thermodynamic variables.
      IF (l_qv_tend) THEN
        
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk
          
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          
          DO jk = istart, iend
            DO jc = i_startidx, i_endidx
              
              qv_tend(jc,jk,jb) = wfac_old * p_latbc_old%qv(jc,jk,jb) + wfac_new * p_latbc_new%qv(jc,jk,jb) &
                &               - p_prog_rcf%tracer(jc,jk,jb,iqv)
              
              ! The actual nudging step of qv has to be postponed to after the nudging 
              ! of the thermodynamic variables, because in its computation, 
              ! we prefer to access the unnudged 'p_prog_rcf%tracer(jc,jk,jb,iqv)'.
              
            ENDDO  !jc
          ENDDO  !jk
        ENDDO  !jb
!$OMP END DO
        
      ENDIF  !IF (l_qv_tend)
      
      !---------------------------------------------------------------
      !              Nudging of thermodynamic variables
      !---------------------------------------------------------------
      
      IF (l_thermdyn .AND. l_hydrostatic) THEN
        
        !---------------------------------------------------------------
        !      CASE 1: Nudge hydrostatic pressure and temperature
        !---------------------------------------------------------------
        
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,pres_tend,temp_tend,tempv_tend,thv_tend,rho_tend) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk
          
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          
          DO jk = istart, iend
            DO jc = i_startidx, i_endidx
              
              ! Given increments:
              ! (We nudge the hydrostatic pressure, instead of the non-hydrostatic pressure = rho * rd * tempv)
              pres_tend = wfac_old * p_latbc_old%pres(jc,jk,jb) + wfac_new * p_latbc_new%pres(jc,jk,jb) &
                &       - p_diag%pres(jc,jk,jb)
              temp_tend = wfac_old * p_latbc_old%temp(jc,jk,jb) + wfac_new * p_latbc_new%temp(jc,jk,jb) &
                &       - p_diag%temp(jc,jk,jb)
              
              ! Transform water vapor and temperature increments into a virtual temperature increment
              ! (see eq. (4) above)
              tempv_tend = ( 1._wp + vtmpc1 * p_prog_rcf%tracer(jc,jk,jb,iqv) ) * temp_tend & 
                &        + p_diag%temp(jc,jk,jb) * vtmpc1 * qv_tend(jc,jk,jb)
              
              ! Transform increments of virtual temperature and pressure into increments 
              ! of density and virtual potential temperature (see eqs. (6) and (8) above)
              rho_tend = ( rrd * pres_tend - p_prog%rho(jc,jk,jb) * tempv_tend ) / p_diag%tempv(jc,jk,jb)
              thv_tend = ( tempv_tend - rcpd * pres_tend / p_prog%rho(jc,jk,jb) ) / p_prog%exner(jc,jk,jb)
              
              ! Nudging step
              ! (Please note: we assume that the nudging update is well-behaved, 
              ! i.e. it will never lead to negative values for the density and the virtual potential temperature. 
              ! For reasons of computational efficiency, we have to refrain from countermeasures (e.g. ... = MAX(..., ...)).
              ! However, if only one of the two, rho or theta_v, would happen to become negative, 
              ! the subsequent computation of the Exner pressure would at least point us to this fact, 
              ! since the program would crash due to a negative argument of the LOG function.)
              p_prog%rho(jc,jk,jb)     = p_prog%rho(jc,jk,jb)     + nudge_coeff_thermdyn(jk) * rho_tend
              p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) + nudge_coeff_thermdyn(jk) * thv_tend
              
              ! Update Exner pressure
              p_prog%exner(jc,jk,jb) = EXP( rd_o_cvd * &
                & LOG( rd_o_p0ref * p_prog%rho(jc,jk,jb) * p_prog%theta_v(jc,jk,jb) ) )
              
            ENDDO  !jc
          ENDDO  !jk
        ENDDO  !jb
!$OMP END DO
        
      ELSEIF (l_thermdyn .AND. .NOT. l_hydrostatic) THEN
        
        !---------------------------------------------------------------
        !    CASE 2: Nudge density and virtual potential temperature 
        !---------------------------------------------------------------
        
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,thv_tend,rho_tend) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk
          
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          
          DO jk = istart, iend
            DO jc = i_startidx, i_endidx
              
              thv_tend = wfac_old * p_latbc_old%theta_v(jc,jk,jb) + wfac_new * p_latbc_new%theta_v(jc,jk,jb) &
                &      - p_prog%theta_v(jc,jk,jb)
              rho_tend = wfac_old * p_latbc_old%rho(jc,jk,jb) + wfac_new * p_latbc_new%rho(jc,jk,jb)         &
                &      - p_prog%rho(jc,jk,jb)
              
              p_prog%rho(jc,jk,jb)     = p_prog%rho(jc,jk,jb)     + nudge_coeff_thermdyn(jk) * rho_tend
              p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) + nudge_coeff_thermdyn(jk) * thv_tend
              p_prog%exner(jc,jk,jb)   = EXP( rd_o_cvd * &
                & LOG( rd_o_p0ref * p_prog%rho(jc,jk,jb) * p_prog%theta_v(jc,jk,jb) ) )
              
            ENDDO  !jc
          ENDDO  !jk
        ENDDO  !jb
!$OMP END DO
        
      ENDIF  !What thermodynamic variables are to be nudged?
      
      !---------------------------------------------------------------
      !                    Nudging of water vapor
      !---------------------------------------------------------------
      
      IF (l_qv) THEN
        
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,qv_update) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk
          
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          
          DO jk = istart, iend
            DO jc = i_startidx, i_endidx
              
              ! Suppress positive nudging tendencies in saturated (cloudy) regions,
              ! in order to avoid runaway effects
              qv_update = MERGE( nudge_coeff_qv(jk) * qv_tend(jc,jk,jb),        &  ! Yes 
                &                MERGE( nudge_coeff_qv(jk) * qv_tend(jc,jk,jb), &  ! {Yes                            }
                &                       0._wp,                                  &  ! {No                             } No
                &                       qv_tend(jc,jk,jb) < 0._wp ),            &  ! {Negative water vapor increment?}
                &                p_prog_rcf%tracer(jc,jk,jb,iqc) < eps_qc       )  ! Cell cloud-free?
              
              p_prog_rcf%tracer(jc,jk,jb,iqv) = p_prog_rcf%tracer(jc,jk,jb,iqv) + qv_update
              
            ENDDO  !jc
          ENDDO  !jk
        ENDDO  !jb
!$OMP END DO
        
      ENDIF  !IF (l_qv)
      
      !...............................................................
      !                Nudging of edge-based variables
      !...............................................................
      
      ! The horizontal loop runs over the entire prognostic domain + halo edges, 
      ! in order to avoid a synchronization afterwards
      rl_start   = grf_bdywidth_e + 1
      rl_end     = min_rledge
      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)
      
      !---------------------------------------------------------------
      !                   Nudging of horizontal wind
      !---------------------------------------------------------------
      
      IF (l_vn) THEN
        
!$OMP DO PRIVATE(jk,je,jb,i_startidx,i_endidx,vn_tend) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk
          
          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          
          DO jk = istart, iend
            DO je = i_startidx, i_endidx
              
              vn_tend = wfac_old * p_latbc_old%vn(je,jk,jb) + wfac_new * p_latbc_new%vn(je,jk,jb) &
                &     - p_prog%vn(je,jk,jb)
              
              p_prog%vn(je,jk,jb) = p_prog%vn(je,jk,jb) + nudge_coeff_vn(jk) * vn_tend
              
            ENDDO  !je
          ENDDO  !jk
        ENDDO  !jb
!$OMP END DO
        
      ENDIF  !IF (l_vn)
      
!$OMP END PARALLEL
      
      ! Clean-up
      IF (ALLOCATED(nudge_coeff_thermdyn)) THEN 
        DEALLOCATE(nudge_coeff_thermdyn, STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Deallocation of nudge_coeff_thermdyn failed!')
      ENDIF
      IF (ALLOCATED(nudge_coeff_vn)) THEN 
        DEALLOCATE(nudge_coeff_vn, STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Deallocation of nudge_coeff_vn failed!')
      ENDIF
      IF (ALLOCATED(nudge_coeff_qv)) THEN 
        DEALLOCATE(nudge_coeff_qv, STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Deallocation of nudge_coeff_qv failed!')
      ENDIF
      IF (ALLOCATED(qv_tend)) THEN 
        DEALLOCATE(qv_tend, STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Deallocation of qv_tend failed!')
      ENDIF
      
    ENDIF  !IF (jg == 1 .AND. nudging_config%ltype(indg_type%globn))

  END SUBROUTINE global_nudging

!---------------------------------------------------------------------------

  !>
  !! This routine tries to analyze the nudging success/impact.
  !!
  !! @par Revision History
  !! This subroutine is a conceptual copy of:
  !! - 'src/atm_dyn_iconam/mo_nh_supervise'
  !! Modifications for global nudging by Sebastian Borchert, DWD (2018-09)
  !!
  SUBROUTINE nudging_diagnostics( p_patch,        & !in
    &                             p_prog,         & !in
    &                             p_metrics,      & !in
    &                             p_diag,         & !in
    &                             p_int_state,    & !in
    &                             p_latbc_old,    & !in
    &                             p_latbc_new,    & !in
    &                             mtime_datetime, & !in
    &                             sim_time,       & !in
    &                             time_config,    & !in
    &                             wfac_old,       & !in
    &                             wfac_new,       & !in
    &                             l_message,      & !in
    &                             nudging_config  ) !inout

    ! In/out variables
    TYPE(t_patch),           INTENT(IN)    :: p_patch
    TYPE(t_nh_prog),         INTENT(IN)    :: p_prog
    TYPE(t_nh_metrics),      INTENT(IN)    :: p_metrics
    TYPE(t_nh_diag),         INTENT(IN)    :: p_diag
    TYPE(t_int_state),       INTENT(IN)    :: p_int_state
    TYPE(t_pi_atm),          INTENT(IN)    :: p_latbc_old, p_latbc_new
    TYPE(datetime), POINTER, INTENT(IN)    :: mtime_datetime        
    REAL(wp),                INTENT(IN)    :: sim_time              
    TYPE(t_time_config),     INTENT(IN)    :: time_config
    REAL(wp),                INTENT(IN)    :: wfac_old, wfac_new
    LOGICAL,                 INTENT(IN)    :: l_message
    TYPE(t_nudging_config),  INTENT(INOUT) :: nudging_config
    
    ! Local variables
    REAL(wp), ALLOCATABLE :: dz(:,:,:),             &
      &                      gpot(:,:),             &
      &                      exner(:,:),            &
      &                      tempv(:,:,:),          &
      &                      divh(:,:,:),           &
      &                      pres_msl(:,:),         &
      &                      pres_msl_driving(:,:), &
      &                      sum_blk(:,:),          &
      &                      sum_buffer_in(:),      &
      &                      sum_buffer_out(:),     &
      &                      output(:)
    REAL(wp) :: thetav1, thetav2, thetav3
    REAL(wp) :: rho1, rho2, rho3
    REAL(wp) :: exner1, exner2, exner3
    REAL(wp) :: inv_cell_number, dpres_msl, dpres_msl_driving
    REAL(wp) :: var_prod, cell_vol
    INTEGER  :: jg, jb, jc, jk
    INTEGER  :: nblks_c, nlev, nlevp1
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: rl_start, rl_end
    INTEGER  :: istat
    INTEGER  :: mpi_comm
    LOGICAL  :: l_stdio_process, l_file_opened, l_close_file
    TYPE(datetime),  POINTER :: mtime_next_call
    TYPE(timedelta), POINTER :: model_time_step
    !
    INTEGER,  PARAMETER :: IPMSL_MEAN     = 1
    INTEGER,  PARAMETER :: IPMSL_MEAN_DRV = 2
    INTEGER,  PARAMETER :: ICELL_NUMBER   = 3
    INTEGER,  PARAMETER :: IDIV_MEAN      = 4
    INTEGER,  PARAMETER :: ICELL_VOL      = 5
    INTEGER,  PARAMETER :: NITEM_1        = 5 ! (= ICELL_VOL)
    !
    INTEGER,  PARAMETER :: IVAR_PMSL      = 1
    INTEGER,  PARAMETER :: IVAR_PMSL_DRV  = 2
    INTEGER,  PARAMETER :: ICOV_PMSL      = 3
    INTEGER,  PARAMETER :: NITEM_2        = 3 ! (= ICOV_PMSL) 
    !
    INTEGER,  PARAMETER :: ISIM_TIME      = 1
    INTEGER,  PARAMETER :: ICORREL        = 2
    INTEGER,  PARAMETER :: IDIV_MEAN_2    = 3
    INTEGER,  PARAMETER :: IDPSDT         = 4
    INTEGER,  PARAMETER :: NITEM_3        = 4 ! (= IDPSDT)
    !
    REAL(wp), PARAMETER :: rd_o_cvd       = 1._wp / cvd_o_rd
    REAL(wp), PARAMETER :: rd_o_p0ref     = rd / p0ref
    REAL(wp), PARAMETER :: eps            = dbl_eps * 1000._wp
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER :: &
      & routine = modname//":nudging_diagnostics"

    !----------------------------------------

    ! - To measure the nudging success, we use the correlation 
    !   between the mean sea-level pressure from ICON on the one hand 
    !   and from the driving model on the other hand.
    !   The correlation is averaged over the entire surface 
    !   (of the mean sea level). For the averaging, we could weight 
    !   each pressure point on the surface with the cell area 
    !   (-> p_patch%cells%area(jc,jb)). However, in order to reduce 
    !   the computational costs, we assume that it is sufficient 
    !   for our purposes to regard the cell areas as being equal. 
    !   So the averaging reduces to summing the correlations and 
    !   dividing the sum by the number of cells (of the horizontal grid).
    !
    ! - To measure the nudging impact, we use the divergence
    !   of the horizontal wind, as it is typically used as a measure 
    !   for gravity wave activity, i.e. atmospheric "noise", 
    !   which might be either intensified or damped by the nudging 
    !   (both cases are undesirable, if simulations with nudging 
    !   are used in the context of gravity wave investigations).
    !   The absolute value of the divergence is averaged 
    !   over the entire model atmosphere. For the averaging, 
    !   one could weight each divergence point with the cell volume 
    !   (-> p_patch%cells%area(jc,jb) * p_metrics%ddqz_z_full(jc,jk,jb) * ...
    !   ... * p_metrics%deepatmo_t1mc(jk,idamtr%t1mc%vol)). 
    !   As in case of the correlation, we neglect the variation of the cell areas.

    ! Domain index
    jg = p_patch%id

    ! Number of grid layers
    nlev = p_patch%nlev
    
    ! Although this query should have been likely covered before calling this subroutine, 
    ! we might feel better, if we make sure that:
    IF (jg == 1 .AND. nudging_config%ltype(indg_type%globn) .AND. nlev >= 3) THEN

      ! Is this process the standard I/O process?
      l_stdio_process = my_process_is_stdio()
      ! ASCII file for output already opened?
      l_file_opened  = nudging_config%diag_kit%lfileopened
      ! Communicator of relevant processes running this subroutine
      mpi_comm = get_my_mpi_work_communicator()
      
      !---------------------------------------------------------------
      !                Prepare output to ASCII file
      !---------------------------------------------------------------
      
      ! If the current call is the first call, we try to open the file 
      ! for the output of the diagnostics
      IF (l_stdio_process .AND. (.NOT. l_file_opened)) THEN
        ! Determine output unit
        nudging_config%diag_kit%fileunit = find_next_free_unit(100,1000)
        ! Try to open the file
        OPEN( UNIT   = nudging_config%diag_kit%fileunit,       &
          &   FILE   = TRIM(nudging_config%diag_kit%filename), &
          &   ACTION = "WRITE",                                &
          &   FORM   = "FORMATTED",                            &
          &   IOSTAT = istat                                   )
        IF (istat/=SUCCESS) CALL finish(routine,'Could not open '//TRIM(nudging_config%diag_kit%filename))
        nudging_config%diag_kit%lfileopened = .TRUE.
        l_file_opened                       = .TRUE.
        ! Write header
        WRITE(nudging_config%diag_kit%fileunit, '(4(a,1x))') 'Simulation_time_(s)', &
          & 'Mean_sea_level_pressure_correlation_(1)', 'Noise_by_<|div|>_(s-1)',    &
          & 'Noise_by_<|dPS/dt|>_(Pa_s-1)'
        IF (l_message) THEN
          WRITE(message_text,'(a)') 'ASCII file '//TRIM(nudging_config%diag_kit%filename)//' successfully opened.'
          CALL message(TRIM(routine), message_text)
        ENDIF
      ENDIF  !IF (l_stdio_process ...)
      
      !---------------------------------------------------------------
      !                     Compute diagnostics
      !---------------------------------------------------------------
      
      nlevp1  = p_patch%nlevp1
      nblks_c = p_patch%nblks_c

      ALLOCATE( dz(3,nproma,nblks_c),             &
        &       gpot(nproma,nblks_c),             &
        &       exner(nproma,nblks_c),            &
        &       tempv(3,nproma,nblks_c),          &
        &       pres_msl(nproma,nblks_c),         &
        &       pres_msl_driving(nproma,nblks_c), &
        &       divh(nproma,nlev,nblks_c),        &
        &       output(NITEM_3),                  &
        &       STAT=istat                        )
      IF (istat /= SUCCESS) CALL finish (routine, 'Allocation of auxiliary fields failed!')

      ! Current model time is first element of output list
      output(ISIM_TIME) = sim_time
      ! dpsdt becomes the last entry of the output list. 
      ! Some notes on dpsdt: 
      ! - Only available, if msg_level >= 11, 
      !   which is checked for in src/namelists/mo_nudging_nml: check_nudging
      ! - Computed in: src/atm_dyn_iconam/mo_nh_supervise: compute_dpsdt, 
      !   which is called by: src/atm_phy_nwp/mo_nh_interface_nwp: nwp_nh_interface
      ! - Only stdio-process has information on dpsdt, but that is ok., 
      !   since the stdio-process writes to the ASCII file and screen
      output(IDPSDT) = nudging_config%dpsdt

      ! Diagnose divergence of horizontal wind 
      ! as measure for atmospheric "noise" (-> gravity wave activity), 
      ! potentially generated or damped by the nudging
      CALL div( vec_e     = p_prog%vn,   & !in
        &       ptr_patch = p_patch,     & !in
        &       ptr_int   = p_int_state, & !in
        &       div_vec_c = divh         ) !out

      ! Diagnose mean sea-level pressure 

      ! Prognostic domain, only
      rl_start   = grf_bdywidth_c+1
      rl_end     = min_rlcell_int
      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
        
        DO jc = i_startidx, i_endidx
          dz(1,jc,jb)    = p_metrics%ddqz_z_full(jc,nlev,jb)
          dz(2,jc,jb)    = p_metrics%ddqz_z_full(jc,nlev-1,jb)
          dz(3,jc,jb)    = 0.5_wp * p_metrics%ddqz_z_full(jc,nlev-2,jb)
          gpot(jc,jb)    = grav * p_metrics%z_ifc(jc,nlevp1,jb)
          exner(jc,jb)   = p_prog%exner(jc,nlev-2,jb)
          tempv(1,jc,jb) = p_diag%tempv(jc,nlev,jb)
          tempv(2,jc,jb) = p_diag%tempv(jc,nlev-1,jb)
          tempv(3,jc,jb) = p_diag%tempv(jc,nlev-2,jb)
        ENDDO  !jc
      ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

      CALL diag_pmsl( p_patch  = p_patch, & !in
        &             dz       = dz,      & !in
        &             gpot     = gpot,    & !in
        &             exner    = exner,   & !in
        &             tempv    = tempv,   & !in
        &             pres_msl = pres_msl ) !out

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,thetav1,thetav2,thetav3, &
!$OMP            exner1,exner2,exner3,rho1,rho2,rho3) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
        
        DO jc = i_startidx, i_endidx
          thetav1 = wfac_old * p_latbc_old%theta_v(jc,nlev,jb)   &
            &     + wfac_new * p_latbc_new%theta_v(jc,nlev,jb)
          thetav2 = wfac_old * p_latbc_old%theta_v(jc,nlev-1,jb) &
            &     + wfac_new * p_latbc_new%theta_v(jc,nlev-1,jb)
          thetav3 = wfac_old * p_latbc_old%theta_v(jc,nlev-2,jb) &
            &     + wfac_new * p_latbc_new%theta_v(jc,nlev-2,jb)
          rho1    = wfac_old * p_latbc_old%rho(jc,nlev,jb)       &
            &     + wfac_new * p_latbc_new%rho(jc,nlev,jb)
          rho2    = wfac_old * p_latbc_old%rho(jc,nlev-1,jb)     &
            &     + wfac_new * p_latbc_new%rho(jc,nlev-1,jb)
          rho3    = wfac_old * p_latbc_old%rho(jc,nlev-2,jb)     &
            &     + wfac_new * p_latbc_new%rho(jc,nlev-2,jb)
          exner1 = EXP( rd_o_cvd * LOG( rd_o_p0ref * rho1 * thetav1 ) )
          exner2 = EXP( rd_o_cvd * LOG( rd_o_p0ref * rho2 * thetav2 ) )
          exner3 = EXP( rd_o_cvd * LOG( rd_o_p0ref * rho3 * thetav3 ) )
          exner(jc,jb)   = exner3
          tempv(1,jc,jb) = thetav1 * exner1
          tempv(2,jc,jb) = thetav2 * exner2
          tempv(3,jc,jb) = thetav3 * exner3
        ENDDO  !jc
        
      ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

      CALL diag_pmsl( p_patch  = p_patch,         & !in
        &             dz       = dz,              & !in
        &             gpot     = gpot,            & !in
        &             exner    = exner,           & !in
        &             tempv    = tempv,           & !in
        &             pres_msl = pres_msl_driving ) !out

      DEALLOCATE(dz, gpot, exner, tempv, STAT=istat)
      IF (istat /= SUCCESS) CALL finish (routine, 'Deallocation of auxiliary fields failed!')

      ALLOCATE(sum_blk(NITEM_1,nblks_c), STAT=istat)
      IF (istat /= SUCCESS) CALL finish (routine, 'Allocation of sum_blk failed!')

      ! Initialization
      sum_blk(:,:) = 0._wp

      ! Block-wise summation for computation of mean
      
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,cell_vol) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
        
        DO jc = i_startidx, i_endidx
          sum_blk(IPMSL_MEAN,jb)     = sum_blk(IPMSL_MEAN,jb) + pres_msl(jc,jb)
          sum_blk(IPMSL_MEAN_DRV,jb) = sum_blk(IPMSL_MEAN_DRV,jb) + pres_msl_driving(jc,jb)
          sum_blk(ICELL_NUMBER,jb)   = sum_blk(ICELL_NUMBER,jb) + 1._wp
        ENDDO  !jc
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            ! Cell volume (divided by cell area and 
            ! including a potential deep-atmosphere modification)
            cell_vol         = p_metrics%ddqz_z_full(jc,jk,jb) * &
              &                p_metrics%deepatmo_t1mc(jk,idamtr%t1mc%vol)
            ! We sum the absolute value of the horizontal wind divergence, 
            ! in order to avoid cancelation of positive and negative values 
            sum_blk(IDIV_MEAN,jb) = sum_blk(IDIV_MEAN,jb) + ABS(divh(jc,jk,jb)) * cell_vol
            sum_blk(ICELL_VOL,jb) = sum_blk(ICELL_VOL,jb) + cell_vol
          ENDDO  !jc
        ENDDO  !jk
      ENDDO  !jb
!$OMP END DO    
!$OMP END PARALLEL

      ALLOCATE(sum_buffer_in(NITEM_1), sum_buffer_out(NITEM_1), STAT=istat)
      IF (istat /= SUCCESS) CALL finish (routine, 'Allocation of sum_buffer_in/out failed!')

      ! Summation over blocks
      sum_buffer_in = SUM(sum_blk, DIM=2)

      ! Averaging over all processes
      sum_buffer_out                 = p_sum(sum_buffer_in, comm = mpi_comm)
      inv_cell_number                = 1._wp / MAX(1._wp, sum_buffer_out(ICELL_NUMBER))
      sum_buffer_out(IPMSL_MEAN)     = inv_cell_number * sum_buffer_out(IPMSL_MEAN)
      sum_buffer_out(IPMSL_MEAN_DRV) = inv_cell_number * sum_buffer_out(IPMSL_MEAN_DRV)
      IF (sum_buffer_out(ICELL_VOL) > eps) THEN
        output(IDIV_MEAN_2) = sum_buffer_out(IDIV_MEAN) / sum_buffer_out(ICELL_VOL)
      ELSE
        output(IDIV_MEAN_2) = 0._wp
      ENDIF      

      DEALLOCATE(sum_blk, STAT=istat)
      IF (istat /= SUCCESS) CALL finish (routine, 'Deallocation of sum_blk failed!')
      ALLOCATE(sum_blk(NITEM_2,nblks_c), STAT=istat)
      IF (istat /= SUCCESS) CALL finish (routine, 'Allocation of sum_blk failed!')

      ! Computation of the correlation

      ! Initialization
      sum_blk(:,:) = 0._wp
      
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,dpres_msl,dpres_msl_driving) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
        
        DO jc = i_startidx, i_endidx
          dpres_msl                 = pres_msl(jc,jb) - sum_buffer_out(IPMSL_MEAN)
          dpres_msl_driving         = pres_msl_driving(jc,jb) - sum_buffer_out(IPMSL_MEAN_DRV)
          sum_blk(IVAR_PMSL,jb)     = sum_blk(IVAR_PMSL,jb) + dpres_msl**2
          sum_blk(IVAR_PMSL_DRV,jb) = sum_blk(IVAR_PMSL_DRV,jb) + dpres_msl_driving**2
          sum_blk(ICOV_PMSL,jb)     = sum_blk(ICOV_PMSL,jb) + dpres_msl * dpres_msl_driving
        ENDDO  !jc
      ENDDO  !jb
!$OMP END DO 
!$OMP END PARALLEL

      DEALLOCATE(sum_buffer_in, sum_buffer_out, STAT=istat)
      IF (istat /= SUCCESS) CALL finish (routine, 'Deallocation of sum_buffer_in/out failed!')
      ALLOCATE(sum_buffer_in(NITEM_2), sum_buffer_out(NITEM_2), STAT=istat)
      IF (istat /= SUCCESS) CALL finish (routine, 'Allocation of sum_buffer_in/out failed!')

      ! Summation over blocks
      sum_buffer_in = SUM(sum_blk, DIM=2)

      ! Summation over all processes
      sum_buffer_out = p_sum(sum_buffer_in, comm = mpi_comm)

      ! Computation of correlation coefficient
      var_prod = SQRT(sum_buffer_out(IVAR_PMSL) * sum_buffer_out(IVAR_PMSL_DRV))
      IF (var_prod > eps) THEN
        output(ICORREL) = sum_buffer_out(ICOV_PMSL) / var_prod
      ELSE
        output(ICORREL) = 0._wp
      ENDIF

      ! If 'io_nml: inextra_2d' is set to 3, we output the mean sea-level pressure of ICON, 
      ! the mean sea-level pressure of the driving data, and the difference between the two
      IF (inextra_2d == 3) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk
          
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          
          DO jc = i_startidx, i_endidx
            p_diag%extra_2d(jc,jb,1) = pres_msl(jc,jb)
            p_diag%extra_2d(jc,jb,2) = pres_msl_driving(jc,jb)
            p_diag%extra_2d(jc,jb,3) = pres_msl(jc,jb) - pres_msl_driving(jc,jb)
          ENDDO  !jc
        ENDDO  !jb
!$OMP END DO 
!$OMP END PARALLEL

      ENDIF  !IF (inextra_2d == 3)
        
      DEALLOCATE( pres_msl, pres_msl_driving, divh,       &
        &         sum_blk, sum_buffer_in, sum_buffer_out, &
        &         STAT=istat                              )
      IF (istat /= SUCCESS) CALL finish (routine, 'Deallocation of auxiliary fields failed!')

      !---------------------------------------------------------------
      !                  Write result to ASCII file
      !---------------------------------------------------------------

      IF (l_stdio_process .AND. l_file_opened) THEN
        WRITE(UNIT=nudging_config%diag_kit%fileunit, FMT='(f20.3,f11.8,e20.12,f12.6)') output(:)
        IF (l_message) THEN
          ! (The global mean of |dPS/dt| is printed by 'src/atm_dyn_iconam/mo_nh_supervise: compute_dpsdt', 
          ! so we do not repeat that here)
          WRITE(message_text,'(a)') 'Ncount = '//TRIM(int2string(nudging_config%diag_kit%ncount)) &
            & //', sim_time = '//TRIM(real2string(output(ISIM_TIME)))                             &
            & //' s, correl = '//TRIM(real2string(output(ICORREL)))                               &
            & //', mean |div| = '//TRIM(real2string(output(IDIV_MEAN_2)))//' s-1'
          CALL message(TRIM(routine), message_text)
        ENDIF
      ENDIF  !IF (l_stdio_process .AND. l_file_opened)

      DEALLOCATE(output, STAT=istat)
      IF (istat /= SUCCESS) CALL finish (routine, 'Deallocation of output failed!')

      !---------------------------------------------------------------
      !                      Close ASCII file
      !---------------------------------------------------------------

      ! If the next call would be equal to or after the stop date, 
      ! we should close the output file
      ! (for safety reasons, we use the interval '2 * dtime * idiagnose')
      mtime_next_call => newDatetime(mtime_datetime)
      model_time_step => time_config%tc_dt_model
      mtime_next_call = mtime_next_call + model_time_step * INT(2 * nudging_config%idiagnose, c_int32_t)
      l_close_file    = ( mtime_next_call >= time_config%tc_stopdate )
      model_time_step => NULL()
      CALL deallocateDatetime(mtime_next_call)
      IF (l_close_file .AND. l_stdio_process .AND. l_file_opened) THEN
        ! Try to close the file
        CLOSE( UNIT   = nudging_config%diag_kit%fileunit, &
          &    IOSTAT = istat                             )
        IF (istat/=SUCCESS) CALL finish(routine,'Could not close '//TRIM(nudging_config%diag_kit%filename))
        nudging_config%diag_kit%lfileopened = .FALSE.
        ! Do not risk that there might be another call, 
        ! which would open the file again
        nudging_config%diag_kit%lfileclosed = .TRUE.
        IF (l_message) THEN
          WRITE(message_text,'(a)') 'ASCII file '//TRIM(nudging_config%diag_kit%filename)//' successfully closed.'
          CALL message(TRIM(routine), message_text)
        ENDIF
      ENDIF  !IF (l_close_file ...)
      ! If the file has been closed, this subroutine should not be called any more. 
      ! But so far, only the stdio-process has this information, so we have to broadcast it.
      CALL p_bcast(t_buffer=nudging_config%diag_kit%lfileclosed, p_source=process_mpi_stdio_id, comm=mpi_comm)
      
    ENDIF ! IF (jg == 1 .AND. nudging_config%ltype(indg_type%globn))

  END SUBROUTINE nudging_diagnostics

!---------------------------------------------------------------------------

  !>
  !! This subroutine diagnoses the mean sea-level pressure.
  !!
  !! @par Revision History
  !! This subroutine copies:
  !! - 'src/atm_dyn_iconam/mo_nh_diagnose_pres_temp: diag_pres'
  !! - 'src/atm_phy_nwp/mo_nh_diagnose_pmsl: diagnose_pmsl_gme'
  !! The subroutines in 'src/atm_phy_nwp/mo_nh_diagnose_pmsl' 
  !! are not suitable for our purposes, unfortunately.
  !!
  SUBROUTINE diag_pmsl( p_patch, & !in
    &                   dz,      & !in
    &                   gpot,    & !in
    &                   exner,   & !in
    &                   tempv,   & !in
    &                   pres_msl ) !out

    ! In/out variables
    TYPE(t_patch), INTENT(IN)  :: p_patch       ! 
    REAL(wp),      INTENT(IN)  :: dz(:,:,:)     ! dz(1:3,jc,jb) from ddqz_z_full(jc,nlev-2:nlev,jb)
    REAL(wp),      INTENT(IN)  :: gpot(:,:)     ! gpot(jc,jb) from z_ifc(jc,nlevp1,jb) * grav
    REAL(wp),      INTENT(IN)  :: exner(:,:)    ! exner(jc,jb) from exner(jc,nlev-2,jb)
    REAL(wp),      INTENT(IN)  :: tempv(:,:,:)  ! tempv(1:3,jc,jb) from tempv(jc,nlev-2:nlev,jb)
    REAL(wp),      INTENT(OUT) :: pres_msl(:,:) ! mean sea-level pressure

    ! Local variables
    REAL(wp) :: pres1, pres_sfc, tstar, tmsl, prt, prtal, alph, pres_ifc1
    REAL(wp) :: dz_o_tempv, rgrav_gpot, tmsl_m_tstar
    INTEGER  :: jb, jc
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: rl_start, rl_end
    REAL(wp), PARAMETER :: cpd_o_rd  = 1._wp / rd_o_cpd
    REAL(wp), PARAMETER :: grav_o_rd = grav / rd
    REAL(wp), PARAMETER :: t_low     = 255.0_wp
    REAL(wp), PARAMETER :: t_high    = 290.5_wp
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER :: &
      & routine = modname//":diag_pmsl"

    !----------------------------------------

    ! Prognostic domain, only
    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)
    
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,dz_o_tempv,rgrav_gpot,tmsl_m_tstar, &
!$OMP            pres1,pres_sfc,tstar,tmsl,prt,prtal,alph,pres_ifc1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      
      DO jc = i_startidx, i_endidx
        dz_o_tempv = dz(1,jc,jb) / tempv(1,jc,jb)
        pres_sfc   = p0ref * EXP( cpd_o_rd * LOG( exner(jc,jb) ) &
          &        + grav_o_rd * ( dz_o_tempv + dz(2,jc,jb) / tempv(2,jc,jb) & 
          &        + dz(3,jc,jb) / tempv(3,jc,jb) ) )
        pres_ifc1  = pres_sfc * EXP( -grav_o_rd * dz_o_tempv )
        pres1      = SQRT(pres_sfc * pres_ifc1)
        tstar      = ( 1._wp - dtdz_standardatm * rd * rgrav * ( pres_sfc / pres1 - 1._wp ) ) * tempv(1,jc,jb)
        IF (tstar < t_low) tstar = 0.5_wp * ( t_low + tstar )
        rgrav_gpot = rgrav * gpot(jc,jb)
        tmsl = tstar - dtdz_standardatm * rgrav_gpot
        IF (tmsl > t_high) THEN
          IF (tstar > t_high) THEN
            tstar = 0.5_wp * ( t_high + tstar )
            tmsl  = tstar
          ELSE
            tmsl  = t_high
          END IF
        END IF
        tmsl_m_tstar = tmsl - tstar
        IF (ABS(tmsl_m_tstar) < 1.E-6_wp) THEN
          alph = 0._wp
        ELSE IF (ABS(rgrav_gpot) > 1.E-4_wp) THEN
          alph = rd * tmsl_m_tstar / gpot(jc,jb)
        ELSE 
          alph = -rd * dtdz_standardatm * rgrav
        END IF
        prt             = gpot(jc,jb) / ( rd * tstar )
        prtal           = prt * alph
        pres_msl(jc,jb) = pres_sfc * EXP( prt * ( 1.0_wp - prtal * ( 0.5_wp - one_third * prtal ) ) )
      ENDDO  !jc
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL
    
  END SUBROUTINE diag_pmsl

END MODULE mo_nudging
