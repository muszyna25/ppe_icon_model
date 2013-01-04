!>
!!
!! @brief prepares and postprocesses the fields for and from nwp physics
!!
!! Depending on the action item different sets of physics will be called:
!! this is
!! 1. condensation only so to apply saturation adjustment where needed
!! 2. 'slow physics' means up to now the whole physical package beside
!!     microphysics.
!! 3. turbulence, microphysics and condensation are considered fast physical package
!! 4. Updating the moist tracers in synchrone time intervalls to
!!    advection and saturation adjustment
!!
!!
!!
!! @par Revision History
!!  first implementation by Kristina Froehlich, DWD (2009-06-12)
!!
!! $Id: n/a$
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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

MODULE mo_nh_interface_nwp

  USE mo_datetime,           ONLY: t_datetime
  USE mo_kind,               ONLY: wp

 ! USE mo_timer,              ONLY: timer_physics, timer_start, timer_stop, &
  USE mo_timer 
  USE mo_exception,          ONLY: message, message_text !, finish
  USE mo_impl_constants,     ONLY: itconv, itccov, itrad, itgscp,         &
    &                              itsatad, itupdate, itturb, itsfc, itradheat, &
    &                              itsso, itgwd, itfastphy, icc,          &
    &                              min_rlcell_int, min_rledge_int, min_rlcell
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_intp_rbf,           ONLY: rbf_vec_interpol_cell
  USE mo_intp,               ONLY: edges2cells_scalar
  USE mo_grf_intp_data_strc, ONLY: t_gridref_state!,t_gridref_single_state, &
  USE mo_model_domain,       ONLY: t_patch
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_nonhydro_types,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydrostatic_config, ONLY: kstart_moist, l_open_ubc, lhdiff_rcf
  USE mo_nwp_lnd_types,      ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_nwp_phy_types,      ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_parallel_config,    ONLY: nproma, p_test_run, use_icon_comm, use_physics_barrier

  USE mo_run_config,         ONLY: ntracer, iqv, iqc, iqi, iqr, iqs,          &
    &                              msg_level, ltimer, timers_level, nqtendphy
  USE mo_io_config,          ONLY: lflux_avg
  USE mo_physical_constants, ONLY: rd, rd_o_cpd, vtmpc1, p0ref, cvd_o_rd, rcvd, cpd

  USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp

  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
  USE mo_lnd_nwp_config,     ONLY: ntiles_total, ntiles_water
  USE mo_cover_koe,          ONLY: cover_koe
  USE mo_satad,              ONLY: satad_v_3D
  USE mo_radiation,          ONLY: radheat, pre_radiation_nwp
  USE mo_radiation_config,   ONLY: irad_aero
  USE mo_nwp_gw_interface,   ONLY: nwp_gwdrag 
  USE mo_nwp_gscp_interface, ONLY: nwp_microphysics
  USE mo_nwp_turb_interface, ONLY: nwp_turbulence
  USE mo_nwp_turb_sfc_interface, ONLY: nwp_turbulence_sfc
  USE mo_nwp_sfc_interface,  ONLY: nwp_surface
  USE mo_nwp_conv_interface, ONLY: nwp_convection
  USE mo_nwp_rad_interface,  ONLY: nwp_radiation
  USE mo_sync,               ONLY: sync_patch_array, sync_patch_array_mult, SYNC_E, &
                                   SYNC_C, SYNC_C1, global_max, global_min, global_sum_array
  USE mo_mpi,                ONLY: my_process_is_mpi_all_parallel, work_mpi_barrier
  USE mo_nwp_diagnosis,      ONLY: nwp_diagnosis
  USE mo_icon_comm_lib,     ONLY: new_icon_comm_variable, delete_icon_comm_variable, &
     & icon_comm_var_is_ready, icon_comm_sync, icon_comm_sync_all, is_ready, until_sync
!  USE mo_communication,      ONLY: time_sync
  USE mo_art_washout_interface,  ONLY:art_washout_interface
  USE mo_art_config,          ONLY:art_config
  USE mo_linked_list,         ONLY: t_var_list

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  REAL(wp), PARAMETER :: rd_o_p0ref = rd / p0ref
  REAL(wp), PARAMETER :: cpd_o_rd = 1._wp / rd_o_cpd

  PUBLIC :: nwp_nh_interface

CONTAINS
  !
  !-----------------------------------------------------------------------
  !
  SUBROUTINE nwp_nh_interface(lcall_phy_jg, lredgrid, dt_loc,      & !input
                            & dtadv_loc, dt_phy_jg,                & !input
                            & p_sim_time, datetime,                & !input
                            & pt_patch, pt_int_state, p_metrics,   & !input
                            & pt_par_patch, pt_par_int_state,      & !input
                            & pt_par_grf_state,                    & !input
                            & ext_data,                            & !input
                            & pt_prog,                             & !inout
                            & pt_prog_now_rcf, pt_prog_rcf,        & !in/inout
                            & pt_diag ,                            & !inout
                            & prm_diag, prm_nwp_tend,lnd_diag,     &
                            & lnd_prog_now, lnd_prog_new,          & !inout
                            & wtr_prog_now, wtr_prog_new,          & !inout
                            & p_prog_list                             ) !in  

    !>
    ! !INPUT PARAMETERS:

    LOGICAL, INTENT(IN)          ::   &             !< physics package time control (switches)
         &                          lcall_phy_jg(:) !< for domain jg
    LOGICAL, INTENT(IN)          :: lredgrid        !< use reduced grid for radiation
    REAL(wp),INTENT(in)          :: dt_loc          !< time step applicable to local grid level
    REAL(wp),INTENT(in)          :: dtadv_loc       !< same for advective time step
    REAL(wp),INTENT(in)          :: dt_phy_jg(:)    !< time interval for all physics
                                                    !< packages on domain jg
    REAL(wp),INTENT(in)          :: p_sim_time

    TYPE(t_datetime),            INTENT(in):: datetime
    TYPE(t_patch),        TARGET,INTENT(in):: pt_patch     !<grid/patch info.
    TYPE(t_patch),        TARGET,INTENT(in):: pt_par_patch !<grid/patch info (parent grid)

    TYPE(t_int_state),    TARGET,INTENT(in):: pt_int_state      !< interpolation state
    TYPE(t_int_state),    TARGET,INTENT(in):: pt_par_int_state  !< " for parent grid
    TYPE(t_gridref_state),TARGET,INTENT(IN):: pt_par_grf_state  !< grid refinement state

    TYPE(t_nh_metrics)   ,       INTENT(in):: p_metrics
    TYPE(t_external_data),       INTENT(inout):: ext_data

    TYPE(t_nh_diag), TARGET, INTENT(inout) :: pt_diag     !<the diagnostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog     !<the prognostic variables
    TYPE(t_nh_prog), TARGET, INTENT(IN)    :: pt_prog_now_rcf !<old state for tke
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog_rcf !<the prognostic variables (with
                                                          !< red. calling frequency for tracers!
    TYPE(t_nwp_phy_diag),       INTENT(inout) :: prm_diag
    TYPE(t_nwp_phy_tend),TARGET,INTENT(inout) :: prm_nwp_tend
    TYPE(t_lnd_prog),           INTENT(inout) :: lnd_prog_now, lnd_prog_new
    TYPE(t_wtr_prog),           INTENT(inout) :: wtr_prog_now, wtr_prog_new
    TYPE(t_lnd_diag),           INTENT(inout) :: lnd_diag

    TYPE(t_var_list), INTENT(in) :: p_prog_list !current prognostic state list

    ! !OUTPUT PARAMETERS:            !<variables induced by the whole physics
    ! Local array bounds:

    INTEGER :: nlev, nlevp1            !< number of full levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    ! Local scalars:

    INTEGER :: jc,jk,jb,jce      !block index
    INTEGER :: jg                !domain id

    LOGICAL :: ltemp, lpres, ltemp_ifc, l_any_fastphys, l_any_slowphys

    INTEGER,  POINTER ::  iidx(:,:,:), iblk(:,:,:), ieidx(:,:,:), ieblk(:,:,:)

    REAL(wp), TARGET :: &                                              !> temporal arrays for 
      & z_ddt_u_tot (nproma,pt_patch%nlev,pt_patch%nblks_c),& 
      & z_ddt_v_tot (nproma,pt_patch%nlev,pt_patch%nblks_c),& !< hor. wind tendencies
      & z_ddt_temp  (nproma,pt_patch%nlev,pt_patch%nblks_c)   !< Temperature tendency
 
    REAL(wp) :: z_exner_sv(nproma,pt_patch%nlev,pt_patch%nblks_c)

    !< vertical interfaces

    REAL(wp) :: z_airmass (nproma,pt_patch%nlev) !< needed for radheat
    REAL(wp) :: zi0       (nproma)     !< solar incoming radiation at TOA   [W/m2]
    REAL(wp) :: zsct ! solar constant (at time of year)
    REAL(wp) :: zcosmu0 (nproma,pt_patch%nblks_c)

    REAL(wp) :: rd_o_cvd
    REAL(wp) :: r_sim_time

    REAL(wp) :: z_qsum       !< summand of virtual increment
    REAL(wp) :: z_ddt_qsum   !< summand of tendency of virtual increment
    REAL(wp) :: rcld(nproma,pt_patch%nlevp1)

    ! auxiliaries for Rayleigh friction computation
    REAL(wp) :: vabs, rfric_fac, ustart, uoffset_q, ustart_q, max_relax

    ! variables for CFL diagnostic
    REAL(wp) :: maxcfl(pt_patch%nblks_c), cflmax, avg_invedgelen(nproma), csfac
    ! Variables for dpsdt diagnostic
    REAL(wp) :: dps_blk(pt_patch%nblks_c), dpsdt_avg
    INTEGER  :: npoints_blk(pt_patch%nblks_c), npoints

    ! variables for extended debug output
    REAL(wp) :: maxtke(pt_patch%nblks_c,pt_patch%nlevp1),tkemax(pt_patch%nlevp1)
    REAL(wp), DIMENSION(pt_patch%nblks_c,pt_patch%nlev) :: maxabs_u, maxabs_v, &
      maxtemp, mintemp, maxqv, minqv, maxqc, minqc, maxtturb, maxuturb, maxvturb
    REAL(wp), DIMENSION(pt_patch%nlev) :: umax, vmax, tmax, tmin, qvmax, qvmin, qcmax, &
      qcmin, tturbmax, uturbmax, vturbmax

    ! communication ids, these do not need to be different variables,
    ! since they are not treated individualy
    INTEGER :: ddt_u_tot_comm, ddt_v_tot_comm, z_ddt_u_tot_comm, z_ddt_v_tot_comm, &
      & tracers_comm, tempv_comm, exner_old_comm

    IF (ltimer) CALL timer_start(timer_physics)

    ! local variables related to the blocking

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    !define pointers
    iidx  => pt_patch%edges%cell_idx
    iblk  => pt_patch%edges%cell_blk
    ieidx => pt_patch%cells%edge_idx
    ieblk => pt_patch%cells%edge_blk

    rd_o_cvd  = 1._wp / cvd_o_rd

    ! factor for sound speed computation
    csfac = rd*cpd*rcvd

    ! Inverse of simulation time
    r_sim_time = 1._wp/MAX(1.e-6_wp, p_sim_time)

    IF (lcall_phy_jg(itsatad) .OR. lcall_phy_jg(itgscp) .OR. &
        lcall_phy_jg(itturb)  .OR. lcall_phy_jg(itsfc)) THEN
      l_any_fastphys = .TRUE.
    ELSE
      l_any_fastphys = .FALSE.
    ENDIF

    IF (lcall_phy_jg(itrad) .OR.  lcall_phy_jg(itconv) .OR. lcall_phy_jg(itccov)  &
       .OR. lcall_phy_jg(itsso) .OR. lcall_phy_jg(itgwd)) THEN
      l_any_slowphys = .TRUE.
    ELSE
      l_any_slowphys = .FALSE.
    ENDIF

    !-------------------------------------------------------------------------
    !>  Update the tracer for every advective timestep,
    !!  all other updates are done in dynamics
    !-------------------------------------------------------------------------

    IF (lcall_phy_jg(itupdate)) THEN

      IF (msg_level >= 15) &
           & CALL message('mo_nh_interface_nwp:', 'update_tracers')

      IF (timers_level > 2) CALL timer_start(timer_update_prog_phy)
      
      CALL nh_update_prog_phy(pt_patch              ,& !in
           &                  dt_phy_jg(itfastphy)  ,& !in
           &                  prm_nwp_tend          ,& !in
           &                  prm_diag              ,& !inout phyfields 
           &                  pt_prog_rcf            )!inout tracer

      IF (timers_level > 2) CALL timer_stop(timer_update_prog_phy)
    
    ENDIF

    IF ( lcall_phy_jg(itturb) .OR. lcall_phy_jg(itconv) .OR. &
         lcall_phy_jg(itsso)  .OR. lcall_phy_jg(itgwd) ) THEN
    
      !-------------------------------------------------------------------------
      !>
      !!   Interpolation from v_n onto u,v =>  Reconstruct u and v
      !!   This is needed for turbulence, convection and SSO/GWdrag
      !!
      !-------------------------------------------------------------------------

      IF (msg_level >= 15) &
           & CALL message('mo_nh_interface_nwp:', 'reconstruct u/v')

      IF (timers_level > 3) CALL timer_start(timer_phys_u_v)
      
      SELECT CASE (pt_patch%cell_type)
      CASE (3)
        CALL rbf_vec_interpol_cell(pt_prog%vn,            & !< normal wind comp.
          &                        pt_patch,              & !< patch
          &                        pt_int_state,          & !< interpolation state
          &                        pt_diag%u, pt_diag%v )   !<  reconstr. u,v wind
      CASE (6)
        CALL edges2cells_scalar(pt_prog%vn,pt_patch, &
          &                     pt_int_state%hex_east ,pt_diag%u)
        CALL edges2cells_scalar(pt_prog%vn,pt_patch, &
          &                     pt_int_state%hex_north,pt_diag%v)
      END SELECT

      IF (timers_level > 3) CALL timer_stop(timer_phys_u_v)

    ENDIF ! diagnose u/v

    IF (l_any_fastphys) THEN

      ! Diagnose temperature if any of the fast physics schemes is called
      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf,    &
           &                              pt_diag, pt_patch,       &
           &                              opt_calc_temp=.TRUE.,    &
           &                              opt_calc_pres=.FALSE.,   &
           &                              opt_rlend=min_rlcell_int )

      IF (msg_level >= 20) THEN ! Initial debug output

        CALL message('mo_nh_interface_nwp:','Initial debug output')

        maxabs_u(:,:) = 0._wp
        maxabs_v(:,:) = 0._wp
        maxtemp(:,:)  = 0._wp
        mintemp(:,:)  = 1.e20_wp
        maxqv(:,:)    = 0._wp
        minqv(:,:)    = 1.e20_wp
        maxqc(:,:)    = 0._wp
        minqc(:,:)    = 1.e20_wp

        rl_start = grf_bdywidth_c+1
        rl_end   = min_rlcell_int

        i_startblk = pt_patch%cells%start_blk(rl_start,1)
        i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx
                maxabs_u(jb,jk) = MAX(maxabs_u(jb,jk),ABS(pt_diag%u(jc,jk,jb)))
                maxabs_v(jb,jk) = MAX(maxabs_v(jb,jk),ABS(pt_diag%v(jc,jk,jb)))
                maxtemp(jb,jk)  = MAX(maxtemp(jb,jk),pt_diag%temp(jc,jk,jb))
                mintemp(jb,jk)  = MIN(mintemp(jb,jk),pt_diag%temp(jc,jk,jb))
                maxqv(jb,jk)    = MAX(maxqv(jb,jk),pt_prog_rcf%tracer(jc,jk,jb,iqv))
                minqv(jb,jk)    = MIN(minqv(jb,jk),pt_prog_rcf%tracer(jc,jk,jb,iqv))
                maxqc(jb,jk)    = MAX(maxqc(jb,jk),pt_prog_rcf%tracer(jc,jk,jb,iqc))
                minqc(jb,jk)    = MIN(minqc(jb,jk),pt_prog_rcf%tracer(jc,jk,jb,iqc))
              ENDDO
            ENDDO

        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        DO jk = 1, nlev
          umax(jk)  = MAXVAL(maxabs_u(:,jk))
          vmax(jk)  = MAXVAL(maxabs_v(:,jk))
          tmax(jk)  = MAXVAL(maxtemp(:,jk))
          tmin(jk)  = MINVAL(mintemp(:,jk))
          qvmax(jk) = MAXVAL(maxqv(:,jk))
          qvmin(jk) = MINVAL(minqv(:,jk))
          qcmax(jk) = MAXVAL(maxqc(:,jk))
          qcmin(jk) = MINVAL(minqc(:,jk))
        ENDDO

        ! Finally take maximum/minimum over all PEs
        umax  = global_max(umax)
        vmax  = global_max(vmax)
        tmax  = global_max(tmax)
        tmin  = global_min(tmin)
        qvmax = global_max(qvmax)
        qvmin = global_min(qvmin)
        qcmax = global_max(qcmax)
        qcmin = global_min(qcmin)

        WRITE(message_text,'(a,i2)') 'max |U|, max |V|, min/max T, min/max QV,&
          & max QC per level in domain ',jg
        CALL message('', TRIM(message_text))
        DO jk = 1, nlev
          WRITE(message_text,'(a,i3,7(a,e12.5))') 'level ',jk,': u =',umax(jk),', v =',vmax(jk), &
            ', t =', tmin(jk),' ', tmax(jk),', qv =', qvmin(jk),' ', qvmax(jk), &
            ', qc =', qcmax(jk)   !,' ',qcmin(jk)
          CALL message('', TRIM(message_text))
        ENDDO

      ENDIF ! debug output for msg_level >= 20

    ENDIF ! fast physics activated

    !!-------------------------------------------------------------------------
    !> Initial saturation adjustment (a second one follows at the end of the microphysics)
    !!-------------------------------------------------------------------------

    IF (lcall_phy_jg(itsatad)) THEN

      IF (msg_level >= 15) CALL message('mo_nh_interface_nwp:', 'satad')

      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL

!$OMP WORKSHARE
        ! Store exner function for sound-wave reduction and open upper boundary condition
        ! this needs to be done for all grid points (including halo points)
        z_exner_sv(:,:,:) = pt_prog%exner(:,:,:)
!$OMP END WORKSHARE

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_qsum) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        !-------------------------------------------------------------------------
        !   call the saturation adjustment
        !-------------------------------------------------------------------------

!#ifdef __BOUNDCHECK
          IF (timers_level > 2) CALL timer_start(timer_satad_v_3D)
          CALL satad_v_3D( &
               & maxiter  = 10                             ,& !> IN
               & tol      = 1.e-3_wp                       ,& !> IN
               & te       = pt_diag%temp       (:,:,jb)    ,& !> INOUT
               & qve      = pt_prog_rcf%tracer (:,:,jb,iqv),& !> INOUT
               & qce      = pt_prog_rcf%tracer (:,:,jb,iqc),& !> INOUT
               & rhotot   = pt_prog%rho        (:,:,jb)    ,& !> IN
               & idim     = nproma                         ,& !> IN
               & kdim     = nlev                           ,& !> IN
               & ilo      = i_startidx                     ,& !> IN
               & iup      = i_endidx                       ,& !> IN
               & klo      = kstart_moist(jg)               ,& !> IN
               & kup      = nlev                            & !> IN
              !& count, errstat,                              !> OUT
               )
          IF (timers_level > 2) CALL timer_stop(timer_satad_v_3D)

        IF (timers_level > 2) CALL timer_start(timer_phys_exner)

        DO jk = kstart_moist(jg), nlev
          DO jc = i_startidx, i_endidx

            ! calculate virtual temperature from condens' output temperature
            ! taken from SUBROUTINE update_tempv_geopot in hydro_atmos/mo_ha_update_diag.f90

            z_qsum=   pt_prog_rcf%tracer (jc,jk,jb,iqc) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqi) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqr) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqs)
            
            pt_diag%tempv(jc,jk,jb) =  pt_diag%temp(jc,jk,jb)                  &
              &                   * ( 1._wp +  vtmpc1                          &
              &                   * pt_prog_rcf%tracer(jc,jk,jb,iqv)  - z_qsum )

            pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                   &
              &                     * pt_prog%rho(jc,jk,jb)*pt_diag%tempv(jc,jk,jb)))

          ENDDO
        ENDDO
        IF (timers_level > 2) CALL timer_stop(timer_phys_exner)
      ENDDO ! nblks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ELSE ! satad turned off

!$OMP PARALLEL WORKSHARE
      ! Store exner function for sound-wave reduction and open upper boundary condition
      z_exner_sv(:,:,:) = pt_prog%exner(:,:,:)
!$OMP END PARALLEL WORKSHARE

    ENDIF ! satad

    !!-------------------------------------------------------------------------
    !>  turbulent transfer and diffusion  and microphysics
    !!
    !!  Because we consider the following physical processes as fast ones
    !!  we allow here the update of prognostic variables inside the subroutines
    !!  This means that the conversion back to the ICON-prognostic variables
    !!  has to be done atferwards
    !!-------------------------------------------------------------------------

    IF (lcall_phy_jg(itgscp) .OR. lcall_phy_jg(itturb) .OR. lcall_phy_jg(itsfc)) THEN

      IF (msg_level >= 15) &
        & CALL message('mo_nh_interface_nwp:', 'diagnose pressure for fast physics')

      !-------------------------------------------------------------------------
      !> temperature and virtual temperature are already up to date:
      !! thus diagnose only pressure on main and interface levels  
      !! =>  opt_calc_pres_nh=.TRUE.
      !-------------------------------------------------------------------------
      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf, &
        & pt_diag, pt_patch,      &
        & opt_calc_temp =.FALSE., &
        & opt_calc_pres =.TRUE.,  &
        & opt_rlend=min_rlcell_int)

    ENDIF

    IF (  lcall_phy_jg(itturb)) THEN
      IF (timers_level > 1) CALL timer_start(timer_nwp_turbulence)
      IF ( atm_phy_nwp_config(jg)%inwp_turb <= 2 ) THEN
        ! Turbulence schemes not including the call to the surface scheme
        CALL nwp_turbulence (  dt_phy_jg(itfastphy),              & !>input
                              & pt_patch, p_metrics,              & !>input
                              & ext_data,                         & !>input
                              & pt_prog,                          & !>inout
                              & pt_prog_now_rcf, pt_prog_rcf,     & !>in/inout
                              & pt_diag ,                         & !>inout
                              & prm_diag,prm_nwp_tend,            & !>inout
                              & wtr_prog_now,                     & !>in
                              & lnd_prog_now,                     & !>inout 
                              & lnd_diag                          ) !>inout
      ELSE
        ! Turbulence schemes including the call to the surface scheme
        CALL nwp_turbulence_sfc (  dt_phy_jg(itfastphy),              & !>input
                                  & pt_patch, p_metrics,              & !>input
                                  & ext_data,                         & !>input
                                  & pt_prog,                          & !>inout
                                  & pt_prog_now_rcf, pt_prog_rcf,     & !>in/inout
                                  & pt_diag ,                         & !>inout
                                  & prm_diag,prm_nwp_tend,            & !>inout 
                                  & lnd_prog_now, lnd_prog_new,       & !>inout 
                                  & lnd_diag                          ) !>inout
      ENDIF
      IF (timers_level > 1) CALL timer_stop(timer_nwp_turbulence)
    ENDIF !lcall(itturb)

    !-------------------------------------------------------------------------
    !  prognostic microphysic and precipitation scheme
    !-------------------------------------------------------------------------
    IF ( lcall_phy_jg(itgscp)) THEN

      IF (msg_level >= 15) &
        & CALL message('mo_nh_interface_nwp:', 'microphysics')

      !> temperature and tracers have been updated by turbulence;
      !! an update of the pressure field is not needed because pressure
      !! is not needed at high accuracy in the microphysics scheme

      IF (timers_level > 1) CALL timer_start(timer_nwp_microphysics)

      CALL nwp_microphysics ( dt_phy_jg(itfastphy),             & !>input
                            & pt_patch, p_metrics,              & !>input
                            & pt_prog,                          & !>inout
                            & pt_prog_rcf,                      & !>inout
                            & pt_diag ,                         & !>inout
                            & prm_diag                          ) !>inout

      IF (timers_level > 1) CALL timer_stop(timer_nwp_microphysics)

    ENDIF

      IF (art_config(jg)%lart) THEN

       CALL art_washout_interface(dt_phy_jg(itfastphy),          & !>in
                  &          pt_patch,                           & !>in
                  &          p_prog_list,                        & !>in
                  &          prm_diag,                           & !>in
                  &          pt_prog%rho,                        & !>in               
                  &          pt_prog_rcf%tracer)                   !>inout             

      ENDIF !lart    

    IF (lcall_phy_jg(itsfc)) THEN

      !> temperature and tracers have been updated by microphysics;
      !! as pressure is needed only for an approximate adiabatic extrapolation
      !! of the temperature at the lowest model level towards ground level,
      !! a recalculation is not required

      CALL nwp_surface    (  dt_phy_jg(itfastphy),              & !>input
                            & pt_patch,                         & !>input
                            & ext_data,                         & !>input
                            & pt_prog_rcf,                      & !>in/inout rcf=reduced calling freq.
                            & pt_diag ,                         & !>inout
                            & prm_diag,                         & !>inout 
                            & lnd_prog_now, lnd_prog_new,       & !>inout
                            & wtr_prog_now, wtr_prog_new,       & !>inout
                            & lnd_diag                          ) !>input

    ENDIF


    IF (lcall_phy_jg(itsatad) .OR. lcall_phy_jg(itgscp) .OR. lcall_phy_jg(itturb)) THEN
      
      IF (timers_level > 1) CALL timer_start(timer_fast_phys)

      ! Remark: in the (unusual) case that satad is used without any other physics,
      ! recalculation of the thermodynamic variables is duplicated here. However,
      ! this is the easiest way to combine minimization of halo communications
      ! with a failsafe flow control

      IF (msg_level >= 15) &
        & CALL message('mo_nh_interface_nwp:', 'recalculate thermodynamic variables')


      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx, i_endidx, z_qsum ) ICON_OMP_DEFAULT_SCHEDULE

      DO jb = i_startblk, i_endblk
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end )

        !-------------------------------------------------------------------------
        !>
        !! re-calculate scalar prognostic variables out of physics variables!
        !!
        !-------------------------------------------------------------------------

        DO jk = 1, nlev
          DO jc =  i_startidx, i_endidx

            z_qsum=   pt_prog_rcf%tracer (jc,jk,jb,iqc) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqi) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqr) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqs)
            pt_diag%tempv(jc,jk,jb) =  pt_diag%temp(jc,jk,jb)          &
&                                  * ( 1._wp +  vtmpc1                 &
&                                  *  pt_prog_rcf%tracer(jc,jk,jb,iqv) &
&                                   - z_qsum )

            pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                   &
              &                     * pt_prog%rho(jc,jk,jb)*pt_diag%tempv(jc,jk,jb)))

            pt_diag%exner_old(jc,jk,jb) = pt_diag%exner_old(jc,jk,jb) + &
              pt_prog%exner(jc,jk,jb) - z_exner_sv(jc,jk,jb)

            pt_prog%theta_v(jc,jk,  jb) = pt_diag%tempv(jc,jk,jb) &
&                                       / pt_prog%exner(jc,jk,jb)

            pt_prog%rhotheta_v(jc,jk,jb) = pt_prog%theta_v(jc,jk,jb) &
&                                        * pt_prog%rho(jc,jk,jb)

          ENDDO
        ENDDO

        ! This loop needs to be split here in order to ensure that the same
        ! compiler optimization is applied as in the next loop below, where
        ! the same calculations are made for the halo points.
        DO jk = kstart_moist(jg), nlev
          DO jc =  i_startidx, i_endidx

            ! finally compute dynamical temperature tendency
            pt_diag%ddt_temp_dyn(jc,jk,jb) = cpd_o_rd*pt_diag%temp(jc,jk,jb)  &
              * pt_diag%exner_dyn_incr(jc,jk,jb)/dt_phy_jg(itfastphy)         &
              / pt_prog%exner(jc,jk,jb) 

            ! reset dynamical exner increment to zero
            ! (it is accumulated over one advective time step in solve_nh)
            pt_diag%exner_dyn_incr(jc,jk,jb) = 0._wp

          ENDDO
        ENDDO

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (timers_level > 1) CALL timer_stop(timer_fast_phys)
   
   ENDIF ! end of fast physics part

    !!-------------------------------------------------------------------------
    !!  slow physics part
    !!-------------------------------------------------------------------------


    IF (l_any_slowphys) THEN

      IF (msg_level >= 15) &
         CALL message('mo_nh_interface', 'diagnose pres/temp for slow physics')

      ! If slow physics is called without fast physics (which should happen
      ! at the initial time step only), temperature needs to be calculated
      ! Otherwise, temperature is up to date
      IF ( .NOT. (lcall_phy_jg(itgscp) .OR. lcall_phy_jg(itturb))) THEN
        ltemp = .TRUE.
      ELSE
        ltemp = .FALSE.
      ENDIF

      ! Pressure should always be rediagnosed for slow physics in order to have
      ! the correct air masses for computing heating rates
      lpres = .TRUE.

      ! Temperature at interface levels is needed if irad_aero = 5 or 6
      ! or if Ritter-Geleyn radiation is called
      IF ( lcall_phy_jg(itrad) .AND. ( irad_aero == 5 .OR. irad_aero == 6 &
           .OR. atm_phy_nwp_config(jg)%inwp_radiation == 2 ) )         THEN 
        ltemp_ifc = .TRUE.
      ELSE
        ltemp_ifc = .FALSE.
      ENDIF

      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf,   &
        &                      pt_diag, pt_patch,                 &
        &                      opt_calc_temp     = ltemp,         &
        &                      opt_calc_pres     = lpres,         &
        &                      lnd_prog          = lnd_prog_new,  &
        &                      opt_calc_temp_ifc = ltemp_ifc,     &
        &                      opt_rlend         = min_rlcell_int )

    ENDIF

    IF ( lcall_phy_jg(itconv)  ) THEN

      IF (msg_level >= 15) &
&           CALL message('mo_nh_interface', 'convection')

      IF (timers_level > 2) CALL timer_start(timer_nwp_convection)
      CALL nwp_convection (  dt_phy_jg(itconv),                 & !>input
                            & pt_patch, p_metrics,              & !>input
                            & ext_data,                         & !>input
                            & pt_prog,                          & !>input
                            & pt_prog_rcf,                      & !>input
                            & pt_diag ,                         & !>inout
                            & prm_diag,prm_nwp_tend             ) !>inout 
      IF (timers_level > 2) CALL timer_stop(timer_nwp_convection)

    ENDIF! convection

    IF ( lcall_phy_jg(itccov) ) THEN

      ! When using a reduced grid for radiation, part of the boundary points need
      ! to be included in order to compute spatial gradients for back-interpolation
      IF (lredgrid) THEN
        rl_start = grf_bdywidth_c-1
      ELSE
        rl_start = grf_bdywidth_c+1
      ENDIF
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

      IF (msg_level >= 15) &
        &  CALL message('mo_nh_interface', 'cloud cover')


      !-------------------------------------------------------------------------
      !> Cloud water distribution: cloud cover, cloud water, cloud ice
      !  inwp_cldcover =    
      !  (0) no clouds
      !  (1) diagnostic cloud cover
      !  (2) prognostic total water variance (not yet started)
      !  (3) clouds as in COSMO
      !  (4) clouds as in turbulence
      !  (5) grid-scale cloud cover [1 or 0]
      !-------------------------------------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,rcld) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)

          rcld = 0.0_wp ! standard deviation of saturation deficit=0 for now, needs to be specified form turbulence

          IF (timers_level > 2) CALL timer_start(timer_cover_koe)
          CALL cover_koe &
&             (kidia  = i_startidx ,   kfdia  = i_endidx  ,       & !! in:  horizonal begin, end indices
&              klon = nproma,  kstart = kstart_moist(jg)  ,       & !! in:  horiz. and vert. vector length
&              klev   = nlev,                                     &
&              icldscheme = atm_phy_nwp_config(jg)%inwp_cldcover ,& !! in:  cloud cover option
&              tt     = pt_diag%temp         (:,:,jb)     ,       & !! in:  temperature at full levels
&              pp     = pt_diag%pres         (:,:,jb)     ,       & !! in:  pressure at full levels
&              ps     = pt_diag%pres_sfc     (:,jb)       ,       & !! in:  surface pressure at full levels
&              t_g    = lnd_prog_new%t_g     (:,jb)       ,       & !! in:  surface temperature
&              pgeo   = p_metrics%geopot_agl (:,:,jb)     ,       & !! in:  geopotential height
&              rho    = pt_prog%rho          (:,:,jb  )   ,       & !! in:  density
&              rcld   = rcld                              ,       & !! in:  standard deviation of saturation deficit
&              ldland = ext_data%atm%llsm_atm_c (:,jb)    ,       & !! in:  land/sea mask
&              ldcum  = prm_diag%locum       (:,jb)       ,       & !! in:  convection on/off
&              kcbot  = prm_diag%mbas_con    (:,jb)       ,       & !! in:  convective cloud base
&              kctop  = prm_diag%mtop_con    (:,jb)       ,       & !! in:  convective cloud top
&              pmfude_rate = prm_diag%con_udd(:,:,jb,3)   ,       & !! in:  convective updraft detrainment rate
&              plu         = prm_diag%con_udd(:,:,jb,7)   ,       & !! in:  updraft condensate 
&              qv     = pt_prog_rcf%tracer   (:,:,jb,iqv) ,       & !! in:  spec. humidity
&              qc     = pt_prog_rcf%tracer   (:,:,jb,iqc) ,       & !! in:  cloud water
&              qi     = pt_prog_rcf%tracer   (:,:,jb,iqi) ,       & !! in:  cloud ice
&              cc_tot = prm_diag%tot_cld     (:,:,jb,icc) ,       & !! out: cloud diagnostics
&              qv_tot = prm_diag%tot_cld     (:,:,jb,iqv) ,       & !! out: qv       -"-
&              qc_tot = prm_diag%tot_cld     (:,:,jb,iqc) ,       & !! out: clw      -"-
&              qi_tot = prm_diag%tot_cld     (:,:,jb,iqi) )         !! out: ci       -"-
          IF (timers_level > 2) CALL timer_stop(timer_cover_koe)

      ENDDO
  
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ENDIF! cloud cover

    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------

    IF ( lcall_phy_jg(itrad) ) THEN

      IF (ltimer) CALL timer_start(timer_nwp_radiation)
      CALL nwp_radiation (lredgrid,              & ! in
           &              p_sim_time,            & ! in
           &              datetime,              & ! in
           &              pt_patch,pt_par_patch, & ! in
           &              pt_par_int_state,      & ! in
           &              pt_par_grf_state,      & ! in
           &              ext_data,              & ! in
           &              lnd_diag,              & ! in
           &              pt_prog,               & ! inout
           &              pt_diag,prm_diag,      & ! inout
           &              lnd_prog_new,          & ! in
           &              wtr_prog_new           ) ! in
      IF (ltimer) CALL timer_stop(timer_nwp_radiation)
     
    ENDIF


    IF ( lcall_phy_jg(itradheat) ) THEN

      IF (msg_level >= 15) &
&           CALL message('mo_nh_interface', 'radiative heating')


      IF (timers_level > 1) CALL timer_start(timer_pre_radiation_nwp)

      CALL pre_radiation_nwp (                      &
        & kbdim      = nproma,                      &
        & p_inc_rad  = dt_phy_jg(itfastphy),        &
        & p_sim_time = p_sim_time,                  &
        & pt_patch   = pt_patch,                    &
        & zsmu0      = zcosmu0,                     &
        & zsct       = zsct )
      IF (timers_level > 1) CALL timer_stop(timer_pre_radiation_nwp)      

      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

      IF (timers_level > 2) CALL timer_start(timer_radheat)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,zi0,z_airmass) ICON_OMP_DEFAULT_SCHEDULE
!
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)

        zcosmu0 (i_startidx:i_endidx,jb) = MAX(zcosmu0(i_startidx:i_endidx,jb),0.0_wp)

        !calculate solar incoming flux at TOA
        zi0 (i_startidx:i_endidx) = zcosmu0(i_startidx:i_endidx,jb) * zsct !zsct by pre_radiation
        prm_diag%flxdwswtoa(i_startidx:i_endidx,jb) = zi0 (i_startidx:i_endidx)
        z_airmass(i_startidx:i_endidx,:) = p_metrics%ddqz_z_full(i_startidx:i_endidx,:,jb) * &
                                           pt_prog%rho(i_startidx:i_endidx,:,jb)

        prm_diag%swflxsfc (:,jb)=0._wp
        prm_diag%lwflxsfc (:,jb)=0._wp
        prm_diag%swflxtoa (:,jb)=0._wp

        IF (atm_phy_nwp_config(jg)%inwp_surface >= 1) THEN

          prm_diag%swflxsfc_t (:,jb,:)=0._wp
          prm_diag%lwflxsfc_t (:,jb,:)=0._wp

          CALL radheat (                   &
          !
          ! input
          ! -----
          !
          & jcs=i_startidx                         ,&! in     start index of inner do loop
          & jce=i_endidx                           ,&! in     end index of inner do loop
          & kbdim=nproma                           ,&! in     loop length and dimension size
          & klev=nlev                              ,&! in     vertical dimension size
          & klevp1=nlevp1                          ,&! in     vertical dimension size
          & ntiles=ntiles_total                    ,&! in     number of tiles of sfc flux fields
          & ntiles_wtr=ntiles_water                ,&! in     number of extra tiles for ocean and lakes
          & pmair=z_airmass                        ,&! in     layer air mass             [kg/m2]
          & pqv=pt_prog_rcf%tracer(:,:,jb,iqv)     ,&! in     specific moisture           [kg/kg]
          & pi0=zi0                                ,&! in     solar incoming flux at TOA  [W/m2]
          & pemiss=ext_data%atm%emis_rad(:,jb)     ,&! in     lw sfc emissivity
          & pqc=prm_diag%tot_cld    (:,:,jb,iqc)   ,&! in     specific cloud water        [kg/kg]
          & pqi=prm_diag%tot_cld    (:,:,jb,iqi)   ,&! in     specific cloud ice          [kg/kg]
          & ppres_ifc=pt_diag%pres_ifc(:,:,jb)     ,&! in     pressure at layer boundaries [Pa]
          & albedo=prm_diag%albvisdif(:,jb),        &! in     grid-box average albedo
          & albedo_t=prm_diag%albvisdif_t(:,jb,:),  &! in     tile-specific albedo
          & lp_count=ext_data%atm%lp_count(jb),     &! in     number of land points
          & gp_count_t=ext_data%atm%gp_count_t(jb,:),&! in   number of land points per tile
          & spi_count =ext_data%atm%spi_count(jb)  ,&! in     number of seaice points
          & idx_lst_lp=ext_data%atm%idx_lst_lp(:,jb), &! in   index list of land points
          & idx_lst_t=ext_data%atm%idx_lst_t(:,jb,:), &! in   index list of land points per tile
          & idx_lst_spi=ext_data%atm%idx_lst_spi(:,jb),&! in  index list of seaice points
          & cosmu0=zcosmu0(:,jb),                   &! in     cosine of solar zenith angle
          & opt_nh_corr=.TRUE.                     ,&! in     switch for NH mode
          & ptsfc=lnd_prog_new%t_g(:,jb)           ,&! in     surface temperature         [K]
          & ptsfc_t=lnd_prog_new%t_g_t(:,jb,:)     ,&! in     tile-specific surface temperature         [K]
          & ptsfctrad=prm_diag%tsfctrad(:,jb)      ,&! in     sfc temp. used for pflxlw   [K]
          & ptrmsw=prm_diag%trsolall (:,:,jb)      ,&! in     shortwave net tranmissivity []
          & pflxlw=prm_diag%lwflxall (:,:,jb)      ,&! in     longwave net flux           [W/m2]
          & opt_use_cv = .TRUE.                    ,&! in     use cv for computing heating rate
          !
          ! output
          ! ------
          !
          & pdtdtradsw=prm_nwp_tend%ddt_temp_radsw(:,:,jb),&! out    rad. heating by SW        [K/s]
          & pdtdtradlw=prm_nwp_tend%ddt_temp_radlw(:,:,jb),&! out    rad. heating by lw        [K/s]
          & pflxsfcsw =prm_diag%swflxsfc (:,jb)   ,&        ! out shortwave surface net flux [W/m2]
          & pflxsfclw =prm_diag%lwflxsfc (:,jb)   ,&        ! out longwave surface net flux  [W/m2]
          & pflxsfcsw_t=prm_diag%swflxsfc_t (:,jb,:)   ,&   ! out tile-specific shortwave surface net flux [W/m2]
          & pflxsfclw_t=prm_diag%lwflxsfc_t (:,jb,:)   ,&   ! out tile-specific longwave surface net flux  [W/m2]
          & pflxtoasw =prm_diag%swflxtoa (:,jb) )           ! out shortwave toa net flux     [W/m2]

        ELSE
          CALL radheat (                   &
          !
          ! input
          ! -----
          !
          & jcs=i_startidx                         ,&! in     start index of inner do loop
          & jce=i_endidx                           ,&! in     end index of inner do loop
          & kbdim=nproma                           ,&! in     loop length and dimension size
          & klev=nlev                              ,&! in     vertical dimension size
          & klevp1=nlevp1                          ,&! in     vertical dimension size
          & ntiles=1                               ,&! in     number of tiles of sfc flux fields
          & ntiles_wtr=0                           ,&! in     number of extra tiles for ocean and lakes
          & pmair=z_airmass                        ,&! in     layer air mass             [kg/m2]
          & pqv=pt_prog_rcf%tracer(:,:,jb,iqv)     ,&! in     specific moisture           [kg/kg]
          & pi0=zi0                                ,&! in     solar incoming flux at TOA  [W/m2]
          & pemiss=ext_data%atm%emis_rad(:,jb)     ,&! in     lw sfc emissivity
          & pqc=prm_diag%tot_cld    (:,:,jb,iqc)   ,&! in     specific cloud water        [kg/kg]
          & pqi=prm_diag%tot_cld    (:,:,jb,iqi)   ,&! in     specific cloud ice          [kg/kg]
          & ppres_ifc=pt_diag%pres_ifc(:,:,jb)     ,&! in     pressure at layer boundaries [Pa]
          & cosmu0=zcosmu0(:,jb),                   &! in     cosine of solar zenith angle
          & opt_nh_corr=.TRUE.                     ,&! in     switch for NH mode
          & ptsfc=lnd_prog_new%t_g(:,jb)           ,&! in     surface temperature         [K]
          & ptsfctrad=prm_diag%tsfctrad(:,jb)      ,&! in     sfc temp. used for pflxlw   [K]
          & ptrmsw=prm_diag%trsolall (:,:,jb)      ,&! in     shortwave net tranmissivity []
          & pflxlw=prm_diag%lwflxall (:,:,jb)      ,&! in     longwave net flux           [W/m2]
          & opt_use_cv = .TRUE.                    ,&! in     use cv for computing heating rate
          !
          ! output
          ! ------
          !
          & pdtdtradsw=prm_nwp_tend%ddt_temp_radsw(:,:,jb),&! out    rad. heating by SW        [K/s]
          & pdtdtradlw=prm_nwp_tend%ddt_temp_radlw(:,:,jb),&! out    rad. heating by lw        [K/s]
          & pflxsfcsw =prm_diag%swflxsfc (:,jb)   ,&        ! out shortwave surface net flux [W/m2]
          & pflxsfclw =prm_diag%lwflxsfc (:,jb)   ,&        ! out longwave surface net flux  [W/m2]
          & pflxtoasw =prm_diag%swflxtoa (:,jb) )           ! out shortwave toa net flux     [W/m2]

        ENDIF

        IF ( p_sim_time > 1.e-6_wp .AND. lflux_avg) THEN

         !sum up for averaged fluxes
          !T.R.: this is not correct for output after 1st timestep,
          !e.g. dt_phy_jg(itradheat) may then be greater than p_sim_time
          !leading to wrong averaging.
         DO jc =  i_startidx, i_endidx

          prm_diag%swflxsfc_a(jc,jb) = ( prm_diag%swflxsfc_a(jc,jb)                     &
                                 &  * (p_sim_time - dt_phy_jg(itfastphy))               &
                                 &  + dt_phy_jg(itfastphy) * prm_diag%swflxsfc(jc,jb))  &
                                 &  * r_sim_time
          prm_diag%lwflxsfc_a(jc,jb) = ( prm_diag%lwflxsfc_a(jc,jb)                     &
                                 &  * (p_sim_time - dt_phy_jg(itfastphy))               &
                                 &  + dt_phy_jg(itfastphy) * prm_diag%lwflxsfc(jc,jb))  &
                                 &  * r_sim_time
          prm_diag%swflxtoa_a(jc,jb) = ( prm_diag%swflxtoa_a(jc,jb)                     &
                                 &  * (p_sim_time - dt_phy_jg(itfastphy))               &
                                 &  + dt_phy_jg(itfastphy) * prm_diag%swflxtoa(jc,jb))  &
                                 &  * r_sim_time
          prm_diag%lwflxtoa_a(jc,jb) = ( prm_diag%lwflxtoa_a(jc,jb)                     &
                                 &  * (p_sim_time - dt_phy_jg(itfastphy))               &
                                 & + dt_phy_jg(itfastphy) * prm_diag%lwflxall(jc,1,jb)) &
                                 &  * r_sim_time
         ENDDO

        ELSEIF ( .NOT. lflux_avg ) THEN

         DO jc =  i_startidx, i_endidx

          prm_diag%swflxsfc_a(jc,jb) = prm_diag%swflxsfc_a(jc,jb)                    &
                                & + dt_phy_jg(itfastphy) * prm_diag%swflxsfc(jc,jb)
          prm_diag%lwflxsfc_a(jc,jb) = prm_diag%lwflxsfc_a(jc,jb)                    &
                                & + dt_phy_jg(itfastphy) * prm_diag%lwflxsfc(jc,jb)
          prm_diag%swflxtoa_a(jc,jb) = prm_diag%swflxtoa_a(jc,jb)                    &
                                & + dt_phy_jg(itfastphy) * prm_diag%swflxtoa(jc,jb)
          prm_diag%lwflxtoa_a(jc,jb) = prm_diag%lwflxtoa_a(jc,jb)                    &
                                & + dt_phy_jg(itfastphy) * prm_diag%lwflxall(jc,1,jb)
         END DO


        END IF

      ENDDO ! blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (timers_level > 2) CALL timer_stop(timer_radheat)
      
    ENDIF  ! inwp_radiation

    !-------------------------------------------------------------------------
    !  Gravity waves drag: orographic and non-orographic
    !-------------------------------------------------------------------------

    IF (lcall_phy_jg(itsso) .OR. lcall_phy_jg(itgwd)) THEN

      IF (msg_level >= 15) &
        &  CALL message('mo_nh_interface', 'gravity waves')

      IF (timers_level > 3) CALL timer_start(timer_sso)

      CALL nwp_gwdrag ( dt_phy_jg(itsso),          & !>input
        &               lcall_phy_jg(itsso),       & !>input
        &               dt_phy_jg(itgwd),          & !>input
        &               lcall_phy_jg(itgwd),       & !>input
        &               pt_patch,p_metrics,        & !>input
        &               ext_data,                  & !>input
        &               pt_diag ,                  & !>inout
        &               prm_diag,prm_nwp_tend      ) !>inout

      IF (timers_level > 3) CALL timer_stop(timer_sso)
    ENDIF ! inwp_sso
    !-------------------------------------------------------------------------
     
    IF (timers_level > 2) CALL timer_start(timer_phys_acc)
    !-------------------------------------------------------------------------
    !>  accumulate tendencies of slow_physics
    !-------------------------------------------------------------------------
    IF (l_any_slowphys .OR. lcall_phy_jg(itradheat)) THEN

      IF (p_test_run) THEN
        z_ddt_u_tot = 0._wp
        z_ddt_v_tot = 0._wp
      ENDIF

      IF (timers_level > 3) CALL timer_start(timer_phys_acc_1)

      ! Coefficients for extra Rayleigh friction
      ustart    = atm_phy_nwp_config(jg)%ustart_raylfric
      uoffset_q = ustart + 40._wp
      ustart_q  = ustart + 50._wp
      max_relax = -1._wp/atm_phy_nwp_config(jg)%efdt_min_raylfric

      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)
      
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx, i_endidx , z_qsum, z_ddt_qsum, vabs, &
!$OMP  rfric_fac) ICON_OMP_DEFAULT_SCHEDULE
!
      DO jb = i_startblk, i_endblk
!
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)


        ! artificial Rayleigh friction: active if GWD or SSO scheme is active
        IF (atm_phy_nwp_config(jg)%inwp_sso > 0 .OR. atm_phy_nwp_config(jg)%inwp_gwd > 0) THEN
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              vabs = SQRT(pt_diag%u(jc,jk,jb)**2 + pt_diag%v(jc,jk,jb)**2)
              rfric_fac = MAX(0._wp, 8.e-4_wp*(vabs-ustart))
              IF (vabs > ustart_q) THEN
                rfric_fac = MIN(1._wp,4.e-4_wp*(vabs-uoffset_q)**2)
              ENDIF
              prm_nwp_tend%ddt_u_raylfric(jc,jk,jb) = max_relax*rfric_fac*pt_diag%u(jc,jk,jb)
              prm_nwp_tend%ddt_v_raylfric(jc,jk,jb) = max_relax*rfric_fac*pt_diag%v(jc,jk,jb)
            ENDDO
          ENDDO
        ENDIF

        ! heating related to momentum deposition by SSO, GWD and Rayleigh friction
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            prm_nwp_tend%ddt_temp_drag(jc,jk,jb) = -rcvd*(pt_diag%u(jc,jk,jb)*             &
                                                   (prm_nwp_tend%ddt_u_sso(jc,jk,jb)+      &
                                                    prm_nwp_tend%ddt_u_gwd(jc,jk,jb)+      &
                                                    prm_nwp_tend%ddt_u_raylfric(jc,jk,jb)) &
                                                   +      pt_diag%v(jc,jk,jb)*             & 
                                                   (prm_nwp_tend%ddt_v_sso(jc,jk,jb)+      &
                                                    prm_nwp_tend%ddt_v_gwd(jc,jk,jb)+      &
                                                    prm_nwp_tend%ddt_v_raylfric(jc,jk,jb)) )
          ENDDO
        ENDDO

        z_ddt_temp(i_startidx:i_endidx,:,jb) =                                               &
   &                                       prm_nwp_tend%ddt_temp_radsw(i_startidx:i_endidx,:,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_radlw(i_startidx:i_endidx,:,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_drag (i_startidx:i_endidx,:,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_pconv(i_startidx:i_endidx,:,jb)


        ! Convert temperature tendency into Exner function tendency
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx

!            z_qsum                   = SUM(pt_prog_rcf%tracer    (jc,jk,jb,iqc:iqs))
            z_qsum=   pt_prog_rcf%tracer (jc,jk,jb,iqc) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqi) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqr) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqs)
            
!            z_ddt_qsum = SUM(prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqc:iqs))

            z_ddt_qsum =   prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqc) &
              &          + prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi)

            pt_diag%ddt_exner_phy(jc,jk,jb) = rd_o_cpd / pt_prog%theta_v(jc,jk,jb)           &
              &                             * (z_ddt_temp(jc,jk,jb)                          &
              &                             *(1._wp + vtmpc1*pt_prog_rcf%tracer(jc,jk,jb,iqv)&
              &                                - z_qsum) + pt_diag%temp(jc,jk,jb)            &
              &           * (vtmpc1 * prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqv)-z_ddt_qsum ))

          ENDDO
        ENDDO

        ! Accumulate wind tendencies of slow physics
        ! Strictly spoken, this would not be necessary if only radiation was called
        ! in the current time step, but the radiation time step should be a multiple 
        ! of the convection time step anyway in order to obtain up-to-date cloud cover fields
        IF (l_any_slowphys) THEN
          z_ddt_u_tot(i_startidx:i_endidx,:,jb) =                   &
   &          prm_nwp_tend%ddt_u_gwd     (i_startidx:i_endidx,:,jb) &
   &        + prm_nwp_tend%ddt_u_raylfric(i_startidx:i_endidx,:,jb) &
   &        + prm_nwp_tend%ddt_u_sso     (i_startidx:i_endidx,:,jb) &
   &        + prm_nwp_tend%ddt_u_pconv  ( i_startidx:i_endidx,:,jb)

          z_ddt_v_tot(i_startidx:i_endidx,:,jb) =                   &
   &          prm_nwp_tend%ddt_v_gwd     (i_startidx:i_endidx,:,jb) &
   &        + prm_nwp_tend%ddt_v_raylfric(i_startidx:i_endidx,:,jb) &
   &        + prm_nwp_tend%ddt_v_sso     (i_startidx:i_endidx,:,jb) &
   &        + prm_nwp_tend%ddt_v_pconv  ( i_startidx:i_endidx,:,jb)
          ENDIF
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      IF (timers_level > 3) CALL timer_stop(timer_phys_acc_1)
    ENDIF ! slow physics tendency accumulation

    !--------------------------------------------------------
    ! Final section: Synchronization of updated prognostic variables,
    !                interpolation of u/v tendendies to edge points,
    !                and diagnostic computations
    !--------------------------------------------------------

    ! Synchronize tracers if any of the updating (fast-physics) processes was active.
    ! In addition, tempv needs to be synchronized, and in case of lhdiff_rcf, also exner_old
    IF (l_any_fastphys) THEN

      IF (timers_level > 3) CALL timer_start(timer_phys_sync_tracers)

      IF (use_icon_comm) THEN ! use communication library

        tracers_comm = new_icon_comm_variable(pt_prog_rcf%tracer, pt_patch%sync_cells_not_in_domain,  &
          & status=is_ready, scope=until_sync, name="pt_prog_rcf%tracer")
        tempv_comm = new_icon_comm_variable(pt_diag%tempv, pt_patch%sync_cells_not_in_domain, &
          & status=is_ready, scope=until_sync, name="pt_diag%tempv")

        IF (lhdiff_rcf) THEN
          exner_old_comm = new_icon_comm_variable(pt_diag%exner_old, &
            & pt_patch%sync_cells_not_in_domain, &
            & status=is_ready, scope=until_sync, name="pt_diag%exner_old")
        ENDIF

      ELSE
        IF (lhdiff_rcf) THEN
          CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer+2, pt_diag%tempv, &
                                     pt_diag%exner_old, f4din=pt_prog_rcf%tracer)
        ELSE
          CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer+1, pt_diag%tempv, f4din=pt_prog_rcf%tracer)
        ENDIF

      ENDIF

      IF (timers_level > 3) THEN
        CALL timer_stop(timer_phys_sync_tracers)
      ENDIF
    ENDIF
          
    !------------------------------------------------------------ 
    ! sync here the slowphys for aggregation
    !-------------------------------------------------------------------
    IF (use_physics_barrier) THEN
      CALL timer_start(timer_barrier)
      CALL work_mpi_barrier()
      CALL timer_stop(timer_barrier)
    ENDIF
    !-------------------------------------------------------------------
    IF (timers_level > 3) CALL timer_start(timer_phys_sync_ddt_u)
    IF (use_icon_comm) THEN
    
      IF (lcall_phy_jg(itturb) ) THEN
        ddt_u_tot_comm = new_icon_comm_variable(prm_nwp_tend%ddt_u_turb, &
          & pt_patch%sync_cells_one_edge_in_domain, status=is_ready, scope=until_sync, &
          & name="prm_nwp_tend%ddt_u_turb")
        ddt_v_tot_comm = new_icon_comm_variable(prm_nwp_tend%ddt_v_turb, &
          & pt_patch%sync_cells_one_edge_in_domain, status=is_ready, scope=until_sync, &
          & name="prm_nwp_tend%ddt_v_turb")
          
        IF ( l_any_slowphys ) THEN
          z_ddt_u_tot_comm = new_icon_comm_variable(z_ddt_u_tot, &
            & pt_patch%sync_cells_one_edge_in_domain, &
            & status=is_ready, scope=until_sync, name="z_ddt_u_tot")
          z_ddt_v_tot_comm = new_icon_comm_variable(z_ddt_v_tot, &
            & pt_patch%sync_cells_one_edge_in_domain, &
            & status=is_ready, scope=until_sync, name="z_ddt_v_tot")
        ENDIF        
      ENDIF
      
       ! sync everything here
      CALL icon_comm_sync_all()

    ELSE
          
      IF ( l_any_slowphys .AND. lcall_phy_jg(itturb) ) THEN
          
        CALL sync_patch_array_mult(SYNC_C1, pt_patch, 4, z_ddt_u_tot, z_ddt_v_tot, &
                                 prm_nwp_tend%ddt_u_turb, prm_nwp_tend%ddt_v_turb)

      ELSE IF (lcall_phy_jg(itturb) ) THEN

        CALL sync_patch_array_mult(SYNC_C1, pt_patch, 2, prm_nwp_tend%ddt_u_turb, &
                                 prm_nwp_tend%ddt_v_turb)
      ENDIF
    ENDIF
    
    IF (timers_level > 3) CALL timer_stop(timer_phys_sync_ddt_u)
    !------------------------------------------------------------
    
      
    !------------------------------------------------------------
    ! compute on the halos
    IF (timers_level > 4) CALL timer_start(timer_phys_acc_par)
    IF (l_any_fastphys) THEN
      IF (my_process_is_mpi_all_parallel() ) THEN

        rl_start = min_rlcell_int-1
        rl_end   = min_rlcell 

        i_startblk = pt_patch%cells%start_blk(rl_start,1)
        i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE

        DO jb = i_startblk, i_endblk
          CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end )

          IF (lhdiff_rcf) THEN
            DO jk = 1, nlev
              DO jc =  i_startidx, i_endidx

                IF (p_metrics%mask_prog_halo_c(jc,jb)) THEN
                  pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                   &
                    &                     * pt_prog%rho(jc,jk,jb)*pt_diag%tempv(jc,jk,jb)))

                  pt_prog%theta_v(jc,jk,  jb) = pt_diag%tempv(jc,jk,jb) &
&                                             / pt_prog%exner(jc,jk,jb)

                  pt_prog%rhotheta_v(jc,jk,jb) = pt_prog%theta_v(jc,jk,jb) &
&                                              * pt_prog%rho(jc,jk,jb)
                ENDIF

              ENDDO
            ENDDO
          ELSE
            DO jk = 1, nlev
              DO jc =  i_startidx, i_endidx

                IF (p_metrics%mask_prog_halo_c(jc,jb)) THEN
                  pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                   &
                    &                     * pt_prog%rho(jc,jk,jb)*pt_diag%tempv(jc,jk,jb)))

                  pt_diag%exner_old(jc,jk,jb) = pt_diag%exner_old(jc,jk,jb) + &
                    pt_prog%exner(jc,jk,jb) - z_exner_sv(jc,jk,jb)

                  pt_prog%theta_v(jc,jk,  jb) = pt_diag%tempv(jc,jk,jb) &
&                                             / pt_prog%exner(jc,jk,jb)

                  pt_prog%rhotheta_v(jc,jk,jb) = pt_prog%theta_v(jc,jk,jb) &
&                                              * pt_prog%rho(jc,jk,jb)
                ENDIF

              ENDDO
            ENDDO
          ENDIF
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ENDIF ! my_process_is_mpi_all_parallel
    ENDIF ! fast-physics synchronization    
    IF (timers_level > 4) CALL timer_stop(timer_phys_acc_par)
    

    ! Initialize fields for runtime diagnostics
    ! In case that average ABS(dpsdt) is diagnosed
    IF (msg_level >= 12) THEN
      dps_blk(:)   = 0._wp
      npoints_blk(:) = 0
    ENDIF

    ! In case that maximum CFL is diagnosed
    IF (msg_level >= 13) THEN
      maxcfl(:) = 0._wp
    ENDIF

    ! In case that turbulence diagnostics are computed
    IF (msg_level >= 18) THEN
      maxtke(:,:)   = 0._wp
      maxtturb(:,:) = 0._wp
      maxuturb(:,:) = 0._wp
      maxvturb(:,:) = 0._wp
    ENDIF

    !-------------------------------------------------------------------------
    !>
    !!    @par Interpolation from  u,v onto v_n
    !!      ddt_vn_phy  =interpol(ddt_u_tot)+interpol(ddt_v_tot)
    !!      Calculate normal velocity at edge midpoints
    !-------------------------------------------------------------------------

    IF (timers_level > 4)  CALL timer_start(timer_phys_acc_2)
!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

    ! exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_e+1
    rl_end   = min_rledge_int

    i_startblk = pt_patch%edges%start_blk(rl_start,1)
    i_endblk   = pt_patch%edges%end_blk(rl_end,i_nchdom)


!$OMP DO PRIVATE(jb,jk,jce,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = i_startblk, i_endblk

      CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF ( l_any_slowphys .AND. lcall_phy_jg(itturb) ) THEN

#ifdef __LOOP_EXCHANGE
        DO jce = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=5
        DO jk = 1, nlev
          DO jce = i_startidx, i_endidx
#endif

            pt_diag%ddt_vn_phy(jce,jk,jb) =   pt_int_state%c_lin_e(jce,1,jb)           &
&                                 * ( z_ddt_u_tot(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
&                                   + z_ddt_v_tot(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
&                                                 + pt_int_state%c_lin_e(jce,2,jb)     &
&                                 * ( z_ddt_u_tot(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                    * pt_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
&                                  +  z_ddt_v_tot(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,2)%v2 )

            pt_prog%vn(jce,jk,jb) = pt_prog%vn(jce,jk,jb) + dtadv_loc * (              &
                                              pt_int_state%c_lin_e(jce,1,jb)           &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
&                       + prm_nwp_tend%ddt_v_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
&                                                 + pt_int_state%c_lin_e(jce,2,jb)     &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                    * pt_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
&                      +  prm_nwp_tend%ddt_v_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                  *  pt_patch%edges%primal_normal_cell(jce,jb,2)%v2 ) )

          ENDDO
        ENDDO

      ELSE IF (lcall_phy_jg(itturb) ) THEN
#ifdef __LOOP_EXCHANGE
        DO jce = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=8
        DO jk = 1, nlev
          DO jce = i_startidx, i_endidx
#endif

            pt_prog%vn(jce,jk,jb) = pt_prog%vn(jce,jk,jb) + dtadv_loc * (              &
                                              pt_int_state%c_lin_e(jce,1,jb)           &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
&                       + prm_nwp_tend%ddt_v_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
&                                                 + pt_int_state%c_lin_e(jce,2,jb)     &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                    * pt_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
&                      +  prm_nwp_tend%ddt_v_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                  *  pt_patch%edges%primal_normal_cell(jce,jb,2)%v2 ) )

          ENDDO
        ENDDO

      ENDIF

    ENDDO
!$OMP END DO


      ! Diagnosis of ABS(dpsdt) if msg_level >= 12
      IF (msg_level >= 12) THEN

        rl_start = grf_bdywidth_c+1
        rl_end   = min_rlcell_int

        i_startblk = pt_patch%cells%start_blk(rl_start,1)
        i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          DO jc = i_startidx, i_endidx
            ! Note: division by time step follows below
            dps_blk(jb) = dps_blk(jb) + &
              ABS(pt_diag%pres_sfc(jc,jb)-pt_diag%pres_sfc_old(jc,jb))
            npoints_blk(jb) = npoints_blk(jb) + 1
            pt_diag%pres_sfc_old(jc,jb) = pt_diag%pres_sfc(jc,jb)
          ENDDO
        ENDDO
!$OMP END DO
      ENDIF

      ! CFL-diagnostic if msg_level >= 13
      IF (msg_level >= 13) THEN

        rl_start = grf_bdywidth_c+1
        rl_end   = min_rlcell_int

        i_startblk = pt_patch%cells%start_blk(rl_start,1)
        i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,avg_invedgelen) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)


          DO jc = i_startidx, i_endidx
            avg_invedgelen(jc) = 3._wp/                                       &
              (pt_patch%edges%dual_edge_length(ieidx(jc,jb,1),ieblk(jc,jb,1))+&
               pt_patch%edges%dual_edge_length(ieidx(jc,jb,2),ieblk(jc,jb,2))+&
               pt_patch%edges%dual_edge_length(ieidx(jc,jb,3),ieblk(jc,jb,3)) )
          ENDDO

          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              maxcfl(jb) = MAX(maxcfl(jb),dt_loc*avg_invedgelen(jc)*( &
                SQRT(pt_diag%u(jc,jk,jb)**2+pt_diag%v(jc,jk,jb)**2)+  &
                SQRT(csfac*pt_diag%temp(jc,jk,jb)) ))
            ENDDO
          ENDDO

        ENDDO
!$OMP END DO

      ENDIF

      ! Extended turbulence diagnostics if msg_level >= 18
      IF (lcall_phy_jg(itturb) .AND. msg_level >= 18) THEN

        rl_start = grf_bdywidth_c+1
        rl_end   = min_rlcell_int

        i_startblk = pt_patch%cells%start_blk(rl_start,1)
        i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,avg_invedgelen) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)


          DO jk = 1, nlevp1
            DO jc = i_startidx, i_endidx
              maxtke(jb,jk) = MAX(maxtke(jb,jk),pt_prog_rcf%tke(jc,jk,jb))
            ENDDO
          ENDDO

          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              maxtturb(jb,jk) = MAX(maxtturb(jb,jk),ABS(prm_nwp_tend%ddt_temp_turb(jc,jk,jb)))
              maxuturb(jb,jk) = MAX(maxuturb(jb,jk),ABS(prm_nwp_tend%ddt_u_turb(jc,jk,jb)))
              maxvturb(jb,jk) = MAX(maxvturb(jb,jk),ABS(prm_nwp_tend%ddt_v_turb(jc,jk,jb)))
            ENDDO
          ENDDO

        ENDDO
!$OMP END DO NOWAIT

      ENDIF

!$OMP END PARALLEL

    IF (timers_level > 4) CALL timer_stop(timer_phys_acc_2)
    IF (timers_level > 3) CALL timer_start(timer_phys_sync_vn)
    IF (lcall_phy_jg(itturb)) CALL sync_patch_array(SYNC_E, pt_patch, pt_prog%vn)
    IF (timers_level > 3) CALL timer_stop(timer_phys_sync_vn)
    IF (timers_level > 2) CALL timer_stop(timer_phys_acc)

      ! dpsdt diagnostic - omitted in the case of a parallization test (p_test_run) because this
      ! is a purely diagnostic quantitiy, for which it does not make sense to implement an order-invariant
      ! summation
      IF (.NOT. p_test_run .AND. msg_level >= 12) THEN
        dpsdt_avg = SUM(dps_blk)
        npoints   = SUM(npoints_blk)
        dpsdt_avg = global_sum_array(dpsdt_avg)
        npoints   = global_sum_array(npoints)
        dpsdt_avg = dpsdt_avg/(REAL(npoints,wp)*dtadv_loc)
        ! Exclude initial time step where pres_sfc_old is zero
        IF (dpsdt_avg < 10000._wp/dtadv_loc) THEN
          WRITE(message_text,'(a,f12.6,a,i3)') 'average |dPS/dt| =',dpsdt_avg,' Pa/s in domain',jg
          CALL message('nwp_nh_interface: ', TRIM(message_text))
        ENDIF
      ENDIF

      IF (msg_level >= 13) THEN ! CFL diagnostic
        cflmax = MAXVAL(maxcfl)
        cflmax = global_max(cflmax) ! maximum over all PEs
        WRITE(message_text,'(a,f12.8,a,i2)') 'maximum horizontal CFL = ', cflmax, ' in domain ',jg
        CALL message('nwp_nh_interface: ', TRIM(message_text))

      ENDIF

      IF (msg_level >= 18 .AND. lcall_phy_jg(itturb)) THEN ! extended turbulence diagnostic
        DO jk = 1, nlevp1
          tkemax(jk) = MAXVAL(maxtke(:,jk))
        ENDDO
        DO jk = 1, nlev
          tturbmax(jk) = MAXVAL(maxtturb(:,jk))
          uturbmax(jk) = MAXVAL(maxuturb(:,jk))
          vturbmax(jk) = MAXVAL(maxvturb(:,jk))
        ENDDO

        ! Take maximum over all PEs
        tkemax   = global_max(tkemax)
        tturbmax = global_max(tturbmax)
        uturbmax = global_max(uturbmax)
        vturbmax = global_max(vturbmax)

        WRITE(message_text,'(a,i2)') 'Extended turbulence diagnostic for domain ',jg
        CALL message('nwp_nh_interface: ', TRIM(message_text))
        WRITE(message_text,'(a)') 'maximum TKE [m**2/s**2] and U,V,T-tendencies/s per level'
        CALL message('', TRIM(message_text))

        DO jk = 1, nlev
          WRITE(message_text,'(a,i3,4(a,e13.5))') 'level ',jk,': TKE =',tkemax(jk), &
            ', utend =',uturbmax(jk),', vtend =',vturbmax(jk),', ttend =',tturbmax(jk)
          CALL message('', TRIM(message_text))
        ENDDO
        jk = nlevp1
        WRITE(message_text,'(a,i3,a,e13.5)') 'level ',jk,': TKE =',tkemax(jk)
        CALL message('', TRIM(message_text))
      ENDIF

   
     CALL nwp_diagnosis(lcall_phy_jg,lredgrid,               & !input
                            & dt_phy_jg,p_sim_time,          & !input
                            & kstart_moist(jg),              & !input
                            & pt_patch, p_metrics,           & !input
                            & pt_prog, pt_prog_rcf,          & !in
                            & pt_diag,                       & !inout
                            & prm_diag,prm_nwp_tend)


    IF (ltimer) CALL timer_stop(timer_physics)


  END SUBROUTINE nwp_nh_interface

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------


  SUBROUTINE nh_update_prog_phy( pt_patch, pdtime, prm_nwp_tend, &
    &                            prm_diag, pt_prog_rcf )

    TYPE(t_patch),       INTENT(IN)   :: pt_patch     !!grid/patch info.
    TYPE(t_nwp_phy_tend),TARGET, INTENT(IN):: prm_nwp_tend   !< atm tend vars
    TYPE(t_nwp_phy_diag),INTENT(INOUT):: prm_diag     !!the physics variables
    TYPE(t_nh_prog),     INTENT(INOUT):: pt_prog_rcf  !!the tracer field at
                                                   !!reduced calling frequency
    REAL(wp),INTENT(in)            :: pdtime

    ! Local array bounds:

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !! blocks
    INTEGER :: i_startidx, i_endidx    !! slices
    INTEGER :: i_nchdom                !! domain index

    ! Local scalars:
    INTEGER  :: nlev        !< number of full levels
    INTEGER  :: jb          !block index
    INTEGER  :: jt          !tracers
    INTEGER  :: jk,jc,jg

    jg = pt_patch%id

    ! number of vertical levels
    nlev = pt_patch%nlev

    ! local variables related to the blocking

    i_nchdom  = MAX(1,pt_patch%n_childdom)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jt,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
           &                       i_startidx, i_endidx, rl_start, rl_end)

! KF fix to positive values
      DO jt=1, nqtendphy  ! qv,qc,qi
        DO jk = kstart_moist(jg), nlev
          DO jc = i_startidx, i_endidx

            pt_prog_rcf%tracer(jc,jk,jb,jt) =MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt)  &
              &                       + pdtime*prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,jt))
          ENDDO
        ENDDO
      ENDDO

!DR additional clipping for qr, qs 
!DR (very small negative values may occur during the transport process (order 10E-15)) 
      DO jt=iqr, iqs  ! qr,qs
        DO jk = kstart_moist(jg), nlev
          DO jc = i_startidx, i_endidx

            pt_prog_rcf%tracer(jc,jk,jb,jt) =MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt))
          ENDDO
        ENDDO
      ENDDO


      prm_diag%rain_con(i_startidx:i_endidx,jb) =                                       &
        &                                  prm_diag%rain_con(i_startidx:i_endidx,jb)    &
        &                                  + pdtime                                     &
        &                                  * prm_diag%rain_con_rate(i_startidx:i_endidx,jb)

      prm_diag%snow_con(i_startidx:i_endidx,jb) =                                       &
        &                                  prm_diag%snow_con(i_startidx:i_endidx,jb)    &
        &                                  + pdtime                                     &
        &                                  * prm_diag%snow_con_rate(i_startidx:i_endidx,jb)

      !for grid scale part: see mo_nwp_gscp_interface/nwp_microphysics
      prm_diag%tot_prec(i_startidx:i_endidx,jb) =                                       &
        &                              prm_diag%tot_prec(i_startidx:i_endidx,jb)        &
        &                              +  pdtime                                        &
        &                              * (prm_diag%rain_con_rate(i_startidx:i_endidx,jb)&
        &                              +  prm_diag%snow_con_rate(i_startidx:i_endidx,jb))

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SUBROUTINE nh_update_prog_phy

END MODULE mo_nh_interface_nwp

