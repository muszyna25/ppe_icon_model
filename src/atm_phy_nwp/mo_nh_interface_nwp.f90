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
MODULE mo_nh_interface_nwp

  USE mo_kind,               ONLY: wp

 ! USE mo_timer,              ONLY: timer_physics, timer_start, timer_stop, &
  USE mo_timer
  
  USE mo_exception,          ONLY: message, message_text !, finish
  USE mo_impl_constants,     ONLY: itconv, itccov, itrad, itgscp,         &
    &                              itsatad, itupdate, itturb, itsfc, itradheat, &
    &                              itsso, itgwd, icc,                           &
    &                              min_rlcell_int, min_rledge_int, min_rlcell
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_interpolation,      ONLY: rbf_vec_interpol_cell, edges2cells_scalar
  USE mo_grf_interpolation,  ONLY: t_gridref_state!,t_gridref_single_state, &
  USE mo_model_domain,       ONLY: t_patch
  USE mo_interpolation,      ONLY: t_int_state
  USE mo_nonhydro_state,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydrostatic_config, ONLY: kstart_moist, l_open_ubc
  USE mo_nwp_lnd_state,      ONLY: t_lnd_prog, t_lnd_diag!, t_lnd_state
  USE mo_ext_data,           ONLY: t_external_data
  USE mo_nwp_phy_state,      ONLY: t_nwp_phy_diag,&
                                 & t_nwp_phy_tend!, prm_diag
  USE mo_parallel_config,  ONLY: nproma, p_test_run
  USE mo_run_config,         ONLY: ntracer, iqv, iqc, iqi, &
       &                           iqr, iqs, msg_level, ltimer, timers_level
  USE mo_physical_constants, ONLY: rd, rd_o_cpd, vtmpc1, p0ref, cvd_o_rd 

  USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp

  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
  USE mo_cover_koe,          ONLY: cover_koe
  USE mo_satad,              ONLY: satad_v_3D
  USE mo_radiation,          ONLY: radheat, pre_radiation_nwp
 ! USE mo_sso_cosmo,          ONLY: sso
  USE mo_nwp_gw_interface,   ONLY: nwp_gwdrag 
  USE mo_nwp_gscp_interface, ONLY: nwp_microphysics
  USE mo_nwp_turb_interface, ONLY: nwp_turbulence
  USE mo_nwp_sfc_interface,  ONLY: nwp_surface
  USE mo_nwp_conv_interface, ONLY: nwp_convection
  USE mo_nwp_rad_interface,  ONLY: nwp_radiation
  USE mo_sync,               ONLY: sync_patch_array, sync_patch_array_mult, &
                                   SYNC_C, SYNC_C1
  USE mo_mpi,                ONLY: my_process_is_mpi_all_parallel
  USE mo_nwp_diagnosis,      ONLY: nwp_diagnosis
!  USE mo_communication,      ONLY: time_sync

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  REAL(wp), PARAMETER :: rd_o_p0ref = rd / p0ref
  REAL(wp), PARAMETER :: cpd_o_rd = 1._wp / rd_o_cpd

  PUBLIC :: nwp_nh_interface, nh_update_prog_phy

CONTAINS
  !
  !-----------------------------------------------------------------------
  !
  SUBROUTINE nwp_nh_interface(lcall_phy_jg,lredgrid,jstep,         & !input
                            & tcall_phy_jg,p_sim_time,             & !input
                            & pt_patch, pt_int_state,p_metrics,    & !input
                            & pt_par_patch, pt_par_int_state,      & !input
                            & pt_par_grf_state,                    & !input
                            & ext_data, mean_charlen,              & !input
                            & pt_prog,                             & !inout
                            & pt_prog_now_rcf, pt_prog_rcf,        & !in/inout
                            & pt_diag ,                            & !inout
                            & prm_diag,prm_nwp_tend,lnd_diag,      &
                            & lnd_prog_now, lnd_prog_new           ) !inout

    !>
    ! !INPUT PARAMETERS:

    LOGICAL, INTENT(IN)          ::   &             !< physics package time control (switches)
         &                          lcall_phy_jg(:) !< for domain jg
    LOGICAL, INTENT(IN)          :: lredgrid        !< use reduced grid for radiation
    INTEGER ,INTENT(in)          :: jstep
    REAL(wp),INTENT(in)          :: tcall_phy_jg(:) !< time interval for all physics
                                                    !< packages on domain jg
    REAL(wp),INTENT(in)          :: p_sim_time
    REAL(wp),INTENT(in),OPTIONAL :: mean_charlen !< characteristic griddistance, needed
                                                 !< by turbulence

    TYPE(t_patch),        TARGET,INTENT(in):: pt_patch     !<grid/patch info.
    TYPE(t_patch),        TARGET,INTENT(in):: pt_par_patch !<grid/patch info (parent grid)

    TYPE(t_int_state),    TARGET,INTENT(in):: pt_int_state      !< interpolation state
    TYPE(t_int_state),    TARGET,INTENT(in):: pt_par_int_state  !< " for parent grid
    TYPE(t_gridref_state),TARGET,INTENT(IN):: pt_par_grf_state  !< grid refinement state

    TYPE(t_nh_metrics)   ,       INTENT(in):: p_metrics
    TYPE(t_external_data),       INTENT(in):: ext_data

    TYPE(t_nh_diag), TARGET, INTENT(inout) :: pt_diag     !<the diagnostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog     !<the prognostic variables
    TYPE(t_nh_prog), TARGET, INTENT(IN)    :: pt_prog_now_rcf !<old state for tke
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog_rcf !<the prognostic variables (with
                                                          !< red. calling frequency for tracers!
    TYPE(t_nwp_phy_diag),       INTENT(inout) :: prm_diag
    TYPE(t_nwp_phy_tend),TARGET,INTENT(inout) :: prm_nwp_tend
    TYPE(t_lnd_prog),           INTENT(inout) :: lnd_prog_now, lnd_prog_new
    TYPE(t_lnd_diag),           INTENT(inout) :: lnd_diag

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

    INTEGER,  POINTER ::  iidx(:,:,:), iblk(:,:,:)

    REAL(wp):: &                                              !> temporal arrays for 
      & z_ddt_u_tot (nproma,pt_patch%nlev,pt_patch%nblks_c),& 
      & z_ddt_v_tot (nproma,pt_patch%nlev,pt_patch%nblks_c),& !< hor. wind tendencies
      & z_ddt_temp  (nproma,pt_patch%nlev,pt_patch%nblks_c)   !< Temperature tendency

    REAL(wp) :: &                                             !> virtual temperature
      &  z_aux_tempv(nproma,pt_patch%nlev  ,pt_patch%nblks_c)
 
    REAL(wp) :: z_exner_sv(nproma,pt_patch%nlev,pt_patch%nblks_c)

    !< vertical interfaces

    REAL(wp) :: z_airmass (nproma,pt_patch%nlev) !< needed for radheat
    REAL(wp) :: zi0       (nproma)     !< solar incoming radiation at TOA   [W/m2]
    REAL(wp) :: zsct ! solar constant (at time of year)
    REAL(wp) :: zcosmu0 (nproma,pt_patch%nblks_c)

    REAL(wp) :: rd_o_cvd

    REAL(wp)  z_qsum       !< summand of virtual increment
    REAL(wp)  z_ddt_qsum   !< summand of tendency of virtual increment

    !KF temporary field
    LOGICAL:: landseemask(nproma,pt_patch%nblks_c)

    REAL(wp) :: rcld(nproma,pt_patch%nlevp1)

!     write(0,*) "Entering nwp_nh_interface"
!     write(0,*) "========================="
    IF (ltimer) CALL timer_start(timer_physics)

    ! local variables related to the blocking

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    !define pointers
    iidx => pt_patch%edges%cell_idx
    iblk => pt_patch%edges%cell_blk

    rd_o_cvd  = 1._wp / cvd_o_rd

    !-------------------------------------------------------------------------
    !>  Update the tracer for every advective timestep,
    !!  all other updates are done in dynamics
    !-------------------------------------------------------------------------

    ! KF synchronize the update with satad

    IF (lcall_phy_jg(itupdate)) THEN

      IF (msg_level >= 12) &
           & CALL message('mo_nh_interface_nwp:', 'update_tracers')

      IF (timers_level > 2) CALL timer_start(timer_update_prog_phy)
      
      CALL nh_update_prog_phy(pt_patch              ,& !in
           &                  tcall_phy_jg(itupdate),& !in
           &                  pt_diag               ,& !in
           &                  prm_diag              ,& !inout phyfields 
           &                  pt_prog_rcf            )!inout tracer

      IF (timers_level > 2) CALL timer_stop(timer_update_prog_phy)
    
    ENDIF

    !!-------------------------------------------------------------------------
    !> ONLY saturation adjustment
    !!-------------------------------------------------------------------------

    IF (lcall_phy_jg(itsatad)) THEN

      IF (msg_level >= 12) &
           & CALL message('mo_nh_interface_nwp:', 'satad')


      IF (timers_level > 2) CALL timer_start(timer_diagnose_pres_temp)
        
      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf,    &
           &                              pt_diag, pt_patch,       &
           &                              opt_calc_temp=.TRUE.,    &
           &                              opt_calc_pres=.TRUE.,    &
           &                              opt_rlend=min_rlcell_int,&
           &                              opt_slev=kstart_moist(jg))

      IF (timers_level > 2) CALL timer_stop(timer_diagnose_pres_temp)

      IF (msg_level >= 15) THEN

        CALL message('mo_nh_interface_nwp:', 'before satad')

        WRITE(message_text,'(a,2E17.9)') ' max/min QV  = ',&
             & MAXVAL(pt_prog_rcf%tracer(:,:,:,iqv)), MINVAL(pt_prog_rcf%tracer (:,:,: ,iqv))
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,2E17.9)') ' max/min QC  = ',&
             & MAXVAL(pt_prog_rcf%tracer(:,:,:,iqc)), MINVAL(pt_prog_rcf%tracer (:,:,: ,iqc))
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,2E17.9)') ' max/min T  = ',&
             & MAXVAL(pt_diag%temp(:,:,:)), MINVAL(pt_diag%temp (:,:,:) )
        CALL message('', TRIM(message_text))

      ENDIF

      !in order to account for mesh refinement
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

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_qsum)
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
               & jdim     = 1                              ,& !> IN
               & kdim     = nlev                           ,& !> IN
               & ilo      = i_startidx                     ,& !> IN
               & iup      = i_endidx                       ,& !> IN
               & jlo      = 1                              ,& !> IN
               & jup      = 1                              ,& !> IN
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

!            z_qsum                = SUM(pt_prog_rcf%tracer (jc,jk,jb,iqc:iqs))
            z_qsum=   pt_prog_rcf%tracer (jc,jk,jb,iqc) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqi) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqr) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqs)
            
            z_aux_tempv(jc,jk,jb) =  pt_diag%temp(jc,jk,jb)                    &
              &                   * ( 1._wp +  vtmpc1                          &
              &                   * pt_prog_rcf%tracer(jc,jk,jb,iqv)  - z_qsum )

            pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                 &
              &                     * pt_prog%rho(jc,jk,jb)*z_aux_tempv(jc,jk,jb)))

            pt_prog%theta_v(jc,jk,jb) = z_aux_tempv  (jc,jk,jb) &
              &                       / pt_prog%exner(jc,jk,jb)

            pt_prog%rhotheta_v(jc,jk,jb) = pt_prog%theta_v(jc,jk,jb) &
              &                          * pt_prog%rho    (jc,jk,jb)

          ENDDO
        ENDDO
        IF (timers_level > 2) CALL timer_stop(timer_phys_exner)
      ENDDO ! nblks

!$OMP END DO
!$OMP END PARALLEL


      IF (msg_level >= 15) THEN

        CALL message('mo_nh_interface_nwp:', 'after satad')

        WRITE(message_text,'(a,2E17.9)') ' max/min QV  = ',&
             & MAXVAL(pt_prog_rcf%tracer(:,:,:,iqv)), MINVAL(pt_prog_rcf%tracer (:,:,: ,iqv) )
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,2E17.9)') ' max/min QC  = ',&
             & MAXVAL(pt_prog_rcf%tracer(:,:,:,iqc)), MINVAL(pt_prog_rcf%tracer (:,:,: ,iqc) )
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,2E17.9)') ' max/min T  = ',&
             & MAXVAL(pt_diag%temp(:,:,:)), MINVAL(pt_diag%temp (:,:,:) )
        CALL message('', TRIM(message_text))

      ENDIF

    ELSE ! satad turned off

!$OMP PARALLEL WORKSHARE
      ! Store exner function for sound-wave reduction and open upper boundary condition
      z_exner_sv(:,:,:) = pt_prog%exner(:,:,:)
!$OMP END PARALLEL WORKSHARE

    ENDIF ! satad

    IF ( lcall_phy_jg(itturb) .OR. lcall_phy_jg(itconv) ) THEN
    
      !-------------------------------------------------------------------------
      !>
      !!   Interpolation from v_n onto u,v =>  Reconstruct u and v
      !!   This is needed for turbulence and convection
      !!
      !-------------------------------------------------------------------------

      IF (msg_level >= 12) &
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

      IF (msg_level >= 15) THEN
        WRITE(message_text,'(a,2E17.9)') 'max/min vn  = ',&
             & MAXVAL (pt_prog%vn(:,:,:)), MINVAL(pt_prog%vn(:,:,:))
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,2E17.9)') 'max/min u  = ',&
             & MAXVAL (pt_diag%u(:,:,:)), MINVAL(pt_diag%u(:,:,:))
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,2E17.9)') 'max/min v  = ',&
             & MAXVAL (pt_diag%v(:,:,:)), MINVAL(pt_diag%v(:,:,:) )
        CALL message('', TRIM(message_text))
      ENDIF

    ENDIF ! diagnose u/v


    !!-------------------------------------------------------------------------
    !>  turbulent transfer and diffusion  and microphysics
    !!
    !!  Because we consider  the followings physical processes as fast ones
    !!  we allow her the update of prognostic variables inside the subroutines
    !!  This means that the conversion back to the ICON-prognostic variables
    !!  has to be done atferwards
    !!-------------------------------------------------------------------------

    IF (lcall_phy_jg(itgscp) .OR. lcall_phy_jg(itturb)) THEN

      IF (msg_level >= 12) &
        & CALL message('mo_nh_interface_nwp:', 'diagnose pres/temp for fast physics')

      !-------------------------------------------------------------------------
      !> diagnose
      !!    - pressure on main and intermediate levels    =>  opt_calc_pres=.TRUE.
      !!    - temperature on main and intermediate levels =>  opt_calc_temp=.TRUE.
      !!    - virtual temperature on main levels          =>  opt_calc_tempv=.TRUE.
      !!    out of prognostic variables
      !!
      !-------------------------------------------------------------------------

      IF (timers_level > 2) CALL timer_start(timer_diagnose_pres_temp)
      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf, &
        & pt_diag, pt_patch,       &
        & opt_calc_temp =.TRUE.,   &
        & opt_calc_pres =.TRUE.,   &
        & opt_calc_tempv=.TRUE.,   &
        & opt_rlend=min_rlcell_int )
      IF (timers_level > 2) CALL timer_stop(timer_diagnose_pres_temp)
    
    ENDIF 

    IF (  lcall_phy_jg(itturb)) THEN

      IF (timers_level > 1) CALL timer_start(timer_nwp_turbulence)
      CALL nwp_turbulence (  tcall_phy_jg(itturb),              & !>input
                            & pt_patch, p_metrics,              & !>input
                            & ext_data, mean_charlen,           & !>input
                            & pt_prog,                          & !>inout
                            & pt_prog_now_rcf, pt_prog_rcf,     & !>in/inout
                            & pt_diag ,                         & !>inout
                            & prm_diag,prm_nwp_tend,            & !>inout 
                            & lnd_prog_now,lnd_diag             )!>inout
      IF (timers_level > 1) CALL timer_stop(timer_nwp_turbulence)
    ENDIF !lcall(itturb)

    !-------------------------------------------------------------------------
    !  prognostic microphysic and precipitation scheme
    !-------------------------------------------------------------------------
    IF ( lcall_phy_jg(itgscp)) THEN

      IF (msg_level >= 12) &
        & CALL message('mo_nh_interface_nwp:', 'microphysics')

      !> temperture and tracers have been update by turbulence therefore
      !! no update from prognostic variables is needed. virtual temperature
      !! no INPUT in microphysics

      IF (timers_level > 2) CALL timer_start(timer_diagnose_pres_temp)
      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf, &
                            &  pt_diag, pt_patch,               &
                            &  opt_calc_temp =.FALSE.,          &
                            &  opt_calc_pres =.TRUE.,           &
                            &  opt_calc_tempv=.FALSE.,          &
                            &  opt_rlend=min_rlcell_int         )
      IF (timers_level > 2) CALL timer_stop(timer_diagnose_pres_temp)

      IF (timers_level > 1) CALL timer_start(timer_nwp_microphysics)
      CALL nwp_microphysics ( tcall_phy_jg(itgscp),             & !>input
                            & pt_patch, p_metrics,              & !>input
                            & pt_prog,                          & !>inout
                            & pt_prog_rcf,                      & !>inout
                            & pt_diag ,                         & !>inout
                            & prm_diag                          ) !>inout 
      IF (timers_level > 1) CALL timer_stop(timer_nwp_microphysics)

    ENDIF


    ! Synchronize tracers if any of the updating processes was active
    IF (timers_level > 1) CALL timer_start(timer_phys_sync_patch)
    CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer, f4din=pt_prog_rcf%tracer, &
                               lpart4d=.TRUE.)
    IF (timers_level > 1) CALL timer_stop(timer_phys_sync_patch)


    IF (lcall_phy_jg(itsatad) .OR. lcall_phy_jg(itgscp) .OR. &
      & lcall_phy_jg(itturb)) THEN
      
      IF (timers_level > 1) CALL timer_start(timer_fast_phys)

      ! Remark: in the (unusual) case that satad is used without any other physics,
      ! recalculation of the thermodynamic variables is duplicated here. However,
      ! this is the easiest way to combine minimization of halo communications
      ! with a failsafe flow control

      IF (msg_level >= 12) &
        & CALL message('mo_nh_interface_nwp:', 'recalculate thermodynamic variables')

      IF (p_test_run) z_aux_tempv(:,:,:) = 0._wp

      !in order to account for mesh refinement
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx, i_endidx, z_qsum )

      DO jb = i_startblk, i_endblk
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end )

        pt_diag%thermal_exp_fastphy(:,jb) = 0._wp

        !-------------------------------------------------------------------------
        !>
        !! Re-calculate scalar prognostic variables out of physics variables!
        !!
        !-------------------------------------------------------------------------

        DO jk = 1, nlev
          DO jc =  i_startidx, i_endidx

            z_qsum=   pt_prog_rcf%tracer (jc,jk,jb,iqc) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqi) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqr) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqs)
            z_aux_tempv(jc,jk,jb) =  pt_diag%temp(jc,jk,jb)            &
&                                  * ( 1._wp +  vtmpc1                 &
&                                  *  pt_prog_rcf%tracer(jc,jk,jb,iqv) &
&                                   - z_qsum )

            pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                 &
              &                     * pt_prog%rho(jc,jk,jb)*z_aux_tempv(jc,jk,jb)))

            pt_diag%exner_old(jc,jk,jb) = pt_diag%exner_old(jc,jk,jb) + &
              pt_prog%exner(jc,jk,jb) - z_exner_sv(jc,jk,jb)

            pt_diag%exner_fphy_incr(jc,jk,jb) =                        &
              -0.5_wp * (pt_prog%exner(jc,jk,jb) - z_exner_sv(jc,jk,jb))

            pt_prog%theta_v(jc,jk,  jb) = z_aux_tempv  (jc,jk,jb) &
&                                       / pt_prog%exner(jc,jk,jb)

            pt_prog%rhotheta_v(jc,jk,jb) = pt_prog%theta_v(jc,jk,jb) &
&                                        * pt_prog%rho(jc,jk,jb)
          ENDDO
        ENDDO

        IF (l_open_ubc) THEN
          ! vertically integrated isobaric expansion needed for open upper boundary condition
          DO jk = 1, nlev
            DO jc =  i_startidx, i_endidx
              pt_diag%thermal_exp_fastphy(jc,jb) = pt_diag%thermal_exp_fastphy(jc,jb) &
                & + cpd_o_rd*(pt_prog%exner(jc,jk,jb) - z_exner_sv(jc,jk,jb))         &
                & / (tcall_phy_jg(itsatad)*pt_prog%exner(jc,jk,jb))                   &
                & * p_metrics%ddqz_z_full(jc,jk,jb)
            ENDDO
          ENDDO
        ENDIF

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      ! Synchronize z_aux_tempv, then recompute thermodynamic variables on halo points
      ! (This is more efficient than synchronizing three variables)
      CALL sync_patch_array(SYNC_C, pt_patch, z_aux_tempv)

      IF (my_process_is_mpi_all_parallel() ) THEN
        rl_start = min_rlcell_int-1
        rl_end   = min_rlcell 

        i_startblk = pt_patch%cells%start_blk(rl_start,1)
        i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx, i_endidx)

        DO jb = i_startblk, i_endblk
          CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end )

          DO jk = 1, nlev
            DO jc =  i_startidx, i_endidx

              IF (p_metrics%mask_prog_halo_c(jc,jb)) THEN
                pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                 &
                  &                     * pt_prog%rho(jc,jk,jb)*z_aux_tempv(jc,jk,jb)))

                pt_diag%exner_old(jc,jk,jb) = pt_diag%exner_old(jc,jk,jb) + &
                  pt_prog%exner(jc,jk,jb) - z_exner_sv(jc,jk,jb)

                pt_diag%exner_fphy_incr(jc,jk,jb) =                        &
                  -0.5_wp * (pt_prog%exner(jc,jk,jb) - z_exner_sv(jc,jk,jb))

                pt_prog%theta_v(jc,jk,  jb) = z_aux_tempv  (jc,jk,jb) &
&                                           / pt_prog%exner(jc,jk,jb)

                pt_prog%rhotheta_v(jc,jk,jb) = pt_prog%theta_v(jc,jk,jb) &
&                                            * pt_prog%rho(jc,jk,jb)
              ENDIF

            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
!$OMP END PARALLEL

      ENDIF


    IF (lcall_phy_jg(itturb)) THEN

      IF (msg_level >= 12) &
        & CALL message('mo_nh_interface_nwp:', 'diagnose pres/temp for fast physics')

      !-------------------------------------------------------------------------
      !> diagnose
      !!    - pressure on main levels    =>  opt_calc_pres=.TRUE.
      !!    - temperature on main levels =>  opt_calc_temp=.TRUE.
      !!    - pressure and temperature on interface levels
      !!
      !!    out of prognostic variables
      !!
      !-------------------------------------------------------------------------

      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf, &
        & pt_diag, pt_patch,       &
        & opt_calc_temp =.TRUE.,   &
        & opt_calc_pres =.TRUE.,   &
        & opt_calc_tempv=.TRUE.,   &
        & opt_rlend=min_rlcell_int )

    ENDIF 


    IF (  lcall_phy_jg(itsfc)) THEN

      CALL nwp_surface    (  tcall_phy_jg(itsfc), jstep,        & !>input
                            & pt_patch,                         & !>input
                            & ext_data,                         & !>input
                            & pt_prog_rcf,     & !>in/inout rcf=reduced calling freq.
                            & pt_diag ,                         & !>inout
                            & prm_diag,                         & !>inout 
                            & lnd_prog_now, lnd_prog_new,       & !>inout
                            & lnd_diag                          ) !>inout
    ENDIF !lcall(itsfc)

    IF (timers_level > 1) CALL timer_stop(timer_fast_phys)
   
   ENDIF ! end of fast physics part

    !!-------------------------------------------------------------------------
    !!  slow physics part
    !!-------------------------------------------------------------------------


    IF (lcall_phy_jg(itrad) .OR.  lcall_phy_jg(itconv) &
&     .OR. lcall_phy_jg(itccov) .OR. lcall_phy_jg(itsso) &
&                                 .OR. lcall_phy_jg(itgwd)) THEN

      IF (msg_level >= 12) &
&           CALL message('mo_nh_interface', 'diagnose pres/temp for slow physics')

      !-------------------------------------------------------------------------
      !> diagnose
      !!    - pressure on main levels   
      !!    - pressure and temperature on interface levels
      !!    - pressure thickness   
      !!    - temperature at interface levels
      !!                                 =>  opt_calc_pres=.TRUE.
      !!    - temperature on main levels =>  opt_calc_temp=.TRUE.
      !!
      !! out of prognostic variables
      !!
      !-------------------------------------------------------------------------

      IF (timers_level > 2) CALL timer_start(timer_diagnose_pres_temp)
      IF ( .NOT. (lcall_phy_jg(itgscp) .OR. lcall_phy_jg(itturb))) THEN

        ! If slow physics is called without fast physics (which should happen
        ! at the initial time step only), temperature and pressure need to be calculated
        CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf, &
          &                               pt_diag, pt_patch,      &
          &                               opt_calc_temp=.TRUE.,   &
          &                               opt_calc_pres=.TRUE.,   &
          &                               opt_rlend=min_rlcell_int)


      ELSE ! diagnose only pressure because temperature is up to date

        CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf, &
          &                               pt_diag, pt_patch,      &
          &                               opt_calc_temp=.FALSE.,  &
          &                               opt_calc_pres=.TRUE.,   &
          &                               opt_rlend=min_rlcell_int)


      ENDIF
      IF (timers_level > 2) CALL timer_stop(timer_diagnose_pres_temp)
    ENDIF ! checking slow physics


    IF ( lcall_phy_jg(itconv)  ) THEN

      IF (msg_level >= 12) &
&           CALL message('mo_nh_interface', 'convection')

      IF (timers_level > 2) CALL timer_start(timer_nwp_convection)
      CALL nwp_convection (  tcall_phy_jg(itconv),              & !>input
                            & pt_patch, p_metrics,              & !>input
!                            & ext_data,                         & !>input
                            & pt_prog,                          & !>input
                            & pt_prog_rcf,                      & !>input
                            & pt_diag ,                         & !>inout
                            & prm_diag,prm_nwp_tend             ) !>inout 
      IF (timers_level > 2) CALL timer_stop(timer_nwp_convection)

    ENDIF! convection

    ! Note: computation of cloud cover is now forced to be synchronized with radiation
    ! because the output of this routine is not used anywhere else
    ! It is also important to call cover_koe even in the case of inwp_cldcover=0 because
    ! qv_tot is zero otherwise
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

      IF (msg_level >= 12) &
&           CALL message('mo_nh_interface', 'cloud cover')

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
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,rcld),SCHEDULE(guided)

      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)

          rcld = 0.0_wp ! standard deviation of saturation deficit=0 for now, needs to be specified form turbulence
          landseemask(:,jb)   = .FALSE. ! has to come from external parameters later on!!!

!#ifdef __BOUNDCHECK
          IF (timers_level > 2) CALL timer_start(timer_cover_koe)
          CALL cover_koe &
&             (kidia  = i_startidx ,   kfdia  = i_endidx  ,       & !! in:  horizonal begin, end indices
&              klon = nproma,  kstart = kstart_moist(jg)  ,       & !! in:  horiz. and vert. vector length
&              klev   = nlev,                                     &
&              icldscheme = atm_phy_nwp_config(jg)%inwp_cldcover ,& !! in:  cloud cover option
&              tt     = pt_diag%temp         (:,:,jb)     ,       & !! in:  temperature at full levels
&              pp     = pt_diag%pres         (:,:,jb)     ,       & !! in:  pressure at full levels
&              ps     = pt_diag%pres_sfc     (:,jb)       ,       & !! in:  surface pressure at full levels
&              t_g    = lnd_prog_now%t_g     (:,jb)       ,       & !! in:  surface temperature
&              pgeo   = p_metrics%geopot_agl (:,:,jb)     ,       & !! in:  geopotential height
&              rho    = pt_prog%rho          (:,:,jb  )   ,       & !! in:  density
&              rcld   = rcld                              ,       & !! in:  standard deviation of saturation deficit
&              ldland = landseemask          (:,jb)       ,       & !! in:  land/sea mask
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
  
!$OMP END DO
!$OMP END PARALLEL

    ENDIF! cloud cover

    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------

    IF ( lcall_phy_jg(itrad) ) THEN

      ! To be on the safe side and to avoid confusion during development,
      ! do the diagnosis of t and p also here.
      ! May be removed when final state is reached.
      IF (timers_level > 2) CALL timer_start(timer_diagnose_pres_temp)
      IF ( atm_phy_nwp_config(jg)%inwp_radiation == 1 ) THEN
        CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf, &
          &                               pt_diag, pt_patch,      &
          &                               opt_calc_temp =.TRUE.,  &
          &                               opt_calc_tempv=.FALSE., &
          &                               opt_calc_pres =.TRUE.,  &
          &                               opt_rlend=min_rlcell_int)

      ELSEIF ( atm_phy_nwp_config(jg)%inwp_radiation == 2 ) THEN
        CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf,   &
          &                               pt_diag, pt_patch,        &
          &                               opt_calc_temp =.TRUE.,    &
          &                               opt_calc_tempv=.FALSE.,   &
          &                               opt_calc_pres =.TRUE.,    &
          &                               lnd_prog = lnd_prog_now,  &
          &                               opt_calc_temp_ifc =.TRUE.,&
          &                               opt_rlend=min_rlcell_int  )
 
      ENDIF
      IF (timers_level > 2) CALL timer_stop(timer_diagnose_pres_temp)
    
      IF (ltimer) CALL timer_start(timer_nwp_radiation)
      CALL nwp_radiation (lredgrid,p_sim_time,   & ! in
           &              pt_patch,pt_par_patch, & ! in
           &              pt_par_int_state,      & ! in
           &              pt_par_grf_state,      & ! in
           &              ext_data,              & ! in
           &              lnd_diag,              & ! in
           &              pt_prog,pt_prog_rcf,   & ! inout
           &              pt_diag,prm_diag,      & ! inout
           &              lnd_prog_now           ) ! in
      IF (ltimer) CALL timer_stop(timer_nwp_radiation)
     
    ENDIF


!GZ: The radheat call has to be done every physics time step once it has reached its final state!!
!    (missing items: computation of surface radiation balance, perhaps also linearized 
!    update of downward LW radiation depending on change of sfc temp since last full rad call)

    IF ( lcall_phy_jg(itradheat) ) THEN

      IF (msg_level >= 12) &
&           CALL message('mo_nh_interface', 'radiative heating')


      IF (timers_level > 1) CALL timer_start(timer_pre_radiation_nwp)
      
      CALL pre_radiation_nwp (                       &
        & kbdim      = nproma,                       &
        & p_inc_rad  = tcall_phy_jg(itradheat),      &
        & p_sim_time = p_sim_time,                   &
        & pt_patch   = pt_patch,                     &
        & zsmu0      = zcosmu0,                      &
        & zsct       = zsct )
      IF (timers_level > 1) CALL timer_stop(timer_pre_radiation_nwp)
      

    !in order to account for mesh refinement
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)
      
    IF (timers_level > 2) CALL timer_start(timer_radheat)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,zi0,z_airmass)
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

!test KF
        prm_diag%swflxsfc (:,jb)=0._wp
        prm_diag%lwflxsfc (:,jb)=0._wp
        prm_diag%swflxtoa (:,jb)=0._wp

!#ifdef __BOUNDCHECK
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
        & pmair=z_airmass                        ,&! in     layer air mass             [kg/m2]
        & pq=pt_prog_rcf%tracer(:,:,jb,iqv)      ,&! in     specific moisture           [kg/kg]
        & pi0=zi0                                ,&! in     solar incoming flux at TOA  [W/m2]
        & pemiss=ext_data%atm%emis_rad(:,jb)     ,&! in     lw sfc emissivity
        & ptsfc=lnd_prog_now%t_g(:,jb)           ,&! in     surface temperature         [K]
        & ptsfctrad=prm_diag%tsfctrad(:,jb)      ,&! in     sfc temp. used for pflxlw   [K]
        & ptrmsw=prm_diag%trsolall (:,:,jb)      ,&! in     shortwave net tranmissivity []
        & pflxlw=prm_diag%lwflxall (:,:,jb)      ,&! in     longwave net flux           [W/m2]
        !
        ! output
        ! ------
        !
        & pdtdtradsw=prm_nwp_tend%ddt_temp_radsw(:,:,jb),&! out    rad. heating by SW        [K/s]
        & pdtdtradlw=prm_nwp_tend%ddt_temp_radlw(:,:,jb),&! out    rad. heating by lw        [K/s]
        & pflxsfcsw =prm_diag%swflxsfc (:,jb)   ,&        ! out shortwave surface net flux [W/m2]
        & pflxsfclw =prm_diag%lwflxsfc (:,jb)   ,&        ! out longwave surface net flux  [W/m2]
        & pflxtoasw =prm_diag%swflxtoa (:,jb) )           ! out shortwave toa net flux     [W/m2]

        IF ( p_sim_time > 1.e-1_wp ) THEN

         !sum up for averaged fluxes
          !T.R.: this is not correct for output after 1st timestep,
          !e.g. tcall_phy_jg(itradheat) may then be greater than p_sim_time
          !leading to wrong averaging.
         DO jc =  i_startidx, i_endidx
          prm_diag%swflxsfc_avg(jc,jb) = ( prm_diag%swflxsfc_avg(jc,jb)                   &
                                 &  * (p_sim_time - tcall_phy_jg(itradheat))              &
                                 &  + tcall_phy_jg(itradheat) * prm_diag%swflxsfc(jc,jb)) &
                                 &  / p_sim_time
          prm_diag%lwflxsfc_avg(jc,jb) = ( prm_diag%lwflxsfc_avg(jc,jb)                   &
                                 &  * (p_sim_time - tcall_phy_jg(itradheat))              &
                                 &  + tcall_phy_jg(itradheat) * prm_diag%lwflxsfc(jc,jb)) &
                                 &  / p_sim_time
          prm_diag%swflxtoa_avg(jc,jb) = ( prm_diag%swflxtoa_avg(jc,jb)                   &
                                 &  * (p_sim_time - tcall_phy_jg(itradheat))              &
                                 &  + tcall_phy_jg(itradheat) * prm_diag%swflxtoa(jc,jb)) &
                                 &  / p_sim_time
          prm_diag%lwflxtoa_avg(jc,jb) = ( prm_diag%lwflxtoa_avg(jc,jb)                   &
                                 &  * (p_sim_time - tcall_phy_jg(itradheat))              &
                                &  + tcall_phy_jg(itradheat) * prm_diag%lwflxall(jc,1,jb)) &
                                &  / p_sim_time
         ENDDO

        END IF

      ENDDO ! blocks

!$OMP END DO
!$OMP END PARALLEL

      IF (timers_level > 2) CALL timer_stop(timer_radheat)
      
      IF (msg_level >= 14) THEN
        WRITE(message_text,'(a,2E17.9)') 'max/min SW transmissivity  = ',&
             & MAXVAL (prm_diag%trsolall(:,:,:)), MINVAL(prm_diag%trsolall(:,:,:))
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,2E17.9)') 'max/min LW net flux  = ',&
             & MAXVAL (prm_diag%lwflxall(:,:,:)), MINVAL(prm_diag%lwflxall(:,:,:))
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,2E17.9)') 'max/min dTdt sw= ',&
          & MAXVAL (prm_nwp_tend%ddt_temp_radsw(:,:,:)),MINVAL(prm_nwp_tend%ddt_temp_radsw(:,:,:))
        CALL message('', TRIM(message_text))
         WRITE(message_text,'(a,2E17.9)') 'max/min dTdt lw= ',&
          & MAXVAL (prm_nwp_tend%ddt_temp_radlw(:,:,:)),MINVAL(prm_nwp_tend%ddt_temp_radlw(:,:,:))
        CALL message('', TRIM(message_text))
      ENDIF
      
    ENDIF  ! inwp_radiation

    !-------------------------------------------------------------------------
    !  Gravity waves drag: orographic and non-orographic
    !-------------------------------------------------------------------------

    IF (lcall_phy_jg(itsso) .OR. lcall_phy_jg(itgwd)) THEN

       IF (msg_level >= 12) &
&           CALL message('mo_nh_interface', 'gravity waves')

        IF (timers_level > 3) CALL timer_start(timer_sso)

       CALL nwp_gwdrag  (  tcall_phy_jg(itsso),             & !>input
                          &tcall_phy_jg(itgwd),             & !> input
                            &   pt_patch,p_metrics,         & !>input
                            &   ext_data,                   & !>input
                            &   pt_prog,                    & !>in
                            &   pt_diag ,                   & !>inout
                            &   prm_diag,prm_nwp_tend      ) !>inout

    IF (timers_level > 3) CALL timer_stop(timer_sso)
    ENDIF ! inwp_sso
     
    IF (ltimer) CALL timer_start(timer_physic_acc)
    !-------------------------------------------------------------------------
    !>  accumulate scalar tendencies of slow_physics
    !-------------------------------------------------------------------------
    IF (lcall_phy_jg(itradheat) .OR. lcall_phy_jg(itconv) & 
&       .OR. lcall_phy_jg(itccov).OR. lcall_phy_jg(itsso)) THEN

      IF (timers_level > 2) CALL timer_start(timer_physic_acc_1)


      !in order to account for mesh refinement
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)
      
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx, i_endidx , z_qsum, z_ddt_qsum)
!
      DO jb = i_startblk, i_endblk
!
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)

        z_ddt_temp(i_startidx:i_endidx,:,jb) =                                               &
   &                                       prm_nwp_tend%ddt_temp_radsw(i_startidx:i_endidx,:,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_radlw(i_startidx:i_endidx,:,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_gwd  (i_startidx:i_endidx,:,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_sso  (i_startidx:i_endidx,:,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_pconv(i_startidx:i_endidx,:,jb)

!   microphysics and turbulent increments are already updated

        pt_diag%ddt_tracer_phy(i_startidx: i_endidx,kstart_moist(jg):,jb,iqv) =       &
          & prm_nwp_tend%ddt_tracer_pconv(i_startidx:i_endidx,kstart_moist(jg):,jb,iqv)

        pt_diag%ddt_tracer_phy(i_startidx: i_endidx,kstart_moist(jg):,jb,iqc) =       &
          & prm_nwp_tend%ddt_tracer_pconv(i_startidx:i_endidx,kstart_moist(jg):,jb,iqc)

        pt_diag%ddt_tracer_phy(i_startidx: i_endidx,kstart_moist(jg):,jb,iqi) =       &
          & prm_nwp_tend%ddt_tracer_pconv(i_startidx:i_endidx,kstart_moist(jg):,jb,iqi)

!        pt_diag%ddt_tracer_phy(i_startidx: i_endidx,kstart_moist(jg):,jb,iqs) =       &
!          & prm_nwp_tend%ddt_tracer_pconv(i_startidx:i_endidx,kstart_moist(jg):,jb,iqs)

!------------------------
!>
!! @par convert temperature tendencies into exner tendencies
!! Since the exner function shows up as @f$\Pi=\frac{T_v}{\theta_v} @f$ this relates
!! to pressure and virtual temperature tendencies
!! @f$ \frac{d \pi}{d t} = \frac{1}{c_{pd} \theta_v \rho} \, \frac{dp}{d t} @f$
!!
!! @f$ \frac{dp}{d t} = (c_p/c_v -1) Q_h + c_p/c_v Q_h @f$ ,
!!
!!  where @f$ Q_h = \frac{d T}{d t} |_{phys} @f$
!!
!! and @f$ Q_m = R_d \,T \,\rho \frac{d \alpha}{d t} @f$.
!!
!! The resulting tendency can be written as
!!
!! @f$ \frac{d \pi}{d t}= \frac{R}{c_{v} \theta_v} \left( \frac{d T}{dt} \\
!!                        + T \frac{d \alpha}{ d t} \right) @f$
!!
!--------------------------------------------------------------------------------

        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx

!            z_qsum                   = SUM(pt_prog_rcf%tracer    (jc,jk,jb,iqc:iqs))
            z_qsum=   pt_prog_rcf%tracer (jc,jk,jb,iqc) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqi) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqr) &
              &     + pt_prog_rcf%tracer (jc,jk,jb,iqs)
            
!            z_ddt_qsum               = SUM(pt_diag%ddt_tracer_phy(jc,jk,jb,iqc:iqs))
            z_ddt_qsum =   pt_diag%ddt_tracer_phy(jc,jk,jb,iqc) &
              &          + pt_diag%ddt_tracer_phy(jc,jk,jb,iqi) &
              &          + pt_diag%ddt_tracer_phy(jc,jk,jb,iqr) &
              &          + pt_diag%ddt_tracer_phy(jc,jk,jb,iqs)

            pt_diag%ddt_exner_phy(jc,jk,jb) = rd_o_cvd / pt_prog%theta_v(jc,jk,jb)           &
              &                             * (z_ddt_temp(jc,jk,jb)                          &
              &                             *(1._wp + vtmpc1*pt_prog_rcf%tracer(jc,jk,jb,iqv)&
              &                                - z_qsum) + pt_diag%temp(jc,jk,jb)             &
              &                 * (vtmpc1 * pt_diag%ddt_tracer_phy(jc,jk,jb,iqv)-z_ddt_qsum ))

          ENDDO
        ENDDO
      ENDDO !blocks

!$OMP END DO
!$OMP END PARALLEL


      IF (msg_level >= 15) THEN

        WRITE(message_text,'(a,2E17.9)') ' max/min z_ddt_temp = ',&
          & MAXVAL( z_ddt_temp(:,:,:)),  &
          & MINVAL( z_ddt_temp(:,:,:))
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,2E17.9)') ' max/min z_ddt_qv  = ',&
          & MAXVAL( prm_nwp_tend%ddt_tracer_pconv(:,:,:,iqv)),  &
          & MINVAL( prm_nwp_tend%ddt_tracer_pconv(:,:,:,iqv))
        CALL message('', TRIM(message_text))

       WRITE(message_text,'(a,2(E17.9,2x))') 'max/min ddt_exner  = ',&
          & MAXVAL (pt_diag%ddt_exner_phy(:,:,:) ), MINVAL(pt_diag%ddt_exner_phy(:,:,:)  )
        CALL message('', TRIM(message_text))
      ENDIF

      IF (timers_level > 2) CALL timer_stop(timer_physic_acc_1)
    ENDIF ! slow physics tendency accumulation

  !--------------------------------------------------------

    !-------------------------------------------------------------------------
    !>  accumulate vector tendencies of all physics
    !-------------------------------------------------------------------------

    IF ( lcall_phy_jg(itconv) .OR. lcall_phy_jg(itturb) ) THEN

      IF (p_test_run) THEN
        z_ddt_u_tot = 0._wp
        z_ddt_v_tot = 0._wp
      ENDIF

      IF (ltimer) CALL timer_start(timer_physic_acc_2)
!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx, i_endidx)
!
      DO jb = i_startblk, i_endblk
!
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)

        z_ddt_u_tot(i_startidx:i_endidx,:,jb) =                                            &
   &                                         prm_nwp_tend%ddt_u_turb (i_startidx:i_endidx,:,jb)&
   &                                       + prm_nwp_tend%ddt_u_gwd  (i_startidx:i_endidx,:,jb)&
   &                                       + prm_nwp_tend%ddt_u_sso  (i_startidx:i_endidx,:,jb)&
   &                                       + prm_nwp_tend%ddt_u_pconv(i_startidx:i_endidx,:,jb)

        z_ddt_v_tot(i_startidx:i_endidx,:,jb) =                                            &
   &                                         prm_nwp_tend%ddt_v_turb (i_startidx:i_endidx,:,jb)&
   &                                       + prm_nwp_tend%ddt_v_gwd  (i_startidx:i_endidx,:,jb)&
   &                                       + prm_nwp_tend%ddt_v_sso  (i_startidx:i_endidx,:,jb)&
   &                                       + prm_nwp_tend%ddt_v_pconv(i_startidx:i_endidx,:,jb)
        ENDDO
!$OMP END DO
!$OMP END PARALLEL

    CALL sync_patch_array_mult(SYNC_C1, pt_patch, 2, z_ddt_u_tot, z_ddt_v_tot)

      !-------------------------------------------------------------------------
      !>
      !!    @par Interpolation from  u,v onto v_n
      !!      ddt_vn_phy  =interpol(ddt_u_tot)+interpol(ddt_v_tot)
      !!      Calculate normal velocity at edge midpoints
      !-------------------------------------------------------------------------

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

      !in order to account for mesh refinement
      rl_start = grf_bdywidth_e+1
      rl_end   = min_rledge_int

      i_startblk = pt_patch%edges%start_blk(rl_start,1)
      i_endblk   = pt_patch%edges%end_blk(rl_end,i_nchdom)


!$OMP DO PRIVATE(jb,jk,jce,i_startidx,i_endidx)

      DO jb = i_startblk, i_endblk

        CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
&                            i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
        DO jce = i_startidx, i_endidx
          DO jk = 1, nlev
#else
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
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      IF (ltimer) CALL timer_stop(timer_physic_acc_2)

      IF (msg_level >= 15) THEN
        WRITE(message_text,'(a,2(E17.9,2x))') 'max/min ddt_vn  = ',&
             & MAXVAL (pt_diag%ddt_vn_phy(:,:,:) ), MINVAL(pt_diag%ddt_vn_phy(:,:,:)  )
        CALL message('', TRIM(message_text))
      ENDIF

    ENDIF ! iaction


    IF (ltimer) CALL timer_stop(timer_physic_acc)

    IF (jstep > 1 .OR. (jstep == 1 .AND. lcall_phy_jg(itupdate))) THEN
     CALL nwp_diagnosis(lcall_phy_jg,lredgrid,jstep,     & !input

                            & tcall_phy_jg,p_sim_time,             & !input
                            & kstart_moist(jg),                    & !input
                            & pt_patch, pt_int_state, p_metrics,   & !input
                            & pt_prog, pt_prog_rcf,                & !in
                            & pt_diag,                            & !inout
                            & prm_diag,prm_nwp_tend)   
    END IF


    IF (ltimer) CALL timer_stop(timer_physics)

  END SUBROUTINE nwp_nh_interface

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------


  SUBROUTINE nh_update_prog_phy(pt_patch, pdtime, pt_diag, prm_diag, pt_prog_rcf)

    TYPE(t_patch),       INTENT(IN)   :: pt_patch     !!grid/patch info.
    TYPE(t_nh_diag),     INTENT(INOUT):: pt_diag      !!the diagnostic variables
    TYPE(t_nwp_phy_diag),INTENT(inout):: prm_diag     !!the physics variables
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
!$OMP DO PRIVATE(jb, jk,jc,jt,i_startidx, i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
           &                       i_startidx, i_endidx, rl_start, rl_end)

! KF fix to positive values
      DO jt=1,ntracer
        DO jk = kstart_moist(jg), nlev
          DO jc = i_startidx, i_endidx
            pt_prog_rcf%tracer(jc,jk,jb,jt) =MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt)  &
 &                                         + pdtime*pt_diag%ddt_tracer_phy(jc,jk,jb,jt))
          ENDDO
        ENDDO
      ENDDO

      prm_diag%rain_con(i_startidx:i_endidx,jb) =                                        &
  &                                          prm_diag%rain_con(i_startidx:i_endidx,jb)   &
  &                                        + pdtime                                          &
  &                                        * prm_diag%tracer_rate(i_startidx:i_endidx,jb,3)
      prm_diag%snow_con(i_startidx:i_endidx,jb) =                                        &
  &                                          prm_diag%snow_con(i_startidx:i_endidx,jb)   &
  &                                        + pdtime                                          &
  &                                        * prm_diag%tracer_rate(i_startidx:i_endidx,jb,4)

      prm_diag%tot_prec(i_startidx:i_endidx,jb) =                                        &
  &                                       prm_diag%tot_prec(i_startidx:i_endidx,jb)      &
  &                                    +  pdtime                                             &
  &                                    * (prm_diag%tracer_rate (i_startidx:i_endidx,jb,3)&
  &                                    +  prm_diag%tracer_rate (i_startidx:i_endidx,jb,4))


    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    IF (msg_level > 15) THEN


        CALL message('mo_nh_update_prog:', 'nh-tendencies')

        WRITE(message_text,'(a,2(E17.9,3x))') ' max/min ddt_qv  = ',&
          & MAXVAL(pt_diag%ddt_tracer_phy(:,:,:,iqv) ), &
          & MINVAL(pt_diag%ddt_tracer_phy(:,:,:,iqv) )
        CALL message('', TRIM(message_text))

        WRITE(message_text,'(a,2E17.9)') 'updated max/min QV  = ',&
          & MAXVAL(pt_prog_rcf%tracer(:,:,:,iqv)), MINVAL(pt_prog_rcf%tracer (:,:,: ,iqv) )
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,2E17.9)') 'updated max/min QC  = ',&
          & MAXVAL(pt_prog_rcf%tracer(:,:,:,iqc)), MINVAL(pt_prog_rcf%tracer (:,:,: ,iqc) )
        CALL message('', TRIM(message_text))
      ENDIF


  END SUBROUTINE nh_update_prog_phy

END MODULE mo_nh_interface_nwp

