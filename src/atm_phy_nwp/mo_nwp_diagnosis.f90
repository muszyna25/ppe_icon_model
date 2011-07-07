!>
!! @brief diagnosis of physics after physic's call 
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! first implementation  by <Pilar Ripodas, DWD> (2011-03)
!!
!! $Id: n/a$
!!
!! @par Copyright
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
MODULE mo_nwp_diagnosis


  USE mo_kind,               ONLY: wp

  USE mo_timer,              ONLY: timer_physics, timer_start, timer_stop
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_impl_constants,     ONLY: itconv, itccov, itrad, itgscp,         &
    &                              itsatad, itupdate, itturb, itradheat,  &
    &                              itsso, icc,                            &
    &                              min_rlcell_int, min_rledge_int, min_rlcell
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_interpolation,      ONLY: t_int_state, rbf_vec_interpol_cell,    &
                                   edges2cells_scalar
  USE mo_grf_interpolation,  ONLY: t_gridref_single_state, &
    &                              t_gridref_state

  USE mo_model_domain,       ONLY: t_patch
  
  USE mo_nonhydro_state,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_state,      ONLY: t_nwp_phy_diag,&
                                 & t_nwp_phy_tend,  prm_diag
  USE mo_parallel_configuration,  ONLY: nproma
  USE mo_run_nml,            ONLY: ntracer, iqv, iqc, iqi, &
       &                              iqr, iqs, msg_level, ltimer
  USE mo_time_config,        ONLY: time_config
  USE mo_physical_constants, ONLY: rd, rd_o_cpd, vtmpc1, p0ref, cvd_o_rd, &
                                   lh_v     => alv      !! latent heat of vapourization

  USE mo_atm_phy_nwp_nml,    ONLY: inwp_cldcover, inwp_radiation,  dt_rad,&
                                   inwp_sso, inwp_turb 
  USE mo_sync,               ONLY: sync_patch_array, sync_patch_array_mult, &
                                   SYNC_C, SYNC_C1
  USE mo_mpi,                ONLY: p_nprocs
  USE mo_parallel_configuration,  ONLY: p_test_run
 


  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC  :: nwp_diagnosis

CONTAINS

  !>
  !! <Short description of the subroutine for listings and indices>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!
  SUBROUTINE nwp_diagnosis(lcall_phy_jg,lredgrid,jstep,     & !input
                            & tcall_phy_jg,p_sim_time,             & !input
                            & kstart_moist,                        & !input
                            & pt_patch, pt_int_state, p_metrics,   & !input
                            & pt_prog, pt_prog_rcf,                & !in
                            & pt_diag,                            & !inout
                            & prm_diag,prm_nwp_tend)    
                            

    !>
    ! !INPUT PARAMETERS:

    LOGICAL, INTENT(IN)          ::   &             !< physics package time control (switches)
         &                          lcall_phy_jg(:) !< for domain jg
    LOGICAL, INTENT(IN)          :: lredgrid        !< use reduced grid for radiation
    INTEGER ,INTENT(in)          :: jstep
    REAL(wp),INTENT(in)          :: tcall_phy_jg(:) !< time interval for all physics
                                                    !< packages on domain jg
    REAL(wp),INTENT(in)          :: p_sim_time

    TYPE(t_patch),        TARGET,INTENT(in):: pt_patch     !<grid/patch info.
    TYPE(t_nh_diag), TARGET, INTENT(inout) :: pt_diag     !<the diagnostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog     !<the prognostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog_rcf !<the prognostic variables (with
                                                          !< red. calling frequency for tracers!


    TYPE(t_int_state),    TARGET,INTENT(in):: pt_int_state      !< interpolation state
    TYPE(t_nh_metrics)   ,       INTENT(in):: p_metrics

    TYPE(t_nwp_phy_diag),       INTENT(inout) :: prm_diag
    TYPE(t_nwp_phy_tend),TARGET,INTENT(inout) :: prm_nwp_tend

    ! Local

    REAL(wp), PARAMETER :: dt_s6avg = 21600._wp   !6 hours in seconds
    LOGICAL :: l_s6avg

    ! Local array bounds:

    INTEGER :: nlev, nlevp1            !< number of full levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    REAL(wp):: z_help, p_sim_time_s6

    INTEGER :: jc,jk,jb,jt      !block index
    INTEGER :: kstart_moist

    INTEGER :: ioverlap(nproma)
    REAL(wp):: cld_cov(nproma)


    ! local variables related to the blocking

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1    


   !in order to account for mesh refinement
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)
    
! if cloud cover is call, vertical integration of cloud content(for iqv, iqc, iqi)

!$OMP PARALLEL

    IF ( lcall_phy_jg(itccov) ) THEN

!$OMP DO PRIVATE(jb, z_help,i_startidx,i_endidx,jc,jk)
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        prm_diag%tot_cld_vi(i_startidx:i_endidx,jb,1:3) = 0.0_wp

        DO jk = kstart_moist, nlev
          DO jc = i_startidx, i_endidx

           z_help = p_metrics%ddqz_z_full(jc,jk,jb) * pt_prog%rho(jc,jk,jb)  
           prm_diag%tot_cld_vi(jc, jb,1:3) = prm_diag%tot_cld_vi(jc, jb,1:3)    + &
                                             z_help * prm_diag%tot_cld(jc,jk,jb,1:3) 
          ENDDO
        ENDDO
      ENDDO ! nblks  
!$OMP END DO

! Calculation of cloud cover (Maximum-Random Overlap) (icc)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk, ioverlap, cld_cov)
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)
          ioverlap(:) = 1
          cld_cov(:)  = 0.0_wp
        DO jk = kstart_moist, nlev
          DO jc = i_startidx, i_endidx
           IF (ioverlap(jc) == 1) THEN
            cld_cov(jc) = MAX(prm_diag%tot_cld(jc,jk,jb,icc), cld_cov(jc))
           ELSE
            cld_cov(jc) = 1._wp - (1._wp-cld_cov(jc)) * (1._wp - prm_diag%tot_cld(jc,jk,jb,icc))
           END IF
           IF (prm_diag%tot_cld(jc,jk,jb,icc) <= 1.e-6_wp) THEN
             ioverlap(jc)=0
           ELSE
             ioverlap(jc)=1
           ENDIF
          ENDDO
        ENDDO
        prm_diag%tot_cld_vi(:,jb,icc) = cld_cov(:)
      ENDDO ! nblks  
!$OMP END DO

    END IF !cloud cover

! average values of the vertically integrated total cloud contents (for iqv, iqc, iqi, icc)
! from the model start


    IF ( p_sim_time .GT. 1.e-1 ) THEN

!$OMP DO PRIVATE(jb, i_startidx,i_endidx,jc)
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)
          DO jc = i_startidx, i_endidx

           prm_diag%tot_cld_vi_avg(jc,jb,1:4) = ( prm_diag%tot_cld_vi_avg(jc,jb,1:4) &
                               &  * (p_sim_time - tcall_phy_jg(itupdate))        &
                               &  + prm_diag%tot_cld_vi(jc,jb,1:4)               &
                               &  * tcall_phy_jg(itupdate) )                 &
                               & / p_sim_time 
          ENDDO
      ENDDO ! nblks     
!$OMP END DO

    END IF

 !! Calculate vertically integrated values of the grid-scale tracers q1, q2 and q3 
 !! and average values of the vertically integrated values from the model start 
 
!$OMP DO PRIVATE(jb, i_startidx,i_endidx,jc,jk,z_help)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        pt_diag%tracer_vi(i_startidx:i_endidx,jb,1:3) = 0.0_wp
        DO jk = kstart_moist, nlev
          DO jc = i_startidx, i_endidx 

            z_help = p_metrics%ddqz_z_full(jc,jk,jb) * pt_prog%rho(jc,jk,jb)  
            pt_diag%tracer_vi(jc, jb,1:3) = pt_diag%tracer_vi(jc, jb,1:3)    + &
                                            z_help * pt_prog_rcf%tracer(jc,jk,jb,1:3) 

          ENDDO
        ENDDO
        IF ( p_sim_time .GT. 1.e-1 ) THEN
         DO jc = i_startidx, i_endidx 
          pt_diag%tracer_vi_avg(jc,jb,1:3) = ( pt_diag%tracer_vi_avg(jc,jb,1:3) &
                              &  * (p_sim_time - tcall_phy_jg(itupdate))    &
                              &  + pt_diag%tracer_vi(jc,jb,1:3)             &
                              &  * tcall_phy_jg(itupdate) )                 &
                              & / p_sim_time
         ENDDO
        END IF
      ENDDO ! nblks   
!$OMP END DO


! average from the model start of total precipitation rate ,
!         convective precipitation rate and grid-scale precipitation rate


    IF ( p_sim_time .GT. 1.e-1 ) THEN

!$OMP DO PRIVATE(jb, i_startidx,i_endidx,jc)
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)
          DO jc = i_startidx, i_endidx

           prm_diag%tot_prec_rate_avg(jc,jb) =  prm_diag%tot_prec(jc,jb) &
                               & / p_sim_time  
           prm_diag%con_prec_rate_avg(jc,jb) =  (prm_diag%rain_con(jc,jb) & 
                               & + prm_diag%snow_con(jc,jb))             &
                               & / p_sim_time 
           prm_diag%gsp_prec_rate_avg(jc,jb) =  (prm_diag%rain_gsp(jc,jb) &
                               & + prm_diag%snow_gsp(jc,jb)) &
                               & / p_sim_time

          ENDDO
      ENDDO ! nblks     
!$OMP END DO
     END IF


! latent heat, sensible heat and evaporation rate at surface. Calculation of average values 
! since model start

    IF ( p_sim_time .GT. 1.e-1 ) THEN

!$OMP DO PRIVATE(jb, i_startidx,i_endidx,jc)
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)
          DO jc = i_startidx, i_endidx

           IF (inwp_turb == 1) THEN
            prm_diag%lhfl_s_avg(jc,jb) = ( prm_diag%lhfl_s_avg(jc,jb)     &
                               &  * (p_sim_time - tcall_phy_jg(itupdate)) &
                               &  + prm_diag%lhfl_s(jc,jb)                &
                               &  * tcall_phy_jg(itupdate) )              &
                               & / p_sim_time 
           ELSEIF (inwp_turb == 2) THEN
            prm_diag%lhfl_s_avg(jc,jb) = ( prm_diag%lhfl_s_avg(jc,jb)     &
                               &  * (p_sim_time - tcall_phy_jg(itupdate)) &
                               &  + prm_diag%qhfl_s(jc,jb)*lh_v           &
                               &  * tcall_phy_jg(itupdate) )              &
                               & / p_sim_time 

           ENDIF
           prm_diag%shfl_s_avg(jc,jb) = ( prm_diag%shfl_s_avg(jc,jb)      &
                               &  * (p_sim_time - tcall_phy_jg(itupdate)) &
                               &  + prm_diag%shfl_s(jc,jb)                &
                               &  * tcall_phy_jg(itupdate) )              &
                               & / p_sim_time 
           IF (inwp_turb == 1) THEN
             prm_diag%qhfl_s_avg(jc,jb) = ( prm_diag%qhfl_s_avg(jc,jb)    &
                               &  * (p_sim_time - tcall_phy_jg(itupdate)) &
                               &  + prm_diag%lhfl_s(jc,jb)/lh_v           &
                               &  * tcall_phy_jg(itupdate) )              &
                               & / p_sim_time

           ELSEIF (inwp_turb == 2) THEN
             prm_diag%qhfl_s_avg(jc,jb) = ( prm_diag%qhfl_s_avg(jc,jb)    &
                               &  * (p_sim_time - tcall_phy_jg(itupdate)) &
                               &  + prm_diag%qhfl_s(jc,jb)                &
                               &  * tcall_phy_jg(itupdate) )              &
                               & / p_sim_time
           ENDIF

          ENDDO
      ENDDO ! nblks     
!$OMP END DO

    END IF
 
! Check if it is 00, 06, 12 or 18 UTC. In this case update the value of 
!    dt_s6avg average variables 

    IF (MOD(p_sim_time + time_config%ini_datetime%daysec, dt_s6avg) == 0 &
      & .AND. p_sim_time > 0.1) THEN
       l_s6avg = .TRUE.
       p_sim_time_s6 = INT( (p_sim_time+time_config%ini_datetime%daysec)/dt_s6avg) * dt_s6avg
    ELSE
       l_s6avg = .FALSE.
    END IF
   ! WRITE(0,*) "diag", p_sim_time, time_config%ini_datetime%daysec, dt_s6avg, & 
   !               & MOD(p_sim_time + time_config%ini_datetime%daysec, dt_s6avg), l_s6avg

   
    IF (l_s6avg ) THEN
!$OMP DO PRIVATE(jb, i_startidx,i_endidx,jc)
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)
          DO jc = i_startidx, i_endidx
            prm_diag%u_10m_s6avg(jc,jb) = ( prm_diag%u_10m_s6avg(jc,jb)   &
                                    & * (p_sim_time_s6 - dt_s6avg)        &
                                    & + prm_diag%u_10m(jc,jb)             &
                                    & * dt_s6avg )                        &
                                    & / p_sim_time_s6
            prm_diag%v_10m_s6avg(jc,jb) = ( prm_diag%v_10m_s6avg(jc,jb)   &
                                    & * (p_sim_time_s6 - dt_s6avg)        &
                                    & + prm_diag%v_10m(jc,jb)             &
                                    & * dt_s6avg )                        &
                                    & / p_sim_time_s6
            prm_diag%t_2m_s6avg(jc,jb) = ( prm_diag%t_2m_s6avg(jc,jb)     &
                                    & * (p_sim_time_s6 - dt_s6avg)        &
                                    & + prm_diag%t_2m(jc,jb)              &
                                    & * dt_s6avg )                        &
                                    & / p_sim_time_s6
            prm_diag%qv_2m_s6avg(jc,jb) = ( prm_diag%qv_2m_s6avg(jc,jb)   &
                                    & * (p_sim_time_s6 - dt_s6avg)        &
                                    & + prm_diag%qv_2m(jc,jb)             &
                                    & * dt_s6avg )                        &
                                    & / p_sim_time_s6


            pt_diag%pres_sfc_s6avg(jc,jb) = ( pt_diag%pres_sfc_s6avg(jc,jb) &
                                    & * (p_sim_time_s6 - dt_s6avg)          &
                                    & + pt_diag%pres_sfc(jc,jb)             &
                                    & * dt_s6avg )                          &
                                    & / p_sim_time_s6
          ENDDO
      ENDDO ! nblks     
!$OMP END DO
    END IF

!$OMP END PARALLEL  


  END SUBROUTINE nwp_diagnosis


END MODULE mo_nwp_diagnosis

