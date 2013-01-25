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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_diagnosis


  USE mo_kind,               ONLY: wp

  USE mo_impl_constants,     ONLY: itccov, itfastphy, icc, min_rlcell_int
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_exception,          ONLY: message, message_text
  USE mo_model_domain,       ONLY: t_patch
  USE mo_run_config,         ONLY: msg_level, iqv, iqc
  USE mo_nonhydro_types,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,      ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_parallel_config,    ONLY: nproma
  USE mo_time_config,        ONLY: time_config
  USE mo_lnd_nwp_config,     ONLY: nlev_soil
  USE mo_physical_constants, ONLY: lh_v     => alv, &      !! latent heat of vapourization
                                   rd, cpd, rcvd
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
  USE mo_io_config,          ONLY: lflux_avg
  USE mo_sync,               ONLY: global_max, global_min

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC  :: nwp_diagnosis, nwp_diag_output_1, nwp_diag_output_2

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
  SUBROUTINE nwp_diagnosis(lcall_phy_jg,lredgrid,         & !input
                            & dt_phy_jg,p_sim_time,       & !input
                            & kstart_moist,               & !input
                            & pt_patch, p_metrics,        & !input
                            & pt_prog, pt_prog_rcf,       & !in
                            & pt_diag,                    & !inout
                            & prm_diag,prm_nwp_tend)    
                            

    !>
    ! !INPUT PARAMETERS:

    LOGICAL, INTENT(IN)          ::   &             !< physics package time control (switches)
         &                          lcall_phy_jg(:) !< for domain jg
    LOGICAL, INTENT(IN)          :: lredgrid        !< use reduced grid for radiation
    REAL(wp),INTENT(in)          :: dt_phy_jg(:)    !< time interval for all physics
                                                    !< packages on domain jg
    REAL(wp),INTENT(in)          :: p_sim_time

    TYPE(t_patch),        TARGET,INTENT(in):: pt_patch     !<grid/patch info.
    TYPE(t_nh_diag), TARGET, INTENT(inout) :: pt_diag     !<the diagnostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog     !<the prognostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog_rcf !<the prognostic variables (with
                                                          !< red. calling frequency for tracers!

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

    REAL(wp):: z_help, p_sim_time_s6, r_sim_time

    INTEGER :: jc,jk,jb,jg      !block index
    INTEGER :: kstart_moist

    INTEGER :: ioverlap(nproma)
    REAL(wp):: cld_cov(nproma)


    ! local variables related to the blocking

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id
    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1    

    ! Inverse of simulation time
    r_sim_time = 1._wp/MAX(1.e-6_wp, p_sim_time)

    ! exclude nest boundary interpolation zone
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)
    
! if cloud cover is call, vertical integration of cloud content(for iqv, iqc, iqi)

!$OMP PARALLEL PRIVATE(l_s6avg,p_sim_time_s6)

    IF ( lcall_phy_jg(itccov) ) THEN

!$OMP DO PRIVATE(jb, z_help,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
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

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk, ioverlap, cld_cov) ICON_OMP_DEFAULT_SCHEDULE
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
              cld_cov(jc) = 1._wp - (1._wp-cld_cov(jc)) * &
                         & (1._wp - prm_diag%tot_cld(jc,jk,jb,icc) )
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


    IF ( p_sim_time > 1.e-6_wp ) THEN

!$OMP DO PRIVATE(jb, i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx

         prm_diag%tot_cld_vi_avg(jc,jb,1:4) = ( prm_diag%tot_cld_vi_avg(jc,jb,1:4) &
                               &  * (p_sim_time - dt_phy_jg(itfastphy))       &
                               &  + prm_diag%tot_cld_vi(jc,jb,1:4)            &
                               &  * dt_phy_jg(itfastphy) )                    &
                               &  * r_sim_time 
        ENDDO
      ENDDO ! nblks     
!$OMP END DO

    END IF

 !! Calculate vertically integrated values of the grid-scale tracers q1, q2, q3, q4 and q5 
 !! and average values of the vertically integrated values from the model start 
 
!$OMP DO PRIVATE(jb, i_startidx,i_endidx,jc,jk,z_help) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      pt_diag%tracer_vi(i_startidx:i_endidx,jb,1:5) = 0.0_wp
      DO jk = kstart_moist, nlev
        DO jc = i_startidx, i_endidx 

          z_help = p_metrics%ddqz_z_full(jc,jk,jb) * pt_prog%rho(jc,jk,jb)  
          pt_diag%tracer_vi(jc, jb,1:5) = pt_diag%tracer_vi(jc, jb,1:5)    + &
                                 z_help * pt_prog_rcf%tracer(jc,jk,jb,1:5) 

        ENDDO
      ENDDO

      IF ( p_sim_time > 1.e-6_wp ) THEN
        DO jc = i_startidx, i_endidx 

          pt_diag%tracer_vi_avg(jc,jb,1:5) = ( pt_diag%tracer_vi_avg(jc,jb,1:5) &
                            &  * (p_sim_time - dt_phy_jg(itfastphy))      &
                            &  + pt_diag%tracer_vi(jc,jb,1:5)             &
                            &  * dt_phy_jg(itfastphy) )                   &
                            &  * r_sim_time
        ENDDO
      END IF
    ENDDO ! nblks   
!$OMP END DO


! average from the model start of total precipitation rate ,
!         convective precipitation rate and grid-scale precipitation rate


    IF ( p_sim_time > 1.e-6_wp ) THEN

!$OMP DO PRIVATE(jb, i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx

          prm_diag%tot_prec_rate_avg(jc,jb) =  prm_diag%tot_prec(jc,jb) &
                               & * r_sim_time
          prm_diag%con_prec_rate_avg(jc,jb) =  (prm_diag%rain_con(jc,jb) & 
                               & + prm_diag%snow_con(jc,jb))             &
                               & * r_sim_time
          prm_diag%gsp_prec_rate_avg(jc,jb) =  (prm_diag%rain_gsp(jc,jb) &
                               & + prm_diag%snow_gsp(jc,jb)) &
                               & * r_sim_time

        ENDDO
      ENDDO ! nblks     
!$OMP END DO
    END IF


    ! latent heat, latent heat from bare soil and sensible heat at surface. 
    ! Calculation of average/accumulated values since model start
    !
    IF ( p_sim_time > 1.e-6_wp .AND. lflux_avg) THEN

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        SELECT CASE (atm_phy_nwp_config(jg)%inwp_turb)
        CASE (1, 2, 3)
          DO jc = i_startidx, i_endidx
            prm_diag%alhfl_s(jc,jb) = ( prm_diag%alhfl_s(jc,jb)         &
                               &  * (p_sim_time - dt_phy_jg(itfastphy)) &
                               &  + prm_diag%lhfl_s(jc,jb)              &!attention to the sign, in the output all fluxes 
                               &  * dt_phy_jg(itfastphy) )              &!must be positive downwards 
                               & * r_sim_time
            prm_diag%alhfl_bs(jc,jb) = ( prm_diag%alhfl_bs(jc,jb)       &
                               &  * (p_sim_time - dt_phy_jg(itfastphy)) &
                               &  + prm_diag%lhfl_bs(jc,jb)             &!attention to the sign, in the output all fluxes
                               &  * dt_phy_jg(itfastphy) )              &!must be positive downwards 
                               & * r_sim_time
            prm_diag%ashfl_s(jc,jb) = ( prm_diag%ashfl_s(jc,jb)         &
                               &  * (p_sim_time - dt_phy_jg(itfastphy)) &
                               &  + prm_diag%shfl_s(jc,jb)              &!attention to the sign, in the output all fluxes
                               &  * dt_phy_jg(itfastphy) )              &!must be positive downwards 
                               & * r_sim_time
          ENDDO  ! jc
          DO jk = 1, nlev_soil
            DO jc = i_startidx, i_endidx
            prm_diag%alhfl_pl(jc,jk,jb) = ( prm_diag%alhfl_pl(jc,jk,jb)       &
                               &  * (p_sim_time - dt_phy_jg(itfastphy)) &
                               &  + prm_diag%lhfl_pl(jc,jk,jb)          &!attention to the sign, in the output all fluxes
                               &  * dt_phy_jg(itfastphy) )              &!must be positive downwards 
                               & * r_sim_time
            ENDDO  ! jc
          ENDDO  ! jk
        CASE (4)
          DO jc = i_startidx, i_endidx
            prm_diag%alhfl_s(jc,jb) = ( prm_diag%alhfl_s(jc,jb)         &
                               &  * (p_sim_time - dt_phy_jg(itfastphy)) &
                               &  + prm_diag%qhfl_s(jc,jb)*lh_v         &
                               &  * dt_phy_jg(itfastphy) )              &
                               & * r_sim_time
            prm_diag%ashfl_s(jc,jb) = ( prm_diag%ashfl_s(jc,jb)         &
                               &  * (p_sim_time - dt_phy_jg(itfastphy)) &
                               &  + prm_diag%shfl_s(jc,jb)              &! it is 0 at the moment, with turb2 the
                               &  * dt_phy_jg(itfastphy) )              &! sensible heat is not output
                               & * r_sim_time
          ENDDO
        END SELECT

      ENDDO ! nblks     
!$OMP END DO

    ELSEIF (.NOT. lflux_avg) THEN
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        SELECT CASE (atm_phy_nwp_config(jg)%inwp_turb)
        CASE (1, 2, 3)
          DO jc = i_startidx, i_endidx
            prm_diag%alhfl_s(jc,jb) =  prm_diag%alhfl_s(jc,jb)       &
                               &  + prm_diag%lhfl_s(jc,jb)           &!attention to the sign, in the output all fluxes 
                               &  * dt_phy_jg(itfastphy)              !must be positive downwards 
            prm_diag%alhfl_bs(jc,jb) =  prm_diag%alhfl_bs(jc,jb)     &
                               &  + prm_diag%lhfl_bs(jc,jb)          &!attention to the sign, in the output all fluxes 
                               &  * dt_phy_jg(itfastphy)              !must be positive downwards 
            prm_diag%ashfl_s(jc,jb) =  prm_diag%ashfl_s(jc,jb)       &
                               &  + prm_diag%shfl_s(jc,jb)           &!attention to the sign, in the output all fluxes 
                               &  * dt_phy_jg(itfastphy)              !must be positive downwards 
          ENDDO
          DO jk = 1, nlev_soil
            DO jc = i_startidx, i_endidx
              prm_diag%alhfl_pl(jc,jk,jb) =  prm_diag%alhfl_pl(jc,jk,jb)&
                               &  + prm_diag%lhfl_pl(jc,jk,jb)          &!attention to the sign, in the output all fluxes 
                               &  * dt_phy_jg(itfastphy)                 !must be positive downwards 
            ENDDO  ! jc
          ENDDO  ! jk
        CASE (4)
          DO jc = i_startidx, i_endidx
            prm_diag%alhfl_s(jc,jb) =  prm_diag%alhfl_s(jc,jb)       &
                               &  + prm_diag%qhfl_s(jc,jb)*lh_v      &
                               &  * dt_phy_jg(itfastphy)
            prm_diag%ashfl_s(jc,jb) =  prm_diag%ashfl_s(jc,jb)       &
                               &  + prm_diag%shfl_s(jc,jb)           &! it is 0 at the moment, with turb2 the
                               &  * dt_phy_jg(itfastphy)              ! sensible heat is not output
          ENDDO
        END SELECT

      ENDDO ! nblks     
!$OMP END DO    

    END IF

! Evaporation rate at surface.  Calculation of 
! average values since model start


    IF ( p_sim_time > 1.e-6_wp ) THEN
!$OMP DO PRIVATE(jb, i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        SELECT CASE (atm_phy_nwp_config(jg)%inwp_turb)
        CASE (1, 2, 3)
          DO jc = i_startidx, i_endidx
            prm_diag%aqhfl_s(jc,jb) = ( prm_diag%aqhfl_s(jc,jb)         &
                               &  * (p_sim_time - dt_phy_jg(itfastphy)) &
                               &  + prm_diag%lhfl_s(jc,jb)/lh_v         & !attention to the sign, in the output all fluxes  
                               &  * dt_phy_jg(itfastphy) )              & !must be positive downwards 
                               & * r_sim_time
          ENDDO
        CASE (4)
          DO jc = i_startidx, i_endidx
             prm_diag%aqhfl_s(jc,jb) = ( prm_diag%aqhfl_s(jc,jb)        &
                               &  * (p_sim_time - dt_phy_jg(itfastphy)) &
                               &  + prm_diag%qhfl_s(jc,jb)              &
                               &  * dt_phy_jg(itfastphy) )              &
                               & * r_sim_time
          ENDDO
        END SELECT

      ENDDO ! nblks     
!$OMP END DO
    END IF


 
! Check if it is 00, 06, 12 or 18 UTC. In this case update the value of 
!    dt_s6avg average variables 

    IF (MOD(p_sim_time + time_config%ini_datetime%daysec, dt_s6avg) == 0._wp &
      & .AND. p_sim_time > 0.01_wp) THEN
      l_s6avg = .TRUE.
      p_sim_time_s6 = REAL(INT( (p_sim_time+time_config%ini_datetime%daysec)/dt_s6avg),wp) &
         &             * dt_s6avg
    ELSE
      l_s6avg = .FALSE.
    END IF
   ! WRITE(0,*) "diag", p_sim_time, time_config%ini_datetime%daysec, dt_s6avg, & 
   !               & MOD(p_sim_time + time_config%ini_datetime%daysec, dt_s6avg), l_s6avg

   
    IF (l_s6avg ) THEN
!$OMP DO PRIVATE(jb, i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
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
!$OMP END DO NOWAIT
    END IF

!$OMP END PARALLEL  


  END SUBROUTINE nwp_diagnosis

  !-------------------------------------------------------------------------
  !>
  !! Extended diagnostics for NWP physics interface - part 1
  !! Was included in mo_nh_interface_nwp before
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2013-01-07)
  !!
  SUBROUTINE nwp_diag_output_1(p_patch, p_diag, p_prog_rcf)

    TYPE(t_patch),   INTENT(in) :: p_patch     !< grid/patch info.
    TYPE(t_nh_diag), INTENT(in) :: p_diag      !< NH diagnostic state
    TYPE(t_nh_prog), INTENT(in) :: p_prog_rcf  !< state for tracer variables


    ! Local variables
    REAL(wp), DIMENSION(p_patch%nblks_c,p_patch%nlev) ::             &
      maxabs_u, maxabs_v, maxtemp, mintemp, maxqv, minqv, maxqc, minqc
    REAL(wp), DIMENSION(p_patch%nlev) ::               &
      umax, vmax, tmax, tmin, qvmax, qvmin, qcmax, qcmin

    ! loop indices
    INTEGER :: jc,jk,jb,jg

    INTEGER :: nlev                    !< number of full levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    CALL message('mo_nwp_diagnosis:','Initial diagnostic output')

    maxabs_u(:,:) = 0._wp
    maxabs_v(:,:) = 0._wp
    maxtemp(:,:)  = 0._wp
    mintemp(:,:)  = 1.e20_wp
    maxqv(:,:)    = 0._wp
    minqv(:,:)    = 1.e20_wp
    maxqc(:,:)    = 0._wp
    minqc(:,:)    = 1.e20_wp

    nlev = p_patch%nlev
    jg   = p_patch%id

    ! Exclude the nest boundary zone
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_nchdom  = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          maxabs_u(jb,jk) = MAX(maxabs_u(jb,jk),ABS(p_diag%u(jc,jk,jb)))
          maxabs_v(jb,jk) = MAX(maxabs_v(jb,jk),ABS(p_diag%v(jc,jk,jb)))
          maxtemp(jb,jk)  = MAX(maxtemp(jb,jk),p_diag%temp(jc,jk,jb))
          mintemp(jb,jk)  = MIN(mintemp(jb,jk),p_diag%temp(jc,jk,jb))
          maxqv(jb,jk)    = MAX(maxqv(jb,jk),p_prog_rcf%tracer(jc,jk,jb,iqv))
          minqv(jb,jk)    = MIN(minqv(jb,jk),p_prog_rcf%tracer(jc,jk,jb,iqv))
          maxqc(jb,jk)    = MAX(maxqc(jb,jk),p_prog_rcf%tracer(jc,jk,jb,iqc))
          minqc(jb,jk)    = MIN(minqc(jb,jk),p_prog_rcf%tracer(jc,jk,jb,iqc))
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

  END SUBROUTINE nwp_diag_output_1


  !-------------------------------------------------------------------------
  !>
  !! Extended diagnostics for NWP physics interface - part 2
  !! Was included in mo_nh_interface_nwp before
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2013-01-07)
  !!
  SUBROUTINE nwp_diag_output_2(p_patch, p_diag, p_prog_rcf, prm_nwp_tend, dtime, lcall_turb)

    TYPE(t_patch), TARGET,INTENT(in) :: p_patch      !< grid/patch info.
    TYPE(t_nh_diag),      INTENT(in) :: p_diag       !< NH diagnostic state
    TYPE(t_nh_prog),      INTENT(in) :: p_prog_rcf   !< state for TKE
    TYPE(t_nwp_phy_tend), INTENT(in) :: prm_nwp_tend !< physics tendencies

    REAL(wp), INTENT(in) :: dtime      ! time step
    LOGICAL,  INTENT(in) :: lcall_turb ! switch if turbulence has been called

    ! Local variables

    ! variables for turbulence diagnostics
    REAL(wp) :: maxtke(p_patch%nblks_c,p_patch%nlevp1),tkemax(p_patch%nlevp1)
    REAL(wp), DIMENSION(p_patch%nblks_c,p_patch%nlev) :: maxtturb, maxuturb, maxvturb
    REAL(wp), DIMENSION(p_patch%nlev) :: tturbmax, uturbmax, vturbmax

    ! variables for CFL diagnostic
    REAL(wp) :: maxcfl(p_patch%nblks_c), cflmax, avg_invedgelen(nproma), csfac

    INTEGER,  POINTER :: ieidx(:,:,:), ieblk(:,:,:)

    ! loop indices
    INTEGER :: jc,jk,jb,jg

    INTEGER :: nlev, nlevp1            !< number of model levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1
    jg     = p_patch%id

    ! Exclude the nest boundary zone
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_nchdom  = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    ! factor for sound speed computation
    csfac = rd*cpd*rcvd

    ! Initialization
    maxcfl(:) = 0._wp

    ! In case that turbulence diagnostics are computed
    IF (msg_level >= 18) THEN
      maxtke(:,:)   = 0._wp
      maxtturb(:,:) = 0._wp
      maxuturb(:,:) = 0._wp
      maxvturb(:,:) = 0._wp
    ENDIF

!$OMP PARALLEL

    ! CFL-diagnostic

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,avg_invedgelen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)


      DO jc = i_startidx, i_endidx
        avg_invedgelen(jc) = 3._wp/                                       &
          (p_patch%edges%dual_edge_length(ieidx(jc,jb,1),ieblk(jc,jb,1))+ &
           p_patch%edges%dual_edge_length(ieidx(jc,jb,2),ieblk(jc,jb,2))+ &
           p_patch%edges%dual_edge_length(ieidx(jc,jb,3),ieblk(jc,jb,3))  )
      ENDDO

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          maxcfl(jb) = MAX(maxcfl(jb),dtime*avg_invedgelen(jc)*( &
            SQRT(p_diag%u(jc,jk,jb)**2+p_diag%v(jc,jk,jb)**2)+   &
            SQRT(csfac*p_diag%temp(jc,jk,jb)) ))
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO

    ! Extended turbulence diagnostics if msg_level >= 18
    IF (lcall_turb .AND. msg_level >= 18) THEN

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO jk = 1, nlevp1
          DO jc = i_startidx, i_endidx
            maxtke(jb,jk) = MAX(maxtke(jb,jk),p_prog_rcf%tke(jc,jk,jb))
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

    ! CFL diagnostic
    cflmax = MAXVAL(maxcfl)
    cflmax = global_max(cflmax) ! maximum over all PEs
    WRITE(message_text,'(a,f12.8,a,i2)') 'maximum horizontal CFL = ', cflmax, ' in domain ',jg
    CALL message('nwp_nh_interface: ', TRIM(message_text))


    IF (msg_level >= 18 .AND. lcall_turb) THEN ! extended turbulence diagnostic
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

  END SUBROUTINE nwp_diag_output_2

END MODULE mo_nwp_diagnosis

