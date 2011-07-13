!OPTION! -cont
!! this command should fix the problem of copying arrays in a subroutine call
!>
!! This module is the interface between nwp_nh_interface to the 
!! microphysical parameterisations:
!! inwp_gscp == 1 == one_moment bulk microphysics by Doms and Schaettler(2004) 
!!                                                and Seifert and Beheng(2006)
!! inwp_gscp == 2 == ...
!!
!! @author Kristina Froehlich, DWD, Offenbach (2010-01-25)
!!
!! @par Revision History
!! Initial Kristina Froehlich, DWD, Offenbach (2010-01-25)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_nwp_gscp_interface

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, message_text, finish
  USE mo_mpi,                  ONLY: p_pe, p_nprocs
  USE mo_parallel_configuration,  ONLY: p_test_run, nproma
  USE mo_sync,                 ONLY: SYNC_C, sync_patch_array_mult

  USE mo_model_domain,         ONLY: t_patch
  USE mo_grf_interpolation,    ONLY: t_gridref_state
  USE mo_impl_constants,       ONLY: min_rlcell_int
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c

  USE mo_nonhydro_state,       ONLY: t_nh_prog, t_nh_diag,&
    &                                t_nh_metrics
  USE mo_nonhydrostatic_nml,   ONLY: kstart_moist
  USE mo_nwp_phy_state,        ONLY: t_nwp_phy_diag,prm_diag,&
    &                                t_nwp_phy_tend
  USE mo_run_config,           ONLY: msg_level, ntracer, iqv, iqc, &
    &                                iqi, iqr, iqs
!  USE mo_atm_phy_nwp_nml,      ONLY: inwp_gscp
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_gscp_cosmo,           ONLY: hydci_pp
  USE mo_satad,                ONLY: satad_v_3D

  IMPLICIT NONE

  PRIVATE

  PUBLIC  ::  nwp_microphysics

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE nwp_microphysics( tcall_gscp_jg,                & !>input
                            &   p_patch,p_metrics,            & !>input
                            &   p_prog,                       & !>inout
                            &   p_prog_rcf,                   & !>inout
                            &   p_diag ,                      & !>inout
                            &   prm_diag                      ) !>inout 



    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !!<grid/patch info.
    TYPE(t_nh_metrics)          ,INTENT(in)   :: p_metrics
    TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog          !<the dyn prog vars
    TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog_rcf      !<call freq
    TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag          !<the dyn diag vars
    TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag        !<the atm phys vars

    REAL(wp),                    INTENT(in)   :: tcall_gscp_jg   !< time interval for 
                                                                 !< microphysics
    ! Local array bounds:

    INTEGER :: nblks_c, nblks_e        !< number of blocks for cells / edges
    INTEGER :: nlev                    !< number of full levels
    INTEGER :: npromz_e, npromz_c      !< length of last block line
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !< blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    ! Local scalars:

    INTEGER :: jc,jb,jg               !<block indeces

 ! local variables related to the blocking

    nblks_c   = p_patch%nblks_int_c
    npromz_c  = p_patch%npromz_int_c
    nblks_e   = p_patch%nblks_int_e
    npromz_e  = p_patch%npromz_int_e
    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! number of vertical levels
    nlev = p_patch%nlev

    ! domain ID
    jg = p_patch%id

    !in order to account for mesh refinement
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx), SCHEDULE(guided)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)

          !>  prognostic microphysics and precipitation scheme from COSMO
          !!  NOTE: since microphysics is a fast process it is
          !!        allowed to do a sequential updating!

        IF( atm_phy_nwp_config(jg)%inwp_gscp == 1 ) THEN


!#ifdef __BOUNDCHECK
          CALL hydci_pp (                                           &
            & ie     =nproma                            ,    & !> in:  actual array size
            & je     =1                                 ,    & !< in:  dummy  array dimension
            & ke     =nlev                              ,    & !< in:  actual array size
            & istart =i_startidx                        ,    & !< in:  start index of calculation
            & iend   =i_endidx                          ,    & !< in:  end index of calculation
            & kstart =kstart_moist(jg)                  ,    & !< in:  vertical start index
            & zdt    =tcall_gscp_jg                     ,    & !< in:  timestep
            & qi0    = atm_phy_nwp_config(jg)%qi0       ,    & 
            & qc0    = atm_phy_nwp_config(jg)%qc0       ,    & 
            & dz     =p_metrics%ddqz_z_full(:,:,jb)     ,    & !< in:  vertical layer thickness
            & t      =p_diag%temp   (:,:,jb)           ,    & !< in:  temp,tracer,...
            & p      =p_diag%pres   (:,:,jb)           ,    & !< in:  full level pres
            & rho    =p_prog%rho    (:,:,jb  )         ,    & !< in:  density
            & qv     =p_prog_rcf%tracer (:,:,jb,iqv)   ,    & !< in:  spec. humidity
            & qc     =p_prog_rcf%tracer (:,:,jb,iqc)   ,    & !< in:  cloud water
            & qi     =p_prog_rcf%tracer (:,:,jb,iqi)   ,    & !< in:  cloud ice
            & qr     =p_prog_rcf%tracer (:,:,jb,iqr)   ,    & !< in:  rain water
            & qs     =p_prog_rcf%tracer (:,:,jb,iqs)   ,    & !< in:  snow
            & prr_gsp=prm_diag%tracer_rate (:,jb,1)     ,    & !< out: precipitation rate of rain
            & prs_gsp=prm_diag%tracer_rate (:,jb,2)     ,    & !< out: precipitation rate of snow
            & idbg=msg_level, l_cv=.TRUE. )
!#else
!          CALL hydci_pp (                                           &
!            & ie     =nproma                            ,    & !> in:  actual array size
!            & je     =1                                 ,    & !< in:  dummy  array dimension
!            & ke     =nlev                              ,    & !< in:  actual array size
!            & istart =i_startidx                        ,    & !< in:  start index of calculation
!            & iend   =i_endidx                          ,    & !< in:  end index of calculation
!            & kstart =kstart_moist(jg)                  ,    & !< in:  vertical start index
!            & zdt    =tcall_gscp_jg                     ,    & !< in:  timestep
!            & dz     =p_metrics%ddqz_z_full(1,1,jb)     ,    & !< in:  vertical layer thickness
!            & t      =p_diag%temp   (1,1,jb)           ,    & !< in:  temp,tracer,...
!            & p      =p_diag%pres   (1,1,jb)           ,    & !< in:  full level pres
!            & rho    =p_prog%rho    (1,1,jb  )         ,    & !< in:  density
!            & qv     =p_prog_rcf%tracer (1,1,jb,iqv)   ,    & !< in:  spec. humidity
!            & qc     =p_prog_rcf%tracer (1,1,jb,iqc)   ,    & !< in:  cloud water
!            & qi     =p_prog_rcf%tracer (1,1,jb,iqi)   ,    & !< in:  cloud ice
!            & qr     =p_prog_rcf%tracer (1,1,jb,iqr)   ,    & !< in:  rain water
!            & qs     =p_prog_rcf%tracer (1,1,jb,iqs)   ,    & !< in:  snow
!            & prr_gsp=prm_diag%tracer_rate (1,jb,1)     ,    & !< out: precipitation rate of rain
!            & prs_gsp=prm_diag%tracer_rate (1,jb,2)     ,    & !< out: precipitation rate of snow
!            & idbg=msg_level, l_cv=.TRUE. )
!!#endif

        ENDIF ! inwp_gscp

        !-------------------------------------------------------------------------
        !>
        !! Calculate surface precipitation
        !!
        !-------------------------------------------------------------------------

        DO jc =  i_startidx, i_endidx

          prm_diag%rain_gsp(jc,jb) = prm_diag%rain_gsp(jc,jb)        &
 &                                      + tcall_gscp_jg              &
 &                                      * prm_diag%tracer_rate (jc,jb,1)

          prm_diag%snow_gsp(jc,jb) = prm_diag%snow_gsp(jc,jb)        &
 &                                      +tcall_gscp_jg               &
 &                                      * prm_diag%tracer_rate (jc,jb,2)

          prm_diag%tot_prec(jc,jb) = prm_diag%tot_prec(jc,jb)            &
 &                                      +  tcall_gscp_jg                 &
 &                                      * ( prm_diag%tracer_rate (jc,jb,1)   &
 &                                      +   prm_diag%tracer_rate (jc,jb,2) )
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL


      IF (msg_level >= 15) THEN
        WRITE(message_text,'(a,2E15.5)') 'max/min QV  = ',&
          & MAXVAL(p_prog_rcf%tracer(:,:,:,iqv)), MINVAL(p_prog_rcf%tracer (:,:,: ,iqv) )
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,2E15.5)') ' max/min QC  = ',&
          & MAXVAL(p_prog_rcf%tracer(:,:,:,iqc)), MINVAL(p_prog_rcf%tracer (:,:,: ,iqc) )
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,2E15.5)') ' max/min QI  = ',&
          & MAXVAL(p_prog_rcf%tracer(:,:,:,iqi)), MINVAL(p_prog_rcf%tracer (:,:,: ,iqi) )
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,2E15.5)') ' max/min T  = ',&
             & MAXVAL(p_diag%temp(:,:,:)), MINVAL(p_diag%temp (:,:,:) )
        CALL message('', TRIM(message_text))
      ENDIF
       
  END SUBROUTINE nwp_microphysics

END MODULE mo_nwp_gscp_interface

