!>
!! @brief turbulent diagnosis for LES physics 
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
!! first implementation  by Anurag Dipankar, MPIM (2014-01)
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

MODULE mo_turbulent_diagnostic


  USE mo_kind,               ONLY: wp
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_model_domain,       ONLY: t_patch
  USE mo_run_config,         ONLY: msg_level, iqv, iqc, iqi, iqr, iqs, dtime
  USE mo_nonhydro_types,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,      ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_lnd_types,      ONLY: t_lnd_prog, t_lnd_diag 
  USE mo_parallel_config,    ONLY: nproma
  USE mo_statistics,         ONLY: levels_horizontal_mean
  USE mo_les_nml,            ONLY: avg_interval_sec, turb_profile_list, &
                                   sampl_freq_sec, turb_tseries_list, &
                                   expname, ldiag_les_out
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_write_netcdf      
  USE mo_impl_constants,     ONLY: min_rlcell, min_rlcell_int
  USE mo_physical_constants,  ONLY: cpd, rcvd, p0ref, grav, rcpd, alv

 
  IMPLICIT NONE

  LOGICAL, ALLOCATABLE :: is_at_full_level(:)

  INTEGER  :: nrec_tseries,   nrec_profile
  INTEGER  :: fileid_tseries, fileid_profile
  INTEGER  :: sampl_freq_step, avg_interval_step
  LOGICAL  :: is_sampling_time, is_writing_time
  REAL(wp) :: time_wt

  !Some indices: think of better way
  INTEGER  :: idx_sgs_th_flx, idx_sgs_qv_flx, idx_sgs_qc_flx
  
  CHARACTER(20) :: tname     = 'time'
  CHARACTER(20) :: tlongname = 'Time'

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC  :: calculate_turbulent_diagnostics, write_vertical_profiles, write_time_series
  PUBLIC  :: init_les_turbulent_output, close_les_turbulent_output
  PUBLIC  :: sampl_freq_step, avg_interval_step, is_sampling_time, is_writing_time
  PUBLIC  :: time_wt, idx_sgs_th_flx, idx_sgs_qv_flx, idx_sgs_qc_flx

CONTAINS

  !>
  !! <Calculates 1D and 0D turbulent diagnostics>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !!
  SUBROUTINE calculate_turbulent_diagnostics(             &
                            & p_patch,                    & !in
                            & p_prog,   p_prog_rcf,       & !in
                            & p_diag,                     & !in
                            & p_prog_land, p_diag_land,   & !in
                            & prm_diag                )     !inout
                            

    !>
    ! !INPUT PARAMETERS:

    TYPE(t_patch),   TARGET, INTENT(in)   :: p_patch    !<grid/patch info.
    TYPE(t_nh_diag), TARGET, INTENT(in)   :: p_diag     !<the diagnostic variables
    TYPE(t_nh_prog), TARGET, INTENT(in)   :: p_prog     !<the prognostic variables
    TYPE(t_nh_prog), TARGET, INTENT(in)   :: p_prog_rcf !<the prognostic variables (with
                                                        !< red. calling frequency for tracers!

    TYPE(t_lnd_prog),        INTENT(in)   :: p_prog_land
    TYPE(t_lnd_diag),        INTENT(in)   :: p_diag_land
    TYPE(t_nwp_phy_diag)   , INTENT(inout):: prm_diag

    ! Local
  
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)  :: var3df, var3dh, theta, w_mc
    REAL(wp), ALLOCATABLE, DIMENSION(:)   :: &
              umean, vmean, thmean, qvmean, qcmean, wmean, outvar, thvmean
    REAL(wp) :: outvar0d

    ! Local array bounds:

    INTEGER :: nlev, nlevp1            !< number of full levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    INTEGER :: jc,jk,jb,jg             !block index
    INTEGER :: nvar, n
    CHARACTER(len=*), PARAMETER :: routine = 'mo_turbulent_diagnostic:calculate_turbulent_diagnostics'

    IF(msg_level>18) & 
      CALL message(routine,'Start!')

    nlev       = p_patch%nlev
    nlevp1     = nlev + 1
    
    !allocation
    ALLOCATE( var3df(nproma,nlev,p_patch%nblks_c), var3dh(nproma,nlevp1,p_patch%nblks_c), &
              theta(nproma,nlev,p_patch%nblks_c),  w_mc(nproma,nlev,p_patch%nblks_c), &
              outvar(nlevp1) )

    i_nchdom  = MAX(1,p_patch%n_childdom)

    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    !Get rho at interface levels
!ICON_OMP_PARALLEL 
!ICON_OMP_DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 1 , nlev
         DO jc = i_startidx, i_endidx
           w_mc(jc,jk,jb) = ( p_prog%w(jc,jk,jb) + p_prog%w(jc,jk+1,jb) ) * 0.5_wp
         END DO
       END DO
    END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

!======================================================================================
                 !Some vertical profiles
!======================================================================================
    
    time_wt = sampl_freq_sec / avg_interval_sec

    nvar = SIZE(turb_profile_list,1)

    !Loop over all variables
    DO n = 1 , nvar

     SELECT CASE (TRIM(turb_profile_list(n)))

     CASE('u')

       ALLOCATE(umean(1:nlev))
       CALL levels_horizontal_mean(p_diag%u, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       umean = outvar(1:nlev)

     CASE('v')

       ALLOCATE(vmean(1:nlev))
       CALL levels_horizontal_mean(p_diag%v, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       vmean = outvar(1:nlev)

     CASE('w')

       ALLOCATE(wmean(1:nlev))
       CALL levels_horizontal_mean(w_mc, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       wmean = outvar(1:nlev)

     CASE('thv')

       ALLOCATE(thvmean(1:nlev))
       CALL levels_horizontal_mean(p_prog%theta_v, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       thvmean = outvar(1:nlev)

     CASE('th') !theta mean

       ALLOCATE(thmean(1:nlev))
!ICON_OMP_PARALLEL 
!ICON_OMP_DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
        DO jk = 1 , nlev
         DO jc = i_startidx, i_endidx
           theta(jc,jk,jb)  = p_diag%temp(jc,jk,jb)/p_prog%exner(jc,jk,jb)
         END DO
        END DO
      END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL

       CALL levels_horizontal_mean(theta, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       thmean = outvar(1:nlev)

     CASE('exner')

       CALL levels_horizontal_mean(p_prog%exner, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('rho')

       CALL levels_horizontal_mean(p_prog%rho, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('qv')

       ALLOCATE(qvmean(1:nlev))
       CALL levels_horizontal_mean(p_prog_rcf%tracer(:,:,:,iqv),  &
             p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       qvmean = outvar(1:nlev) 

     CASE('qc')

       ALLOCATE(qcmean(1:nlev))
       CALL levels_horizontal_mean(p_prog_rcf%tracer(:,:,:,iqc),  &
            p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       qcmean = outvar(1:nlev)
       
     CASE('wu')

       IF(ALLOCATED(wmean).AND.ALLOCATED(umean))THEN
!ICON_OMP_PARALLEL 
!ICON_OMP_DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (w_mc(jc,jk,jb)-wmean(jk))*(p_diag%u(jc,jk,jb)-umean(jk))*p_prog%rho(jc,jk,jb) 
            END DO
          END DO
        END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
       ELSE
         CALL finish(routine,'put <wu> after <w> and <u> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('wv')!At full levels

       IF(ALLOCATED(wmean).AND.ALLOCATED(vmean))THEN
!ICON_OMP_PARALLEL 
!ICON_OMP_DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (w_mc(jc,jk,jb)-wmean(jk))*(p_diag%v(jc,jk,jb)-vmean(jk))*p_prog%rho(jc,jk,jb) 
            END DO
          END DO
        END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
       ELSE
         CALL finish(routine,'put <wv> after <w> and <v> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('wth')

       IF(ALLOCATED(wmean).AND.ALLOCATED(thmean))THEN
!ICON_OMP_PARALLEL 
!ICON_OMP_DO PRIVATE(jb,jc,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (w_mc(jc,jk,jb)-wmean(jk))*(theta(jc,jk,jb)-thmean(jk))*p_prog%rho(jc,jk,jb) 
            END DO
          END DO
        END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
       ELSE
         CALL finish(routine,'put <wth> after <w> and <th> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       outvar = outvar * cpd           
  
     CASE('wthv')

       IF(ALLOCATED(wmean).AND.ALLOCATED(thvmean))THEN
!ICON_OMP_PARALLEL 
!ICON_OMP_DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (w_mc(jc,jk,jb)-wmean(jk))*(p_prog%theta_v(jc,jk,jb)-thvmean(jk))*p_prog%rho(jc,jk,jb) 
            END DO
          END DO
        END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
       ELSE
         CALL finish(routine,'put <wthv> after <w> and <thv> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       outvar = outvar * cpd           

     CASE('wqv')

       IF(ALLOCATED(wmean).AND.ALLOCATED(qvmean))THEN
!ICON_OMP_PARALLEL 
!ICON_OMP_DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (w_mc(jc,jk,jb)-wmean(jk))*(p_prog_rcf%tracer(jc,jk,jb,iqv)-qvmean(jk))*p_prog%rho(jc,jk,jb) 
            END DO
          END DO
        END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
       ELSE
         CALL finish(routine,'put <wqv> after <w> and <qv> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       outvar = outvar * alv

     CASE('wqc')

       IF(ALLOCATED(wmean).AND.ALLOCATED(qcmean))THEN
!ICON_OMP_PARALLEL 
!ICON_OMP_DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (w_mc(jc,jk,jb)-wmean(jk))*(p_prog_rcf%tracer(jc,jk,jb,iqc)-qcmean(jk))*p_prog%rho(jc,jk,jb) 
            END DO
          END DO
        END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
       ELSE
         CALL finish(routine,'put <wqc> after <w> and <qc> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       outvar = outvar * alv

     CASE('ww')

       IF(ALLOCATED(wmean))THEN
!ICON_OMP_PARALLEL 
!ICON_OMP_DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (w_mc(jc,jk,jb)-wmean(jk))**2 
            END DO
          END DO
        END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
       ELSE
         CALL finish(routine,'put <ww> after <w> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('thth')

       IF(ALLOCATED(thmean))THEN
!ICON_OMP_PARALLEL 
!ICON_OMP_DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (theta(jc,jk,jb)-thmean(jk))**2 
            END DO
          END DO
        END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
       ELSE
         CALL finish(routine,'put <thth> after <th> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('qvqv')

       IF(ALLOCATED(qvmean))THEN
!ICON_OMP_PARALLEL 
!ICON_OMP_DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (p_prog_rcf%tracer(jc,jk,jb,iqv)-qvmean(jk))**2 
            END DO
          END DO
        END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
       ELSE
         CALL finish(routine,'put <qvqv> after <qv> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('qcqc')

       IF(ALLOCATED(qcmean))THEN
!ICON_OMP_PARALLEL 
!ICON_OMP_DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (p_prog_rcf%tracer(jc,jk,jb,iqc)-qcmean(jk))**2 
            END DO
          END DO
        END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
       ELSE
         CALL finish(routine,'put <qcqc> after <qc> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('uu')

       IF(ALLOCATED(umean))THEN
!ICON_OMP_PARALLEL 
!ICON_OMP_DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (p_diag%u(jc,jk,jb)-umean(jk))**2 
            END DO
          END DO
        END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
       ELSE
         CALL finish(routine,'put <uu> after <u> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('vv')

       IF(ALLOCATED(vmean))THEN
!ICON_OMP_PARALLEL 
!ICON_OMP_DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (p_diag%v(jc,jk,jb)-vmean(jk))**2 
            END DO
          END DO
        END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
       ELSE
         CALL finish(routine,'put <vv> after <v> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('kh')

       CALL levels_horizontal_mean(prm_diag%tkvh, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlevp1))

     CASE('km')

       CALL levels_horizontal_mean(prm_diag%tkvm, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlevp1))

     CASE DEFAULT !In case calculations are performed somewhere else
      
       outvar = 0._wp
       
     END SELECT

     !Calculate time mean
     IF(is_at_full_level(n))THEN
       prm_diag%turb_diag_1dvar(1:nlev,n) = prm_diag%turb_diag_1dvar(1:nlev,n)+outvar(1:nlev)*time_wt
     ELSE
       prm_diag%turb_diag_1dvar(1:nlevp1,n) = prm_diag%turb_diag_1dvar(1:nlevp1,n)+outvar(1:nlevp1)*time_wt
     END IF

    END DO!nvar

!======================================================================================
       !Some time series
!======================================================================================

    nvar = SIZE(turb_tseries_list,1)

    !Loop over all variables
    DO n = 1 , nvar

     SELECT CASE (TRIM(turb_tseries_list(n)))

     CASE('ccover')
       CALL levels_horizontal_mean(prm_diag%clct, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('shflx')
       CALL levels_horizontal_mean(prm_diag%shfl_s, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('lhflx')
       CALL levels_horizontal_mean(prm_diag%lhfl_s, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('ustress')
       CALL levels_horizontal_mean(prm_diag%umfl_s, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('vstress')
       CALL levels_horizontal_mean(prm_diag%vmfl_s, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('tsfc')
       CALL levels_horizontal_mean(p_prog_land%t_g, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('qsfc')
       CALL levels_horizontal_mean(p_diag_land%qv_s, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('hbl')
       CALL levels_horizontal_mean(prm_diag%z_pbl, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('psfc')
       CALL levels_horizontal_mean(p_diag%pres_sfc, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     END SELECT  

     prm_diag%turb_diag_0dvar(n) = outvar0d 
 
    END DO
 
     
    DEALLOCATE( umean, vmean, thmean, qvmean, qcmean, wmean, outvar, var3df, var3dh, theta, w_mc )

    IF(msg_level>18) & 
      CALL message(routine,'Over!')

  END SUBROUTINE calculate_turbulent_diagnostics

!===========================================================================
!>
  !! <write out profile>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !!
  SUBROUTINE write_vertical_profiles(outvar, sim_time)
    REAL(wp), INTENT(IN) :: outvar(:,:), sim_time

    INTEGER :: nvar, n
 
    !Write profiles
    nvar = SIZE(turb_profile_list,1)

    !Loop over all variables
    IF( my_process_is_stdio() )THEN

      !First write time
      CALL writevar_nc(fileid_profile, tname, sim_time, nrec_profile) 

      DO n = 1 , nvar       
       CALL writevar_nc(fileid_profile, TRIM(turb_profile_list(n)), outvar(:,n), nrec_profile) 
      END DO

    END IF

  END SUBROUTINE write_vertical_profiles

!===========================================================================
!>
  !! <write out time series>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !!
  SUBROUTINE write_time_series(outvar, sim_time)
    REAL(wp), INTENT(IN) :: outvar(:), sim_time

    INTEGER :: nvar, n

    !Write time series 
    nvar = SIZE(turb_tseries_list,1)

    !Loop over all variables
    IF( my_process_is_stdio() )THEN

      !First write time
      CALL writevar_nc(fileid_tseries, tname, sim_time, nrec_tseries) 

      DO n = 1 , nvar
        CALL writevar_nc(fileid_tseries, turb_tseries_list(n), outvar(n), nrec_tseries) 
      END DO

    END IF

  END SUBROUTINE write_time_series

!===========================================================================
!>
  !! <initialize turbulent output>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !!
  SUBROUTINE init_les_turbulent_output(p_patch, p_metrics, time, ldelete)
   TYPE(t_patch),   TARGET, INTENT(in)   :: p_patch    !<grid/patch info.
   TYPE(t_nh_metrics)     , INTENT(in)   :: p_metrics
   REAL(wp), INTENT(IN)                  :: time
   LOGICAL, INTENT(IN), OPTIONAL         :: ldelete
  
   CHARACTER (40), ALLOCATABLE, DIMENSION(:) :: dimname, dimlongname, dimunit
   CHARACTER (LEN=80)                        :: longname, unit
   REAL(wp), ALLOCATABLE                     :: dimvalues(:,:)
   INTEGER,  ALLOCATABLE                     :: dimsize(:)
   INTEGER :: n, nlev, nlevp1, nvar
   REAL(wp) :: z_mc_avg(p_patch%nlev), z_ic_avg(p_patch%nlev+1)
   CHARACTER(len=*), PARAMETER :: routine = 'mo_turbulent_diagnostic:init_les_turbulent_output'
 
   !Check if diagnostics are to be calculated or not
   IF(.NOT.ldiag_les_out)THEN
    sampl_freq_step   = 0
    avg_interval_step = 0
    RETURN
   END IF

   IF(msg_level>18)CALL message(routine,'Start!')

   !Sampling and output frequencies in terms of time steps
   sampl_freq_step   = MAX(1,NINT(sampl_freq_sec/dtime))
   avg_interval_step = MAX(1,NINT(avg_interval_sec/dtime))

   nlev   = p_patch%nlev
   nlevp1 = nlev + 1

   !Dimensions
   ALLOCATE( dimname(2), dimlongname(2), dimunit(2), dimsize(2), dimvalues(nlevp1,2) )

   !Calculate average height
   CALL levels_horizontal_mean(p_metrics%z_mc, p_patch%cells%area, p_patch%cells%owned,  z_mc_avg)
   CALL levels_horizontal_mean(p_metrics%z_ifc, p_patch%cells%area, p_patch%cells%owned, z_ic_avg)

   !open profile file
   IF( my_process_is_stdio() ) &
     CALL open_nc(TRIM(expname)//'_profile.nc', fileid_profile, nrec_profile, time, ldelete)
 
   !addvar
   nvar = SIZE(turb_profile_list,1)

   ALLOCATE(is_at_full_level(nvar))

   DO n = 1 , nvar

     is_at_full_level(n) = .TRUE.

     SELECT CASE (TRIM(turb_profile_list(n)))
    
     CASE('u')
      longname = 'zonal wind'
      unit     = 'm/s'
     CASE('v')
      longname = 'meridional wind'
      unit     = 'm/s'
     CASE('w')
      longname = 'vertical wind'
      unit     = 'm/s'
     CASE('th') !theta mean
      longname = 'potential temperature'
      unit     = 'K'
     CASE('thv') !thetav mean
      longname = 'virtual potential temperature'
      unit     = 'K'
     CASE('exner')
      longname = 'exner function'
      unit     = ''
     CASE('rho')
      longname = 'density'
      unit     = 'kg/m3'
     CASE('qv')
      longname = 'specific humidity'
      unit     = 'kg/kg'
     CASE('qc')
      longname = 'cloud water'
      unit     = 'kg/kg'
     CASE('wu')
       longname = 'resolved zonal wind flux'
       unit     = 'm2/s2'
     CASE('wv')
       longname = 'resolved meridional wind flux'
       unit     = 'm2/s2'
     CASE('wth')
       longname = 'resolved potential temperature flux'
       unit     = 'W/m2'
     CASE('wthv')
       longname = 'resolved virtual potential temperature flux'
       unit     = 'W/m2'
     CASE('wqv')
       longname = 'resolved specific humidity flux'
       unit     = 'W/m2'
     CASE('wqc')
       longname = 'resolved cloud water flux'
       unit     = 'W/m2'
     CASE('ww')
       longname = 'resolved vertical velocity variance'
       unit     = 'm2/s2'
     CASE('thth')
       longname = 'resolved potential temperature variance'
       unit     = 'K2'
     CASE('qvqv')
       longname = 'resolved specific humidity variance'
       unit     = 'kg2/kg2'
     CASE('qcqc')
       longname = 'resolved cloud water variance'
       unit     = 'kg2/kg2'
     CASE('uu')
       longname = 'resolved zonal wind variance'
       unit     = 'm2/s2'
     CASE('vv')
       longname = 'resolved meridional wind variance'
       unit     = 'm2/s2'
     CASE('kh')
       longname = '(mass) eddy diffusivity'
       unit     = 'kg/ms'
       is_at_full_level(n) = .FALSE.
     CASE('km')
       longname = '(mass) eddy viscosity'
       unit     = 'kg/ms'
       is_at_full_level(n) = .FALSE.
     CASE('wthd') !diffuse flux
       longname = 'subgrid potential temperature flux'
       unit     = 'W/m2'
       is_at_full_level(n) = .FALSE.
       idx_sgs_th_flx = n
     CASE('wqvd') !diffuse flux
       longname = 'subgrid specific humidity flux'
       unit     = 'W/m2'
       is_at_full_level(n) = .FALSE.
       idx_sgs_qv_flx = n
     CASE('wqcd') !diffuse flux
       longname = 'subgrid cloud water flux'
       unit     = 'W/m2'
       is_at_full_level(n) = .FALSE.
       idx_sgs_qc_flx = n
     END SELECT

     dimname(2) = tname
     dimlongname(2) = tlongname
     dimunit(1) = 'm'
     dimunit(2) = 's'
     dimvalues = 0._wp
     IF(is_at_full_level(n))THEN
      dimname(1) = 'zf'
      dimlongname(1) = 'Full level height'
      dimsize = (/nlev,0/)
      dimvalues(1:nlev,1) = z_mc_avg(1:nlev)     
     ELSE
      dimname(1) = 'zh'
      dimlongname(1) = 'Half level height'
      dimsize = (/nlevp1,0/)
      dimvalues(1:nlevp1,1) = z_ic_avg(1:nlevp1)     
     END IF

     IF( my_process_is_stdio() ) &
       CALL addvar_nc(fileid_profile, TRIM(turb_profile_list(n)), TRIM(longname), TRIM(unit), &
                      dimname, dimlongname, dimunit, dimsize, dimvalues)

    END DO!nvar
    

    !deallocate
    DEALLOCATE( dimname, dimlongname, dimunit, dimsize, dimvalues )
    ALLOCATE( dimname(1), dimlongname(1), dimunit(1) )


   !open time series file
   IF( my_process_is_stdio() ) &
      CALL open_nc(TRIM(expname)//'_tseries.nc', fileid_tseries, nrec_tseries, time, ldelete)
 
   !addvar
   nvar = SIZE(turb_tseries_list,1)

   DO n = 1 , nvar

     SELECT CASE (TRIM(turb_tseries_list(n)))
 
     CASE('ccover')
       longname = 'cloud cover'
       unit     = ' '
     CASE('shflx')
       longname = 'surface sensible heat flux'
       unit     = 'W/m2'
     CASE('lhflx')
       longname = 'surface latent heat flux'
       unit     = 'W/m2'
     CASE('ustress')
       longname = 'surface zonal stress'
       unit     = 'Kg/ms2'
     CASE('vstress')
       longname = 'surface meridional stress'
       unit     = 'Kg/ms2'
     CASE('tsfc')
       longname = 'surface temperature'
       unit     = 'K'
     CASE('qsfc')
       longname = 'surface humidity'
       unit     = 'kg/kg'
     CASE('hbl')
       longname = 'boundary layer height'
       unit     = 'm'
     CASE('psfc')
       longname = 'surface pressure'
       unit     = 'Pa'
     END SELECT  
  
     dimname(1) = tname
     dimlongname(1) = tlongname
     dimunit(1) = 's'
     IF( my_process_is_stdio() ) &
        CALL addvar_nc(fileid_tseries, TRIM(turb_tseries_list(n)), TRIM(longname), TRIM(unit), &
                       dimname, dimlongname, dimunit)

   END DO!nvar

   DEALLOCATE( dimname, dimlongname, dimunit )

   IF(msg_level>18)CALL message(routine,'Over!')

  END SUBROUTINE init_les_turbulent_output

!===========================================================================
!>
  !! <close turbulent output>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !!
  SUBROUTINE close_les_turbulent_output

   IF(.NOT.ldiag_les_out)THEN
    RETURN
   END IF

   IF( my_process_is_stdio() ) THEN
     CALL close_nc(fileid_profile) 
     CALL close_nc(fileid_tseries) 
   END IF
   DEALLOCATE(is_at_full_level)

  END SUBROUTINE close_les_turbulent_output

!===========================================================================

END MODULE mo_turbulent_diagnostic

