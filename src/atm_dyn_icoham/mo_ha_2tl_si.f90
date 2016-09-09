!>
!! module for two-time-level semi-implicit time stepping scheme
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! Initial version by Hui Wan, MPI-M (2009-11)
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
MODULE mo_ha_2tl_si

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_physical_constants,  ONLY: rd, rcpd
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: nlev, nlevp1, msg_level
  USE mo_dynamics_config,     ONLY: lshallow_water 
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp,                ONLY: cells2edges_scalar
  USE mo_math_gradients,     ONLY: grad_fd_norm
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_ha_dynamics,         ONLY: dyn_temp, continuity,                   &
    &                               energy_conversion_terms
  USE mo_ha_dynamics_adv,     ONLY: temp_adv_vertical, temp_adv_horizontal, &
    &                               vn_adv_vertical,   vn_adv_horizontal
  USE mo_ha_diag_util,        ONLY: update_diag_state
  USE mo_vertical_coord_table,ONLY: vct_b
  USE mo_exception,           ONLY: message,finish
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, EULER_FORWARD, AB2
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_ha_2tl_si_solver,    ONLY: solver
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, sync_patch_array
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, timer_step_2tl_si

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: step_2tl_si

  CONTAINS
  !>
  !! main subroutine of the two-time-level semi-implicit time stepping scheme
  !!
  !! <Describe the function of the routine and algorithm(s) used in the routine>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !! Initial version by Hui Wan, MPI-M (2009-11-06)
  !!
  SUBROUTINE step_2tl_si( si_expl_scheme, si_2tls, si_rtol, &
                          ltheta_dyn,                       &
                          p_dtime, p_patch, p_int_state,    &
                          p_now, p_ext_data,                &
                          p_new, p_diag, p_tend_dyn )

  CHARACTER(len=*),PARAMETER :: routine = 'mo_ha_2tl_si:step_2tl_si'

  !! arguments

  INTEGER ,INTENT(IN) :: si_expl_scheme
  REAL(wp),INTENT(IN) :: si_2tls, si_rtol
  LOGICAL, INTENT(IN) :: ltheta_dyn
  REAL(wp),INTENT(IN) :: p_dtime

  TYPE(t_patch), TARGET, INTENT(INOUT) :: p_patch
  TYPE(t_int_state),     INTENT(IN)    :: p_int_state
  TYPE(t_hydro_atm_prog),INTENT(IN)    :: p_now
  TYPE(t_external_data), INTENT(INOUT) :: p_ext_data   !< external data
  TYPE(t_hydro_atm_prog),INTENT(INOUT) :: p_new
  TYPE(t_hydro_atm_diag),INTENT(INOUT) :: p_diag
  TYPE(t_hydro_atm_prog),INTENT(INOUT) :: p_tend_dyn

  !! local variables

  REAL(wp),DIMENSION(nproma,nlev,p_patch%nblks_c) :: z_dtemp_expl, z_dtemp
  REAL(wp),DIMENSION(nproma,nlev,p_patch%nblks_e) :: z_dvn_expl,   z_dvn, z_rhs
  REAL(wp),DIMENSION(nproma,     p_patch%nblks_c) :: z_dps_expl,   z_dps

  REAL(wp),DIMENSION(nproma,nlev,p_patch%nblks_e) :: z_cgradps,z_tmp_e
  REAL(wp),DIMENSION(nproma,nlev,p_patch%nblks_c) :: z_tmp_c

  INTEGER :: nblks_c, nblks_e, npromz_e
  INTEGER :: jb, jk, jbs,is,ie

  REAL(wp) :: zconst, z_dtsq

! for the GMRES solver

  INTEGER, PARAMETER :: nmax_iter = 100
  LOGICAL  :: lmaxiter
  INTEGER  :: niter
  REAL(wp) :: z_residual(nmax_iter)
  CHARACTER(LEN=MAX_CHAR_LENGTH) :: string

!----------------------------------------------------------------------------
   IF (ltimer) CALL timer_start(timer_step_2tl_si)

   nblks_c   = p_patch%nblks_c
   nblks_e   = p_patch%nblks_e
   npromz_e  = p_patch%npromz_e

!----------------------------------------------------------
! 1. Calculate right-hand side of the increment equations
!----------------------------------------------------------
 IF (ltheta_dyn) CALL finish(TRIM(routine),'theta dyn. not implimented')

   CALL tend_expl( si_expl_scheme, p_dtime,                &! in
                   p_now, p_patch, p_int_state,            &! in
                   p_ext_data,                             &! in
                   p_diag, p_tend_dyn,                     &! inout
                   z_dvn_expl, z_dtemp_expl, z_dps_expl  )  ! out

! 1.3 Calculate the prefactor of the pressure gradient term.
!     It is needed by subroutine pgrad_vn

   jbs = p_patch%cells%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,nblks_c
         CALL get_indices_c( p_patch, jb,jbs,nblks_c, is,ie, 2)
         z_tmp_c(is:ie,:,jb) = p_diag%tempv(is:ie,:,jb)   &
                              /p_diag%pres_mc(is:ie,:,jb)
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   CALL cells2edges_scalar (z_tmp_c, p_patch,    &! in
                            p_int_state%c_lin_e, &! in
                            z_tmp_e              )! out

   jbs = p_patch%edges%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk,zconst) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_e
      CALL get_indices_e( p_patch, jb,jbs,nblks_e, is,ie, 2)

      DO jk = 1,nlev
         zconst = rd*0.5_wp*(vct_b(jk)+vct_b(jk+1))
         z_cgradps(is:ie,jk,jb) = zconst*z_tmp_e(is:ie,jk,jb)
      ENDDO
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

! 1.3 The changes in wind consist two parts:
!     - the explicit part from the velocity equation, and
!     - the changes due to z_dtps

   CALL pgrad_vn( z_dtemp_expl, z_dps_expl,                     &
                  p_diag%rdlnpr_c, p_diag%rdalpha_c, z_cgradps, &
                  p_patch, z_tmp_e )

   jbs = p_patch%edges%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_e
      CALL get_indices_e( p_patch, jb,jbs,nblks_e, is,ie, 2)
      z_rhs(is:ie,:,jb) = -p_dtime*si_2tls*z_tmp_e(is:ie,:,jb) &
                         + z_dvn_expl(is:ie,:,jb)
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!! CALL sync_patch_array(SYNC_E, p_patch, z_rhs)

!---------------------------------------------------------
! 2. Solve for wind increment
!---------------------------------------------------------

   z_dtsq = p_dtime*p_dtime*si_2tls*si_2tls

!$OMP PARALLEL
!$OMP WORKSHARE
   z_dvn(:,:,:) = 0._wp  ! first guess
!$OMP END WORKSHARE
!$OMP END PARALLEL

   CALL solver( z_dvn(:,:,:),       &! x. Input is the first guess
               lhs_dvn_eqn,         &! sbr. calculating l.h.s.
               p_now%temp,          &! used for calculating l.h.s.
               p_diag%rdelp_c,      &! used for calculating l.h.s.
               p_diag%delp_e,       &! used for calculating l.h.s.
               p_diag%rdlnpr_c,     &! used for calculating l.h.s.
               p_diag%rdalpha_c,    &! used for calculating l.h.s.
               z_cgradps,           &! used for calculating l.h.s.
               p_patch,             &! used for calculating l.h.s.
               p_int_state,         &! interpolation state
               nblks_e,             &! number of blocks
               npromz_e,            &! length of last block
               z_dtsq,              &! used for calculating l.h.s.
               z_rhs(:,:,:),        &! right hand side as input
               si_rtol,             &! relative tolerance
               .FALSE.,             &! NOT absolute tolerance
               nmax_iter,           &! max. # of iterations to do
               lmaxiter,            &! out: .true. = not converged
               niter,               &! out: # of iterations done
               z_residual           &! out: the residual (array)
               )

   IF (lmaxiter) THEN
      CALL finish('GMRES solver: ','NOT YET CONVERGED !!')
   ENDIF
   IF (msg_level >= 1) THEN
     WRITE(string,'(a,i4,a,e20.10)') 'GMRES solver: iteration ', niter,  &
                                   ', residual = ', ABS(z_residual(niter))
     CALL message(TRIM(routine),TRIM(string))
   ENDIF

!---------------------------------------------------------
! 3. Update prognostic variables
!---------------------------------------------------------
! 3.1 Wind

  jbs = p_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = jbs,nblks_e
     CALL get_indices_e( p_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1)
     p_new%vn(is:ie,:,jb) = p_now%vn(is:ie,:,jb) + z_dvn(is:ie,:,jb)
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

! 3.2 Temperature and surface pressure

  CALL conteq_vn( z_dvn, p_diag%delp_e, p_now%temp, p_diag%rdelp_c, &
                  p_diag%rdlnpr_c, p_diag%rdalpha_c,                &
                  p_patch, p_int_state, z_dtemp, z_dps )

  jbs = p_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = jbs,nblks_c
     CALL get_indices_c( p_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)

     p_new%temp(is:ie,:,jb) = p_now%temp(is:ie,:,jb)    &
                            + z_dtemp_expl(is:ie,:,jb)  &
                            + si_2tls*p_dtime*z_dtemp(is:ie,:,jb)

     p_new%pres_sfc(is:ie,jb) = p_now%pres_sfc(is:ie,jb)    &
                              + z_dps_expl(is:ie,jb)        &
                              + si_2tls*p_dtime*z_dps(is:ie,jb)
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  IF (ltimer) CALL timer_stop(timer_step_2tl_si)

  END SUBROUTINE step_2tl_si

  !>
  !! Calculate the explicit terms of the governing equations
  !! for the 2-time-level semi-implicit correction scheme
  !!
  SUBROUTINE tend_expl( si_expl_scheme,                          &! in
                        p_dtime, pt_now, pt_patch, pt_int_state, &! in
                        pt_ext_data,                             &! in
                        pt_diag, pt_tend_save,                   &! inout
                        p_dvn, p_dtemp, p_dps )                   ! out
  !! Arguments

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_ha_2tl_si:tend_expl'

  INTEGER ,INTENT(IN) :: si_expl_scheme 
  REAL(wp),INTENT(IN) :: p_dtime

  TYPE(t_patch), TARGET, INTENT(INOUT) :: pt_patch
  TYPE(t_int_state),TARGET, INTENT(IN) :: pt_int_state

  TYPE(t_hydro_atm_prog),INTENT(IN)    :: pt_now        !< prognostic variables
  TYPE(t_external_data),   INTENT(INOUT) :: pt_ext_data   !< external data

  TYPE(t_hydro_atm_diag),INTENT(INOUT) :: pt_diag       !< diagnostic variables
  TYPE(t_hydro_atm_prog),INTENT(INOUT) :: pt_tend_save  !< for the simple
  !< Euler forward scheme, this tendency state stores the tendencies of the
  !< prognostic variables at the current time step;
  !< For the 2nd order Adams-Bashforth scheme
  !< it contains only the "slow" component.
  !< The input is at time step n-1; The output is at step n.

  REAL(wp),INTENT(OUT) :: p_dvn   (:,:,:)
  REAL(wp),INTENT(OUT) :: p_dtemp (:,:,:)
  REAL(wp),INTENT(OUT) :: p_dps   (:,:)

  !! Local variables

  INTEGER :: nblks_c, nblks_e, jb,jbs,is,ie

  REAL(wp) :: z_mdiv_int (nproma,nlevp1,pt_patch%nblks_c)
  REAL(wp) :: z_mdiv     (nproma,nlev  ,pt_patch%nblks_c)

  REAL(wp) :: z_ddt_temp_slow (nproma,nlev,pt_patch%nblks_c)
  REAL(wp) :: z_ddt_temp_fast (nproma,nlev,pt_patch%nblks_c)
  REAL(wp) :: z_ddt_vn_slow   (nproma,nlev,pt_patch%nblks_e)
  REAL(wp) :: z_ddt_vn_fast   (nproma,nlev,pt_patch%nblks_e)

  INTEGER      :: ischeme
  LOGICAL,SAVE :: lfirst_step = .TRUE.

  CHARACTER(LEN=MAX_CHAR_LENGTH) :: string

! Dimension parameters

   nblks_c   = pt_patch%nblks_c
   nblks_e   = pt_patch%nblks_e

! Start calculation

  IF (lfirst_step) THEN
     ischeme = EULER_FORWARD
     lfirst_step  = .FALSE.
  ELSE
     ischeme = si_expl_scheme
  ENDIF

  WRITE(string,'(a,i2)') 'si slow comp =',ischeme
  CALL message(TRIM(routine),TRIM(string))

  SELECT CASE(ischeme)
  CASE(EULER_FORWARD)
!=============================================================================
! Use the Eurler forward scheme for the slow components
!=============================================================================

   CALL dyn_temp( pt_patch, pt_int_state, pt_now, pt_ext_data,   &! in
                  pt_diag, pt_tend_save                        )  ! out

! Multiply by time step to get the increments
!---------------------------------------------
! Surface pressure

!$OMP PARALLEL PRIVATE(jbs)
   jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_c
      CALL get_indices_c( pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)
      p_dps(is:ie,jb) = p_dtime*pt_tend_save%pres_sfc(is:ie,jb)
   ENDDO
!$OMP END DO

!-------------
! Temperature
  IF (.NOT.lshallow_water) THEN
     jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
     DO jb = jbs,nblks_c
       CALL get_indices_c( pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)
       p_dtemp(is:ie,:,jb) = p_dtime*pt_tend_save%temp(is:ie,:,jb)
     ENDDO
!$OMP END DO
  ENDIF

!----------
! Velocity

   jbs = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_e
      CALL get_indices_e( pt_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1)
      p_dvn(is:ie,:,jb) = p_dtime*pt_tend_save%vn(is:ie,:,jb)
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  CASE(AB2)
!=============================================================================
! Use the 2nd order Adams-Bashforth scheme for the slow components
!=============================================================================
! 1. Update the diagnostic state vector. This includes the calculation of
! pressure and related quantities, vorticity and divergence, u- and v-wind,
! virtual temperature and geopotential.

   CALL update_diag_state( pt_now, pt_patch, pt_int_state, pt_ext_data, &
     &                     pt_diag )

! 2. The continuity equation does not contain any slow component.
! Evaluate the surface pressure tendency at current time level.

   CALL continuity( pt_now%vn, pt_diag%delp_e,               &! in
                    pt_patch, pt_int_state, .TRUE.,          &! in
                    z_mdiv, z_mdiv_int, pt_diag%mass_flux_e, &! out
                    p_dps, pt_diag%weta                    )  ! out

!$OMP PARALLEL PRIVATE(jbs)
   jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_c
      CALL get_indices_c( pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)
      p_dps(is:ie,jb) = p_dtime*p_dps(is:ie,jb)
   ENDDO
!$OMP END DO

! 3. Initialize velocity and temperature tendencies

   jbs = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1)
      z_ddt_vn_slow(is:ie,:,jb) = 0._wp
      z_ddt_vn_fast(is:ie,:,jb) = 0._wp
   ENDDO
!$OMP END DO

   IF (.NOT.lshallow_water) THEN
      jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,nblks_c
         CALL get_indices_c( pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)
         z_ddt_temp_slow(is:ie,:,jb) = 0._wp
!        z_ddt_temp_fast does not need initialization
      ENDDO
!$OMP END DO
   ENDIF
!$OMP END PARALLEL

! 4. Calculate tendency induced by slow processes

   IF (.NOT.lshallow_water) THEN

! 4.1 Vertical advection of momentum

     CALL vn_adv_vertical( pt_now%vn, pt_diag%weta, pt_diag%delp_e, &! in
                           pt_patch,  pt_int_state,                 &! in
                           z_ddt_vn_slow                  )          ! inout

! 4.2 Vertical and horizontal advection of temperature

     CALL temp_adv_vertical( pt_now%temp, pt_diag%weta,       &! in
                             pt_diag%rdelp_c, pt_patch,       &! in
                             z_ddt_temp_slow                )  ! inout

     CALL temp_adv_horizontal( pt_now%temp, pt_diag%mass_flux_e,   &! in
                               pt_diag%rdelp_c, z_mdiv,            &! in
                               pt_patch, pt_int_state,             &! in
                               z_ddt_temp_slow                  )   ! inout
   ENDIF

! 4.3 Horizontal advection and Coriolis force

   CALL vn_adv_horizontal( pt_now%vn,                           &! in
                           pt_diag%rel_vort,                    &! in
                           pt_diag%delp_c,                      &! in
                           pt_patch, pt_int_state,              &! in
                           z_ddt_vn_slow,                       &! inout
                           pt_diag%e_kin, pt_diag%vt,           &! inout
                           pt_diag%delp_v)                       ! inout

! 5. Energy conversion terms include both fast and slow components

   CALL energy_conversion_terms( pt_diag, pt_patch, pt_int_state, &! inout,in,in
                                 z_mdiv, z_mdiv_int,              &! in
                                 z_ddt_temp_slow,                 &! inout
                                 z_ddt_vn_fast,                   &! out
                                 .TRUE.,                          &! in
                                 z_ddt_temp_fast               )   ! out

! 6. Merge the tendency terms, and save the "slow" tendencies of step n

!$OMP PARALLEL PRIVATE(jbs)
   IF (.NOT.lshallow_water) THEN
      jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,nblks_c
        CALL get_indices_c( pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)

        p_dtemp(is:ie,:,jb) = p_dtime* &
                            (-0.5_wp*pt_tend_save%temp(is:ie,:,jb) &
                             +1.5_wp*z_ddt_temp_slow(is:ie,:,jb)           &
                             +       z_ddt_temp_fast(is:ie,:,jb)          )

        pt_tend_save%temp(is:ie,:,jb) = z_ddt_temp_slow(is:ie,:,jb)
      ENDDO
!$OMP END DO

      !LL : this is already calculated on the halos
      CALL sync_patch_array( SYNC_C, pt_patch, pt_tend_save%temp )
   ENDIF

   jbs = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_e
      CALL get_indices_e( pt_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1)

      p_dvn(is:ie,:,jb) = p_dtime* &
                        ( -0.5_wp*pt_tend_save%vn(is:ie,:,jb)  &
                          +1.5_wp*z_ddt_vn_slow(is:ie,:,jb)    &
                          +       z_ddt_vn_fast(is:ie,:,jb)   )

      pt_tend_save%vn(is:ie,:,jb) = z_ddt_vn_slow(is:ie,:,jb)
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   CALL sync_patch_array( SYNC_E, pt_patch, pt_tend_save%vn )

  END SELECT

! Synchronize tendencies

   CALL sync_patch_array( SYNC_C, pt_patch, p_dps  )
   CALL sync_patch_array( SYNC_C, pt_patch, p_dtemp)
   CALL sync_patch_array( SYNC_E, pt_patch, p_dvn  )

  END SUBROUTINE tend_expl
  !>
  !!
  !! @par Revision History
  !! Initial version by Hui Wan (2009-11.10)
  !!
  SUBROUTINE conteq_vn (p_vn,p_delp_e,p_temp,p_rdelp_c,p_rdlnpr,p_rdalpha, &
                        p_patch, p_int_state, p_dtemp, p_dps)

  TYPE(t_patch),    INTENT(INOUT) :: p_patch
  TYPE(t_int_state),INTENT(IN) :: p_int_state

  REAL(wp),INTENT(IN)  :: p_vn      (nproma,nlev,p_patch%nblks_e)
  REAL(wp),INTENT(IN)  :: p_delp_e  (nproma,nlev,p_patch%nblks_e)
  REAL(wp),INTENT(IN)  :: p_temp    (nproma,nlev,p_patch%nblks_c)
  REAL(wp),INTENT(IN)  :: p_rdelp_c (nproma,nlev,p_patch%nblks_c)
  REAL(wp),INTENT(IN)  :: p_rdlnpr  (nproma,nlev,p_patch%nblks_c)
  REAL(wp),INTENT(IN)  :: p_rdalpha (nproma,nlev,p_patch%nblks_c)
  REAL(wp),INTENT(OUT) :: p_dtemp   (nproma,nlev,p_patch%nblks_c)
  REAL(wp),INTENT(OUT) :: p_dps     (nproma,     p_patch%nblks_c)

  REAL(wp) :: z_mdiv     (nproma,nlev  ,p_patch%nblks_c)
  REAL(wp) :: z_mdiv_int (nproma,nlevp1,p_patch%nblks_c)
  REAL(wp), POINTER :: z_mflux(:, :, :)

  INTEGER  :: nblks_c, npromz_c
  INTEGER  :: jb, nlen, jk, return_status

!---------------------------------------
! Dimension parameters

  nblks_c  = p_patch%nblks_c
  npromz_c = p_patch%npromz_c
  
  ALLOCATE( z_mflux(nproma, nlev, p_patch%nblks_e),stat=return_status)
  IF (return_status > 0) &
    CALL finish ("conteq_vn", 'ALLOCATE z_mflux failed')
  
! Mass divergence and its vertical integral

  ! LL: continuity should not need any communication

! (.FALSE. in the argument means do not diagnose vertical velocity)
  CALL continuity( p_vn, p_delp_e, p_patch, p_int_state, .FALSE.,   &! in
                   z_mdiv, z_mdiv_int, z_mflux, p_dps              ) ! out

! Tendency of temperature
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1,nblks_c

    IF (jb /= nblks_c) THEN
        nlen = nproma
    ELSE
        nlen = npromz_c
    ENDIF

    DO jk = 1,nlev
       p_dtemp(1:nlen,jk,jb) = -rcpd*p_temp(1:nlen,jk,jb)  &
                         *p_rdelp_c(1:nlen,jk,jb)          &
                         *( p_rdlnpr(1:nlen,jk,jb)*z_mdiv_int(1:nlen,jk,jb) &
                           +p_rdalpha(1:nlen,jk,jb)*z_mdiv(1:nlen,jk,jb)   )
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
!---------------------------------------
  DEALLOCATE( z_mflux )
  
  END SUBROUTINE conteq_vn
  !>
  !!
  !! @par Revision History
  !! Initial version by Hui Wan (2009-11.10)
  !!
  SUBROUTINE pgrad_vn (p_dtemp, p_dps, p_rdlnpr, p_rdalpha, p_cgradps, &
                        p_patch, p_dvn)

  TYPE(t_patch),INTENT(IN) :: p_patch

  REAL(wp),INTENT(IN)  :: p_dtemp   (nproma,nlev,p_patch%nblks_c)
  REAL(wp),INTENT(IN)  :: p_dps     (nproma,1,   p_patch%nblks_c)
  REAL(wp),INTENT(IN)  :: p_rdlnpr  (nproma,nlev,p_patch%nblks_c)
  REAL(wp),INTENT(IN)  :: p_rdalpha (nproma,nlev,p_patch%nblks_c)
  REAL(wp),INTENT(IN)  :: p_cgradps (nproma,nlev,p_patch%nblks_e)
  REAL(wp),INTENT(OUT) :: p_dvn     (nproma,nlev,p_patch%nblks_e)

  REAL(wp) :: z_tmp_e  (nproma,nlev,p_patch%nblks_e)
  REAL(wp) :: z_tmp_c  (nproma,nlev,p_patch%nblks_c)
  REAL(wp) :: z_gradps (nproma,1   ,p_patch%nblks_e)

  INTEGER  :: nblks_c, npromz_c
  INTEGER  :: nblks_e, npromz_e
  INTEGER  :: jb, nlen, jk, jkm

!---------------------------------------

  nblks_c  = p_patch%nblks_c
  npromz_c = p_patch%npromz_c
  nblks_e  = p_patch%nblks_e
  npromz_e = p_patch%npromz_e

!------------------------------------------------------
! gradient of geopotential change
!------------------------------------------------------
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jkm) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1,nblks_c

     IF (jb /= nblks_c) THEN
         nlen = nproma
     ELSE
         nlen = npromz_c
     ENDIF

     z_tmp_c(1:nlen,nlev,jb) = 0._wp

     DO jk = nlev,2,-1
        jkm = jk - 1
        z_tmp_c(1:nlen,jkm,jb) =  z_tmp_c(1:nlen,jk,jb)  &
                                + p_rdlnpr(1:nlen,jk,jb) &
                                 *p_dtemp(1:nlen,jk,jb)
     ENDDO

     DO jk = nlev,1,-1
        z_tmp_c(1:nlen,jk,jb) =   z_tmp_c(1:nlen,jk,jb)   &
                                + p_rdalpha(1:nlen,jk,jb) &
                                 *p_dtemp(1:nlen,jk,jb)
     ENDDO
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!------------------------------------------------------
! LL: here requires z_tmp_c, and p_dps on halo cells sharing an edge
!     with owned cells

  
  CALL grad_fd_norm( z_tmp_c, p_patch, z_tmp_e, opt_rlstart=4 )

!------------------------------------------------------
! gradient of surface pressure change
!------------------------------------------------------

  CALL grad_fd_norm( p_dps, p_patch, z_gradps,  &
                     opt_slev=1, opt_elev=1, opt_rlstart=4 )
!------------------------------------------------------
! Combine two parts
!------------------------------------------------------
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1,nblks_e

     IF (jb /= nblks_e) THEN
         nlen = nproma
     ELSE
         nlen = npromz_e
     ENDIF

     DO jk = 1,nlev
        p_dvn(1:nlen,jk,jb) =  z_tmp_e(1:nlen,jk,jb)   &
                             + p_cgradps(1:nlen,jk,jb) &
                              *z_gradps(1:nlen,1 ,jb)
     ENDDO
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE pgrad_vn
  
  !>
  !! @par Revision History
  !! Initial version by Hui Wan (2009-11.10)
  !!
  FUNCTION lhs_dvn_eqn ( p_dvn, &
                         p_temp,p_rdelp_c,p_delp_e,p_rdlnpr,p_rdalpha, &
                         p_cgradps, p_patch, p_int_state, p_coeff )    &
  RESULT(p_lhs)

  !! Input information

  REAL(wp),INTENT(IN) :: p_dvn(:,:,:)
  REAL(wp),INTENT(IN) :: p_temp(:,:,:), p_rdelp_c(:,:,:)
  REAL(wp),INTENT(IN) :: p_rdlnpr(:,:,:), p_rdalpha(:,:,:)
  REAL(wp),INTENT(IN) :: p_delp_e(:,:,:)
  REAL(wp),INTENT(IN) :: p_cgradps(:,:,:)
  REAL(wp),INTENT(IN) :: p_coeff

  TYPE(t_patch), INTENT(INOUT)  :: p_patch
  TYPE(t_int_state), INTENT(IN) :: p_int_state

  !! Result

  REAL(wp) :: p_lhs(SIZE(p_dvn,1),SIZE(p_dvn,2),SIZE(p_dvn,3))

  !! Local variables

  REAL(wp) :: z_dtemp(SIZE(p_temp,1),SIZE(p_temp,2),SIZE(p_temp,3))
  REAL(wp) :: z_dps  (SIZE(p_temp,1),               SIZE(p_temp,3))

  INTEGER :: jb, nblks_e, npromz_e, nlen

!----------------------------------------------------

  nblks_e  = p_patch%nblks_e
  npromz_e = p_patch%npromz_e

  ! LL: conteq_vn should not need communication
  CALL conteq_vn( p_dvn, p_delp_e, p_temp, p_rdelp_c,        &
                  p_rdlnpr, p_rdalpha, p_patch, p_int_state, &
                  z_dtemp, z_dps )

  ! LL: pgrad_vn requires one communication point
  CALL pgrad_vn( z_dtemp, z_dps, p_rdlnpr, p_rdalpha, p_cgradps, &
                 p_patch, p_lhs )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1,nblks_e

     IF (jb /= nblks_e) THEN
         nlen = nproma
     ELSE
         nlen = npromz_e
     ENDIF
     p_lhs(1:nlen,:,jb) = p_dvn(1:nlen,:,jb)       &
                        + p_coeff*p_lhs(1:nlen,:,jb)
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END FUNCTION lhs_dvn_eqn

END MODULE mo_ha_2tl_si

