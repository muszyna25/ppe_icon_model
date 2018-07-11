!>
!! Module for explicit Runge-Kutta method.
!!
!! Following pages 9-10 in
!!  Giorgetta et al. (2009): Conservative space and time regularizations
!!  for the ICON model. Reports on Earth System Science 67, Max Planck
!!  Institute for Meteorology, Hamburg, Germany. Available at
!!  <http://www.mpimet.mpg.de/wissenschaft/publikationen/erdsystemforschung.html#c2612>,
!! we can summarize the explicit Runge-Kutta method as follows:
!!
!! The ordinary differential equation to be solved reads  dU/dt = S(U) .
!! To advance the numerical solution from time level t^n to time level
!! t^n+1 = t^n + dt, we use the following algorithm:
!!
!!  - let U^(0) = U^n
!!  - for i = 1,...,s compute
!!    U^(i) = sum_{k=0}^{i-1} [ a_ik * U^(k) + dt*b_ik * S(U^(k)) ]
!!  - set U^n+1 = U^(s).
!!
!! The Runge-Kutta method coded in this module follows the above formalism.
!! To reduce memory usage, we assume that
!! the coefficient matrices A = (a_ik) has non-zero entries only
!! in the first column and the last row, and along the diagonal, i.e.,
!!
!!     | a_10     0     0      0    0   ...   0         0     |
!!     | a_20    a_21   0      0    0   ...   0         0     |
!!     | a_30     0    a_32    0    0   ...   0         0     |
!! A = | a_40     0     0     a_43  0   ...   0         0     |  ;
!!     | .        .     .      .    .         0         0     |
!!     | .        0     0      0    0   ...   0         0     |
!!     | a_s-1,0  0     0      0    0   ... a_s-1,s-2   0     |
!!     | a_s0    a_s1  a_s2   a_s3 a_s4 ... a_s,s-2   a_s,s-1 |
!!
!! The matrix B = (b_ik) has non-zero entries
!! only in the last row and along the diagonal, i.e.,
!!
!!     | b_10     0     0      0    0   ...   0         0     |
!!     |  0      b_21   0      0    0   ...   0         0     |
!!     |  0       0    b_32    0    0   ...   0         0     |
!! B = |  0       0     0     b_43  0   ...   0         0     |  .
!!     |  .       .     .      .    .         0         0     |
!!     |  .       0     0      0    0   ...   0         0     |
!!     |  0       0     0      0    0   ... b_s-1,s-2   0     |
!!     | b_s0    b_s1  b_s2   b_s3 b_s4 ... b_s,s-2   b_s,s-1 |
!!
!! This assumption is valid, e.g., for the standard RK4
!! method which has
!!
!!     | 1  0  0  0 |          | 1/2   0    0   0  |
!! A = | 1  0  0  0 |  ,   B = |  0   1/2   0   0  | ,
!!     | 1  0  0  0 |          |  0    0    1   0  |
!!     | 1  0  0  0 |          | 1/6  1/3  1/3 1/6 |
!!
!! and for the SSPRK(5,4) scheme which was proposed by Spiteri and Ruuth (2002,
!! SIAM J. Numer. Anal) and used in Giorgetta et al. (2009).
!!
!! @author Hui Wan, MPI-M
!!
!!
!! @par Revision History
!! Initial version by Hui Wan (2009-08-03)
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
MODULE mo_ha_rungekutta

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: finish
  USE mo_impl_constants,   ONLY: SUCCESS, RK4, SSPRK54
  USE mo_model_domain,     ONLY: t_patch
  USE mo_ext_data_types,   ONLY: t_external_data
  USE mo_intp_data_strc,   ONLY: t_int_state
  USE mo_parallel_config,  ONLY: nproma
  USE mo_run_config,       ONLY: nlev, ltransport
  USE mo_icoham_dyn_types, ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_ha_prog_util,     ONLY: copy_prog_state
  USE m_dyn,               ONLY: dyn_theta
  USE mo_ha_dynamics,      ONLY: dyn_temp
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,      ONLY: get_indices_c, get_indices_e
  USE mo_timer,            ONLY: ltimer, timer_start, timer_stop,&
    & timer_RK_tend, timer_RK_update, timer_step_RK

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_RungeKutta, step_RungeKutta, cleanup_RungeKutta
!   PUBLIC :: RK_update_prog

  !! Module variables

  INTEGER :: nstage

  REAL(wp),ALLOCATABLE :: a_sk(:) !< the last row of matrix A
  REAL(wp),ALLOCATABLE :: b_sk(:) !< the last row of matrix B

  REAL(wp),ALLOCATABLE :: a_i0(:) !< the first column of matrix A

  REAL(wp),ALLOCATABLE :: a_dg(:) !< the diagonal of matrix A
  REAL(wp),ALLOCATABLE :: b_dg(:) !< the diagonal of matrix B

  REAL(wp),ALLOCATABLE :: accm(:) !< flag indicating whether the
                                  !< previous stage exists.
                                  !<  = 0  for stage 0
                                  !<  = 1  for the other stages
  REAL(wp),ALLOCATABLE :: c_sk(:)

  !! Parameters of the standard 4th order Runge-Kutta method

  REAL(wp),PARAMETER :: a10_RK4 = 1._wp
  REAL(wp),PARAMETER :: a20_RK4 = 1._wp, a21_RK4 = 0._wp
  REAL(wp),PARAMETER :: a30_RK4 = 1._wp, a32_RK4 = 0._wp

  REAL(wp),PARAMETER :: a40_RK4 = 1._wp, a41_RK4 = 0._wp
  REAL(wp),PARAMETER :: a42_RK4 = 0._wp, a43_RK4 = 0._wp

  REAL(wp),PARAMETER :: b10_RK4 = 0.5_wp, b21_RK4 = 0.5_wp, b32_RK4 = 1._wp
  REAL(wp),PARAMETER :: b40_RK4 = 1._wp/6._wp
  REAL(wp),PARAMETER :: b41_RK4 = 1._wp/3._wp
  REAL(wp),PARAMETER :: b42_RK4 = 1._wp/3._wp
  REAL(wp),PARAMETER :: b43_RK4 = 1._wp/6._wp

  REAL(wp),PARAMETER :: c40_RK4 = b40_RK4
  REAL(wp),PARAMETER :: c41_RK4 = b41_RK4
  REAL(wp),PARAMETER :: c42_RK4 = b42_RK4
  REAL(wp),PARAMETER :: c43_RK4 = b43_RK4

  !! Parameters of the SSPRK(5,4) method

  REAL(wp),PARAMETER :: a10=1.0_wp
  REAL(wp),PARAMETER :: a20=0.44437049406734_wp , a21=1.0_wp-a20
  REAL(wp),PARAMETER :: a30=0.62010185138540_wp , a32=1.0_wp-a30
  REAL(wp),PARAMETER :: a40=0.17807995410773_wp , a43=1.0_wp-a40

  REAL(wp),PARAMETER :: a50=0.00683325884039_wp , a51=0._wp
  REAL(wp),PARAMETER :: a52=0.51723167208978_wp , a53=0.12759831133288_wp
  REAL(wp),PARAMETER :: a54=1.0_wp-(a50+a52+a53)

  REAL(wp),PARAMETER :: b10=0.39175222700392_wp
  REAL(wp),PARAMETER :: b21=0.36841059262959_wp
  REAL(wp),PARAMETER :: b32=0.25189177424738_wp
  REAL(wp),PARAMETER :: b43=0.54497475021237_wp

  REAL(wp),PARAMETER :: b50=0._wp, b51=0._wp, b52=0._wp
  REAL(wp),PARAMETER :: b53=0.08460416338212_wp
  REAL(wp),PARAMETER :: b54=1.0_wp-((a52*a21+a53*a32*a21+a54*a43*a32*a21)*b10 &
                                   +(    a52+    a53*a32+    a54*a43*a32)*b21 &
                                   +(                a53+        a54*a43)*b32 &
                                   +                                 a54 *b43 &
                                   + b53 )

  REAL(wp),PARAMETER :: c50=b50+ (a51+a52*a21+a53*a32*a21+a54*a43*a32*a21)*b10
  REAL(wp),PARAMETER :: c51=b51+ (    a52    +a53*a32    +a54*a43*a32    )*b21
  REAL(wp),PARAMETER :: c52=b52+ (            a53        +a54*a43        )*b32
  REAL(wp),PARAMETER :: c53=b53+                          a54             *b43
  REAL(wp),PARAMETER :: c54=b54

CONTAINS

  !>
  !! SUBROUTINE init_RungeKutta

  SUBROUTINE init_RungeKutta( p_method )

  INTEGER,INTENT(in) :: p_method   !< the chosen method
  INTEGER            :: ist        !< tmp. variable for checking
                                   !< the memory allocation status

  SELECT CASE(p_method)
  CASE(RK4) ! standard RK4 method

  ! The standard RK4 method is 4-stage

    nstage = 4

  ! allocate memory for the coefficients

    ALLOCATE( a_i0(1:nstage-1)                                      &
            , a_dg(1:nstage-1), b_dg(1:nstage-1)                    &
            , a_sk(0:nstage-1), b_sk(0:nstage-1), c_sk(0:nstage-1)  &
            , accm(0:nstage-1), STAT=ist )

    IF (ist/=SUCCESS) THEN
      CALL finish('init_RungeKutta (std. RK4)', 'allocation failed')
    ENDIF

  ! the coefficients
  ! Note that a_i0(1) is set to zero because a_i0(1) and a_dg(1)
  ! are the same entry.

    a_i0 = (/0._wp,   a20_RK4, a30_RK4/)

    a_dg = (/a10_RK4, a21_RK4, a32_RK4/)
    b_dg = (/b10_RK4, b21_RK4, b32_RK4/)

    a_sk = (/a40_RK4, a41_RK4, a42_RK4, a43_RK4/)
    b_sk = (/b40_RK4, b41_RK4, b42_RK4, b43_RK4/)
    c_sk = (/c40_RK4, c41_RK4, c42_RK4, c43_RK4/)

    accm = (/0._wp,   1._wp,   1._wp,   1._wp/)

  CASE(SSPRK54) ! SSP RK(5,4)

  ! This method has 5 stages

    nstage = 5

  ! allocate memory for the coefficients

    ALLOCATE( a_i0(1:nstage-1)                                      &
            , a_dg(1:nstage-1), b_dg(1:nstage-1)                    &
            , a_sk(0:nstage-1), b_sk(0:nstage-1), c_sk(0:nstage-1)  &
            , accm(0:nstage-1), STAT=ist )

    IF (ist/=SUCCESS) THEN
       CALL finish('init_RungeKutta (SSP)', 'allocation failed')
    ENDIF

  ! the coefficients
  ! Note that a_i0(1) is set to zero because a_i0(1) and a_dg(1)
  ! are the same entry.

    a_i0 = (/0._wp,a20,a30,a40/)

    a_dg = (/a10,a21,a32,a43/)
    b_dg = (/b10,b21,b32,b43/)

    a_sk = (/a50,a51,a52,a53,a54/)
    b_sk = (/b50,b51,b52,b53,b54/)
    c_sk = (/c50,c51,c52,c53,c54/)

    accm = (/0._wp,1._wp,1._wp,1._wp,1._wp/)

  CASE DEFAULT
    CALL finish('init_RungeKutta','unknow choice of method')
  END SELECT

  END SUBROUTINE init_RungeKutta

  !>
  !! SUBROUTINE cleanup_RungeKutta
  !!
  !! Short description:
  !! Release memory used by the Runge-Kutta coefficients.
  !!
  !! @par Revision History
  !! Initial version by Hui Wan (2009-08-03)

  SUBROUTINE cleanup_RungeKutta

  INTEGER :: ist

    DEALLOCATE( a_i0, a_dg, b_dg, a_sk, b_sk, c_sk, accm, STAT=ist )

    IF (ist/=SUCCESS) THEN
       CALL finish('cleanup_RungeKutta', 'Deallocation failed')
    ENDIF

  END SUBROUTINE cleanup_RungeKutta

  !>
  !! SUBROUTINE step_RungeKutta
  !!
  !! @par Revision History
  !! Initial version by Hui Wan (2009-08-04)
  !!
  SUBROUTINE step_RungeKutta( pdtime, ltheta_dyn      &!< in
                            , curr_patch, p_int_state &!< in
                            , p_now, p_ext_data       &!< in
                            , p_stg, p_new            &!< in,tmp,inout
                            , p_diag, p_tend_dyn      &!< inout
                            )

  !! input and output variables

   REAL(wp),INTENT(IN) :: pdtime     !< time step in seconds
   LOGICAL, INTENT(IN) :: ltheta_dyn

   TYPE(t_patch), TARGET, INTENT(INOUT) :: curr_patch
   TYPE(t_int_state), TARGET, INTENT(IN)  :: p_int_state

   TYPE(t_hydro_atm_prog),INTENT(IN)    :: p_now
   TYPE(t_external_data), INTENT(INOUT) :: p_ext_data !< external data
   TYPE(t_hydro_atm_prog),INTENT(INOUT) :: p_stg
   TYPE(t_hydro_atm_prog),INTENT(INOUT) :: p_new

   TYPE(t_hydro_atm_diag),INTENT(INOUT) :: p_diag
   TYPE(t_hydro_atm_prog),INTENT(INOUT) :: p_tend_dyn !adiabatic tendencies

   REAL(wp) :: z_mflxh_me (nproma,nlev,curr_patch%nblks_e)
   REAL(wp) :: z_mflxv_ic (nproma,nlev,curr_patch%nblks_c)

   INTEGER :: ii, kk               !< stage indices


    IF (ltimer) CALL timer_start(timer_step_RK)
     !--------------------------------------------------------------------------
     ! If tracer transport is on, initialize mass fluxes

     IF (ltransport) THEN
!$OMP PARALLEL WORKSHARE
       z_mflxh_me(:,:,:) = 0._wp
       z_mflxv_ic(:,:,:) = 0._wp
!$OMP END PARALLEL WORKSHARE
     ENDIF

     ! U^(0) = U^n
     CALL copy_prog_state( p_now, p_stg, ltheta_dyn, .FALSE. )

     ! For stage ii, let kk = ii -1, do:
     !  1) compute tendency S(U^(kk));
     !  2) accumulate the contribution of U^(kk) and S(U^(kk)) to U^n+1;
     !  3) calculate U^(ii) from U^(0), U^(kk), and S(U^(kk)).
     ! Because a) the prognostic state vector and the tendency vector have
     ! different components, and b) we have to treat boundaries separately
     ! in the case of grid refinement,
     ! steps 2) and 3) can not be done by simply using array syntax.
     ! Therefore we use the subroutine RK_update_prog
     !----------------------
     ! stages 1 to nstage-1
     !----------------------
     StageLoop: DO ii = 1,nstage-1

       kk = ii-1

       ! compute tendencies
       IF (ltimer) CALL timer_start(timer_RK_tend)
       IF (ltheta_dyn) THEN
         CALL dyn_theta( curr_patch, p_int_state, p_ext_data, p_stg,  & ! in
                         p_diag, p_tend_dyn )                           ! inout
       ELSE
         CALL dyn_temp(  curr_patch, p_int_state, p_stg, p_ext_data,  & ! in
                         p_diag, p_tend_dyn )                           ! inout
       ENDIF
       IF (ltimer) CALL timer_stop(timer_RK_tend)

      ! p_new = p_new + a_s,i-1*U^(i-1) + b_s,i-1*dt*S(U^(i-1))

       CALL RK_update_prog( p_new                        & ! inout
                          , accm(kk),        p_new       & ! in, in
                          , a_sk(kk),        p_stg       & ! in, in
                          , b_sk(kk)*pdtime, p_tend_dyn  & ! in, in
                          , curr_patch,      ltheta_dyn  ) ! in

      ! U^(i) = a_i0*U^(0) + a_i,i-1*U^(i-1) + b_i,i-1*dt*S(U^(i-1))

       CALL RK_update_prog( p_stg                       & ! inout
                          , a_i0(ii),        p_now      & ! in, in
                          , a_dg(ii),        p_stg      & ! in, in
                          , b_dg(ii)*pdtime, p_tend_dyn & ! in, in
                          , curr_patch,      ltheta_dyn ) ! in

      ! If tracer transport is on, accumulate mass flux

      IF (ltransport) THEN
!$OMP PARALLEL WORKSHARE
        z_mflxh_me(:,1:nlev,:) = z_mflxh_me(:,1:nlev,:) + c_sk(kk)*p_diag%mass_flux_e(:,1:nlev,:)
        z_mflxv_ic(:,2:nlev,:) = z_mflxv_ic(:,2:nlev,:) + c_sk(kk)*p_diag%weta(:,2:nlev,:)
!$OMP END PARALLEL WORKSHARE
      ENDIF

     ENDDO StageLoop

     !-----------------
     ! The last stage
     !-----------------
     ! Similar to the previous stages but only do steps 1) and 2) listed above.

       kk = nstage -1

       ! compute tendencies
       IF (ltimer) CALL timer_start(timer_RK_tend)
       IF (ltheta_dyn) THEN
         CALL dyn_theta( curr_patch, p_int_state, p_ext_data, p_stg, & ! in
                         p_diag, p_tend_dyn )                          ! out
       ELSE
         CALL dyn_temp(  curr_patch, p_int_state, p_stg, p_ext_data, & ! in
                         p_diag, p_tend_dyn )                          ! out
       ENDIF
       IF (ltimer) CALL timer_stop(timer_RK_tend)

      ! p_new = p_new + a_s,s-1*U^(s-1) + b_s,s-1*dt*S(U^(s-1))

       CALL RK_update_prog( p_new                        & ! inout
                          , accm(kk),        p_new       & ! in, in
                          , a_sk(kk),        p_stg       & ! in, in
                          , b_sk(kk)*pdtime, p_tend_dyn  & ! in, in
                          , curr_patch,      ltheta_dyn  ) ! in

      ! If tracer transport is on, accumulate mass flux

      IF (ltransport) THEN
!$OMP PARALLEL WORKSHARE
        p_diag%mass_flux_e(:,1:nlev,:) =   z_mflxh_me(:,1:nlev,:)                  &
                                       & + c_sk(kk)*p_diag%mass_flux_e(:,1:nlev,:)
        p_diag%weta       (:,2:nlev,:) =   z_mflxv_ic(:,2:nlev,:)                  &
                                       & + c_sk(kk)*p_diag%weta(:,2:nlev,:)
!$OMP END PARALLEL WORKSHARE
      ENDIF

    IF (ltimer) CALL timer_stop(timer_step_RK)
    
  END SUBROUTINE step_RungeKutta

  !>
  !! SUBROUTINE RK_update_prog
  !!
  !! @par Revision History
  !! Initial version by Hui Wan (2009-08-05)
  !!
  SUBROUTINE RK_update_prog( &
             p_new       &!< inout: new values of the prog. variables
           , pc1         &!< in: coefficient of p_old1
           , p_old1      &!< in: prog. variables at old time level 1
           , pc2         &!< in: coefficient of p_old2
           , p_old2      &!< in: prog. variables at old time level 2
           , pdtime      &!< in: time step
           , p_tend      &!< in: tendency of the prog. variables
           , curr_patch  &!< in: patch info.
           , ltheta_dyn  &!< in
           )

  !! input and output variables
  !! The intent property of p_new is "inout" instead of "out" because
  !! the tracer concentrations should be kept intact.

   TYPE(t_hydro_atm_prog),INTENT(INOUT) :: p_new
   TYPE(t_hydro_atm_prog),INTENT(IN)    :: p_old1
   TYPE(t_hydro_atm_prog),INTENT(IN)    :: p_old2
   TYPE(t_hydro_atm_prog),INTENT(IN)    :: p_tend
   REAL(wp),              INTENT(IN)    :: pc1
   REAL(wp),              INTENT(IN)    :: pc2
   REAL(wp),              INTENT(IN)    :: pdtime
   TYPE(t_patch),         INTENT(IN)    :: curr_patch
   LOGICAL,               INTENT(IN)    :: ltheta_dyn

  !! local variables

   INTEGER :: i_startblk, i_startidx, i_endidx, i_endblk
   INTEGER :: jb, jc, jk, je

   IF (ltimer) CALL timer_start(timer_RK_update)
   !--------------------------------------------------------
   ! update cell-based variables in the interior
   !--------------------------------------------------------
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)
   i_startblk = curr_patch%cells%start_blk(grf_bdywidth_c+1,1)
   i_endblk   = curr_patch%nblks_c
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk, i_endblk

     CALL get_indices_c(curr_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, grf_bdywidth_c+1)

      DO jc = i_startidx, i_endidx
         p_new%pres_sfc(jc,jb) = pc1*p_old1%pres_sfc(jc,jb)  &
                               + pc2*p_old2%pres_sfc(jc,jb)  &
                            + pdtime*p_tend%pres_sfc(jc,jb)
      ENDDO

      IF (ltheta_dyn) THEN
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            p_new%theta(jc,jk,jb) = pc1*p_old1%theta(jc,jk,jb)  &
                                  + pc2*p_old2%theta(jc,jk,jb)  &
                               + pdtime*p_tend%temp(jc,jk,jb)
          ENDDO
        ENDDO
      ELSE
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            p_new%temp(jc,jk,jb) = pc1*p_old1%temp(jc,jk,jb)  &
                                 + pc2*p_old2%temp(jc,jk,jb)  &
                              + pdtime*p_tend%temp(jc,jk,jb)
          ENDDO
        ENDDO
      ENDIF

   ENDDO
!$OMP END DO

   !--------------------------------------------------------
   ! update edge-based variable in the interior
   !--------------------------------------------------------

   i_startblk = curr_patch%edges%start_blk(grf_bdywidth_e+1,1)
   i_endblk   = curr_patch%nblks_e
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk, i_endblk

     CALL get_indices_e(curr_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, grf_bdywidth_e+1)

      DO jk = 1, nlev
         DO je = i_startidx, i_endidx
            p_new%vn(je,jk,jb) = pc1*p_old1%vn(je,jk,jb) &
                               + pc2*p_old2%vn(je,jk,jb) &
                            + pdtime*p_tend%vn(je,jk,jb)
         ENDDO ! edge loop
      ENDDO ! level loop

   ENDDO ! block loop
!$OMP END DO

   !--------------------------------------------------------
   ! update cell-based variables at boundaries
   !--------------------------------------------------------

   i_startblk = curr_patch%cells%start_blk(1,1)
   i_endblk   = curr_patch%cells%end_blk(grf_bdywidth_c,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk, i_endblk

     CALL get_indices_c(curr_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, 1, grf_bdywidth_c)

      DO jc = i_startidx, i_endidx
         p_new%pres_sfc(jc,jb) = pc1*p_old1%pres_sfc(jc,jb)  &
                               + pc2*p_old2%pres_sfc(jc,jb)  &
                            + pdtime*p_tend%pres_sfc(jc,jb)
      ENDDO

      IF (ltheta_dyn) THEN
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            p_new%theta(jc,jk,jb) = pc1*p_old1%theta(jc,jk,jb)  &
                                  + pc2*p_old2%theta(jc,jk,jb)  &
                                + pdtime*p_tend%temp(jc,jk,jb)
          ENDDO
        ENDDO
      ELSE
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            p_new%temp(jc,jk,jb) = pc1*p_old1%temp(jc,jk,jb)  &
                                 + pc2*p_old2%temp(jc,jk,jb)  &
                              + pdtime*p_tend%temp(jc,jk,jb)
          ENDDO
        ENDDO
      ENDIF

   ENDDO
!$OMP END DO

   !--------------------------------------------------------
   ! update edge-based variable at boundaries
   !--------------------------------------------------------

   i_startblk = curr_patch%edges%start_blk(1,1)
   i_endblk   = curr_patch%edges%end_blk(grf_bdywidth_e,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk, i_endblk

     CALL get_indices_e(curr_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, 1, grf_bdywidth_e)

      DO jk = 1, nlev
         DO je = i_startidx, i_endidx
            p_new%vn(je,jk,jb) = pc1*p_old1%vn(je,jk,jb) &
                               + pc2*p_old2%vn(je,jk,jb) &
                            + pdtime*p_tend%vn(je,jk,jb)
         ENDDO
      ENDDO

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   IF (ltimer) CALL timer_stop(timer_RK_update)
   
   END SUBROUTINE RK_update_prog

END MODULE mo_ha_RungeKutta

