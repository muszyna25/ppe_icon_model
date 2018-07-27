!>
!!  Semi-implicit correction for the leapfrog time stepping scheme.
!!
!!  There are two options:
!! <ol>
!! <li> directly solve a 3D Helmholtz equation;
!! <li> decompose the 3D equation => solve 2D equations of the fastest modes
!!     => re-couple the solutions (coded after GME of DWD).
!! </ol>
!!
!! @par Revision History
!!  Original version (solving 2D equations of divergence)
!!    by Hui Wan (MPI-M, 2007-08-08).
!!    Subroutines <i>conteq</i> and <i>pgrad</i> were adapted from ECHAM5.3.01.
!!    The solver used was SOR.
!!  Modifications by Hui Wan, MPI-M (2007-12-15)
!!  - implemented the GMRES solver.
!!  Modifications by Hui Wan, MPI-M (2008-01-25)
!!  - Recoded after GME to solve 2D equations of the second temproal derivative
!!    of divergence.
!!  Modifications by Hui Wan, MPI-M (2008-01-31)
!!  - added the option of directly solving the 3D equation of
!!    the second temproal derivative of divergence.
!!  Modifications by Hui Wan, MPI-M (2008-02-06)
!!  - merged the two options into one module so that they can be chosen
!!    through a namelist variable.
!!  Modifications by Hui Wan, MPI-M (2008-02-13)
!!  - sorting of the eigenvalues added.
!!  Modifications by Hui Wan, MPI-M (2008-04-05)
!!  - ldivavg, divavg_c0 and divavg_cj renamed ldiv_avg, div_avg_c0
!!    and div_avg_c2, respectively.
!!  - Offcentering parameters asi and bsi removed from namelist.
!!    They are now private module variables. The values are calculated from
!!    namelist variables si_offctr_temp and si_offctr_vn, respectively.
!!  - relative tolerance rtol renamed si_rtol.
!!  Modifications by Hui Wan, MPI-M (2008-05-29)
!!  - checked the combination of ldiv_avg = .true. and lsi_3d = .false.
!!  Restructuring by Almut Gassman, MPI-M (2008-09-18)
!!  Modifications by Marco Giorgetta, MPI-M (2009-02-14)
!!  - Instead of passing the identity operator as preconditioner to "gmres"
!!    no argument is passed. The dummy argument for the preconditioner in
!!    "gmres" is now an optional argument.
!!  - A non-identity preconditioner can be passed to "gmres" as last
!!    argument.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_si_correction
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
! !USES
  USE mo_kind,               ONLY: wp, dp
  USE mo_exception,          ONLY: message, finish
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_physical_constants, ONLY: rcpd, rd, grav
  USE mo_math_gradients,     ONLY: grad_fd_norm
  USE mo_math_divrot,        ONLY: div, div_avg
  USE mo_math_laplace,       ONLY: nabla2_scalar, nabla2_scalar_avg
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_model_domain,       ONLY: t_patch
  USE mo_grid_config,        ONLY: l_limited_area
  USE mo_dynamics_config,    ONLY: idiv_method, sw_ref_height
  USE mo_parallel_config,    ONLY: nproma
  USE mo_run_config,         ONLY: msg_level, nlev, nlevm1, nlevp1
  USE mo_icoham_dyn_types,   ONLY: t_hydro_atm_prog
  USE mo_vertical_coord_table,ONLY: delpr, rdelpr
  USE mo_eta_coord_diag,     ONLY: half_level_pressure, auxhyb
  USE mo_gmres,               ONLY: gmres
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_sync,               ONLY: SYNC_C, SYNC_E, sync_patch_array
  
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop,&
    & timer_si_correction

  IMPLICIT NONE

  PRIVATE

  REAL(wp),PARAMETER :: tr = 300._wp      ! reference temperature
  REAL(wp),PARAMETER :: rdtr = rd*tr      ! rd*(reference temperature).
  LOGICAL,PARAMETER  :: loperm = .false.  ! used by subroutine conteq.

! !ARRAYS

! GZ: module arrays with runtime-determined shape are not allowed
  REAL(wp), ALLOCATABLE :: &
            & ralphr(:),   & ! *constant array for use by pgrad.
            & rlnmar(:),   & ! *constant array for use by pgrad.
            & aktlrd(:),   & ! *constant array for use by conteq.
            & altrcp(:)      ! *constant array for use by conteq.

  REAL(wp), ALLOCATABLE :: bb(:,:) ! transpose of the gravity wave matrix

  REAL(dp), ALLOCATABLE :: eigenval(:) ! eigenvalues of bb-tranpose

  REAL(dp), ALLOCATABLE ::    &
            & eigenvec(:,:),  & ! eigenvectors of bb-tranpose,
            & reigenvec(:,:)    ! reverse of the eignenvector matrix

  INTEGER  :: nmodes  ! number of modes whose phase speed is larger than si_cmin

  REAL(wp) :: asi_n, asi_o ! implicit weights on new and old time step

  REAL(wp) :: apr != 80000._wp:   reference surface pressure
  REAL(wp) :: rpr != 1._wp/apr:   reciprocal of reference surface pressure.

  PUBLIC  :: init_si_params, si_correction

  PRIVATE :: conteq, pgrad, add_si_correction_2d, add_si_correction_3d, &
             sort_eigen

  CONTAINS

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
!

  !>
  !! Initializes constants for the vertical part of the semi-implicit scheme.
  !!
  !! Compute constants used in the implementation of the
  !! vertical part of the semi-implicit time scheme.
  !!
  !! Input is from module *mo_hyb_params* and from this module.
  !! Output is array *bb*, its eigenvalues and eigenvectors in this module.
  !!
  !! @par Revision History
  !!  A. J. Simmons, ECMWF, November 1981, original source
  !!  L. Kornblueh, MPI, May 1998, f90 rewrite
  !!  U. Schulzweida, MPI, May 1998, f90 rewrite
  !!  I. Kirchner, MPI, December 2000, time control
  !! @par
  !!  H. Wan, MPI, May 2007, adapted for ICON hydrostatic dynamical core:
  !!   - moved the declaration and allocation of some arrays here
  !!     from mo_hyb_params;
  !!   - added calculation of the eigenvectors of structure matrix.
  !!   - changed the name from inhysi to init_si_params.
  !!
  SUBROUTINE init_si_params( lsi_3d, si_offctr, si_cmin, lshallow_water )

  LOGICAL, INTENT(IN) :: lsi_3d
  REAL(wp),INTENT(IN) :: si_offctr, si_cmin
  LOGICAL, INTENT(IN) :: lshallow_water
  INTEGER  :: jk, ist

! !for building up the structure matrix
  REAL(wp) :: z_div(nlev,nlev), z_dtemp (nlev,nlev)
  REAL(wp) :: z_dlnps(nlev), z_rlnpr(nlev)
  REAL(wp) :: z_ph(nlevp1), z_lnph(nlevp1)
  REAL(wp) :: z_pra(1)

! !for finding the eigenvectors and the inverse matrix
  REAL(dp) :: z_wrk2(nlev,nlev), &
              z_wrk3(nlev,nlev), &
              z_vl  (nlev,nlev), &
              z_wrk1(6*nlev),    & ! temporary arrays
              z_ei  (nlev) ! imag. part of eigenvalues. should be all zeros
  INTEGER  :: ipivot(nlev)

  REAL(wp) :: z_maxi       ! for checking whether the eigenvalues are real

  CHARACTER(len=MAX_CHAR_LENGTH) :: routine = 'mo_si_correction:init_si_params'
  CHARACTER(len=MAX_CHAR_LENGTH) :: string
!-------------------------------------------------------------------------

  asi_n = si_offctr*2._wp
  asi_o = 2.0_wp-asi_n

!-- 0. Allocate memory for time-independent parameters

  ALLOCATE (ralphr(nlev), STAT=ist)
  IF(ist/=SUCCESS)THEN
     CALL finish (TRIM(routine), ' allocation of ralphr failed')
  ENDIF

  ALLOCATE (rlnmar(nlev), STAT=ist)
  IF(ist/=SUCCESS)THEN
     CALL finish (TRIM(routine), ' allocation of rlnmar failed')
  ENDIF

  ALLOCATE (aktlrd(nlev), STAT=ist)
  IF(ist/=SUCCESS)THEN
     CALL finish (TRIM(routine), ' allocation of aktlrd failed')
  ENDIF

  ALLOCATE (altrcp(nlev), STAT=ist)
  IF(ist/=SUCCESS)THEN
     CALL finish (TRIM(routine), ' allocation of altrcp failed')
  ENDIF

  ALLOCATE (eigenval(nlev), STAT=ist)
  IF(ist/=SUCCESS)THEN
     CALL finish (TRIM(routine), ' allocation of eigenval failed')
  ENDIF

  ALLOCATE (bb(nlev,nlev), STAT=ist)
  IF(ist/=SUCCESS)THEN
     CALL finish (TRIM(routine), ' allocation of bb failed')
  ENDIF

  ALLOCATE (eigenvec(nlev,nlev), STAT=ist)
  IF(ist/=SUCCESS)THEN
     CALL finish (TRIM(routine), ' allocation of eigenvec failed')
  ENDIF

  ALLOCATE (reigenvec(nlev,nlev), STAT=ist)
  IF(ist/=SUCCESS)THEN
     CALL finish (TRIM(routine), ' allocation of reigenvec failed')
  ENDIF

IF (.NOT. lshallow_water) THEN

!-- 1. Set local values

  apr = 80000.0_wp
  rpr = 1.0_wp/apr
  z_pra(1) = apr

!-- 2. Calculate delpr, ralphr and rlnmar

  CALL half_level_pressure(z_pra,1,1,z_ph)
  CALL auxhyb(z_ph,1,1,delpr,rdelpr,z_lnph,z_rlnpr,ralphr)

  rlnmar(1) = 0._wp
  DO jk = 2, nlev
    rlnmar(jk) = z_rlnpr(jk) - ralphr(jk)
  END DO

!-- 3. Calculate aktlrd and altrcp

  DO jk = 1, nlev
    aktlrd(jk) = -rcpd*tr*z_rlnpr(jk)*rdelpr(jk)
    altrcp(jk) = -rcpd*tr*ralphr(jk)
  END DO

ENDIF

!-- 4. Calculate gravity-wave matrix bb

  IF (.NOT.lsi_3d) THEN

   IF (lshallow_water) THEN

     ! Shallow water mode
     !===================

     eigenval(1)    = sw_ref_height * grav
     eigenvec(1,1)  = 1.0_wp
     reigenvec(1,1) = 1.0_wp
     apr            = 1.0_wp
     rpr            = 1.0_wp
     nmodes         = 1

   ELSE

     !hydrostatic mode
     !===================

     z_div(:,:) = 0._wp

     DO jk = 1, nlev
        z_div(jk,jk) = -1._wp
     END DO

     CALL conteq( z_div,nlev,nlev,loperm,    & ! input
                  z_dtemp,z_dlnps )            ! output

     CALL pgrad ( z_dtemp,z_dlnps,nlev,nlev, & ! input
                  bb )                         ! output

!-- 6. The eigenvectors and the inverse

!-- 6.1 The eigenvalues and eigenvectors

    z_wrk2(:,:) = TRANSPOSE( bb(:,:) )

    CALL dgeev (&       ! LAPACK eigenvalue solver A*X = lamdba*X
         'N',&          ! calculate no left eigenvectors
         'V',&          ! calculate right eigenvectors
         nlev,&         ! order of matrix A
         z_wrk2,&       ! inp: structure matrix A, out: overwritten
         nlev,&         ! leading dimension of A
         eigenval,&     ! real part of eigenvalues
         z_ei,&         ! imaginary part of eigenvalues
         z_vl,&         ! contains nothing
         nlev,&         ! leading dimension of left eigenvectors
         eigenvec,&     ! right eigenvectors in columns
         nlev,&         ! leading dimension of right eigenvectors
         z_wrk1,&       ! work space
         6*nlev,&       ! dimension of workspace
         ist)           ! error return code

  !---------------------------------------
  ! quality control
  !---------------------------------------

    IF (ist /= 0) THEN
      CALL finish (TRIM(routine),'Calculation of eigenvectors/values failed.')
    ENDIF

    z_maxi = MAXVAL(ABS(z_ei(1:nlev)))  !eigenvalue with largest imaginary part
    IF ( z_maxi > 0._wp ) THEN
         CALL finish (TRIM(routine),&
         'Found complex eigenvalues. Failure in vertical structure matrix')
    ENDIF

  !--------------------------------------------------------------------
  ! sort the eigenvectors according to the corresponding eigenvalues.
  !--------------------------------------------------------------------

    CALL sort_eigen( eigenval, eigenvec, nlev )

  !--------------------------------------------------------------------
  ! determine how many modes to solve
  !--------------------------------------------------------------------

    nmodes = 0
    DO jk=1,nlev
       WRITE(string,'(a,i4,a,f10.5)')   &
                    'mode',jk,' phase speed =',SQRT(eigenval(jk))
       CALL message(TRIM(routine),TRIM(string))
       IF ( SQRT(eigenval(jk)) > si_cmin ) nmodes = nmodes + 1
    ENDDO

    WRITE(string,'(a,f6.2,a,i4)')   &
                 'Number of modes with phase speed > ',si_cmin,' is ', nmodes
    CALL message(TRIM(routine),TRIM(string))

!-- 6.2 The inverse of the eigenvector matrix

  !--------------------------------------------
  ! the unit matrix
  !--------------------------------------------

   z_wrk3 = 0._wp
   DO jk=1,nlev
      z_wrk3(jk,jk) = 1._wp
   END DO

  !----------------------------------------------
  ! the working array containing the eigenvectors
  !----------------------------------------------

   z_wrk2(:,:) = eigenvec(:,:)

  !--------------------------------------------
  ! calculate the inverse
  !--------------------------------------------

   CALL dgesv (& ! BLAS linear equation solver
        nlev,&   ! number of linear equations NEQ [int]
        nlev,&   ! number of right hand sides NRHS [int]
        z_wrk2,& ! in: coefficients matrix A(LDA,NEQ),
                 ! out: l-u-factorization [real]
        nlev,&   ! leading dimension LDA [int]
        ipivot,& ! out: permutation matrix IPIV(NEQ) [int]
        z_wrk3,& ! in: right hand side X(LDX,RHS), out: solution matrix [real]
        nlev,&   ! leading dimension LDX [int]
        ist)     ! error code [int]

  !---------------------------------------
  ! quality control
  !---------------------------------------

    IF (ist /= 0) THEN
      CALL finish (TRIM(routine), &
           'Inverting eigenvector matrix, failure in dgesv')
    ENDIF

  !---------------------------------------
  ! save the result
  !---------------------------------------

    ! Transpose eigenvector matrices for simplification of runtime matrix operations
    reigenvec(:,:) = TRANSPOSE(z_wrk3(:,:))    ! inverse of the eigenvector matrix
    z_wrk2(:,:)    = eigenvec(:,:)
    eigenvec(:,:)  = TRANSPOSE(z_wrk2(:,:))

   ENDIF ! hydrostatic mode
  ENDIF !only if choose to solve decomposed 2D equations


  CALL message(TRIM(routine), &
       'initialization for the semi-implicit correction finished')

  END SUBROUTINE init_si_params
!-------------------------------------------------------------------------
!
!

  !>
  !!  Sort the eigenvalues the largest to the smallest.
  !!
  !!  Rearrange the eigenvectors accordingly.
  !!
  !! @par Revision History
  !!  Original version  by Hui Wan, MPI-M (2008-02-12).
  !!
  SUBROUTINE sort_eigen( p_eigenval, p_eigenvec, k_dim )


   INTEGER,INTENT(IN) :: k_dim

! !INPUT AND OUTPUT:

   REAL(wp), INTENT(INOUT) :: p_eigenval( k_dim )
   REAL(wp), INTENT(INOUT) :: p_eigenvec( k_dim, k_dim )

! !WORKING VARIABLES

   INTEGER  :: jk         !for the loop over all modes
   INTEGER  :: imaxa(1)   !subscript of the max. value, array
   INTEGER  :: imax       !subscript of the max. value

   REAL(wp) :: z_vec( k_dim, k_dim )
   REAL(wp) :: z_sca( k_dim )
!-------------------------------------------------------------------------
!
   IF (MINVAL(p_eigenval(:)) < 0.0_wp) &
      CALL finish('init_si_params','found negative eigenvalue!')
!
!  copy the eigenvalues and eigenvector matrix to temporary array
!
   z_sca(:)   = p_eigenval(:)
   z_vec(:,:) = p_eigenvec(:,:)
!
!  loop over all modes
!
   DO jk=1,k_dim
!
!     find location (i.e., subscript) of the largest eigenvalue
!
      imaxa = MAXLOC(z_sca)  ! looks awkward but assigning imax = MAXLOC(z_sca)
      imax  = imaxa(1)       ! causes a compiler error
!
!     re-arrange the eigenvalue/vector array
!
      p_eigenval(jk)   = z_sca(imax)
      p_eigenvec(:,jk) = z_vec(:,imax)
!
!     mask the already sorted value
!
      z_sca(imax) = -999._wp
   ENDDO

  END SUBROUTINE sort_eigen
!-------------------------------------------------------------------------
!
!

  !>
  !!
  !! @par Revision History
  !!  Original version  by Hui Wan, MPI-M (2007-08-11).
  !!  Code restructuring by Almut Gassmann, MPI-M (2008-09-18)
  !!
  SUBROUTINE si_correction( lsi_3d, si_coeff, si_rtol, &! in
                            lshallow_water, p_dtime,   &! in
                            pt_patch,                  &! in
                            pt_int_state,              &! in
                            pt_prog_old,               &! in
                            pt_prog_now,               &! in
                            pt_prog_new                )! in


    LOGICAL, INTENT(IN) :: lsi_3d
    REAL(wp),INTENT(IN) :: si_coeff, si_rtol
    LOGICAL, INTENT(IN) :: lshallow_water
    REAL(wp),INTENT(IN) :: p_dtime      !< time step in seconds
    TYPE(t_patch), TARGET, INTENT(INOUT) :: pt_patch
    TYPE(t_int_state),INTENT(IN) :: pt_int_state !< horizontal interpolation coeff.

    TYPE(t_hydro_atm_prog) :: pt_prog_old
    TYPE(t_hydro_atm_prog) :: pt_prog_now
    TYPE(t_hydro_atm_prog) :: pt_prog_new

    IF (ltimer) call timer_start(timer_si_correction)
   
    IF (lsi_3d) THEN
      CALL add_si_correction_3d( si_coeff, si_rtol, p_dtime,            &
                                 pt_patch, pt_int_state,                &
                                 pt_prog_old, pt_prog_now, pt_prog_new  )
    ELSE
      CALL add_si_correction_2d( si_coeff, si_rtol, lshallow_water, p_dtime,&
                                 pt_patch, pt_int_state,                &
                                 pt_prog_old, pt_prog_now, pt_prog_new  )
    ENDIF
    
    IF (ltimer) call timer_stop(timer_si_correction)

  END SUBROUTINE si_correction

  !>
  !!
  SUBROUTINE add_si_correction_2d( si_coeff, si_rtol, lshallow_water,&
                                   p_dtime, pt_patch, pt_int_state,  &
                                   pt_old, pt_now, pt_new)

! !DESCRIPTION
!  This subroutine added the semi-implicit correction to the
!  estimated values of V, T and ps at time n+1 given by the
!  leapfrog scheme.
!  After the GME model, the central task is to solve a set of decoupled 2D
!  Helmholtz equations of ** the second temporal derivative ** of divergence.
!
! !REVISION HISTORY
!  Original version by Hui Wan (MPI-M, 2008-01)
!  Code resturcturing by Almut Gassmann, MPI-M, (2008-09-18)
!

  !0!CHARACTER(len=*), PARAMETER :: routine = 'mo_si_correction:add_si_correction_2d'

  REAL(wp),INTENT(IN) :: si_coeff, si_rtol
  LOGICAL, INTENT(IN) :: lshallow_water
  REAL(wp),INTENT(IN) :: p_dtime      ! time step in seconds
  TYPE(t_patch), TARGET, INTENT(INOUT) :: pt_patch   ! single patch
  TYPE(t_int_state),INTENT(IN) :: pt_int_state ! single interpolation state

  TYPE(t_hydro_atm_prog),INTENT(IN) :: pt_old  ! prognostic variables at step n-1
  TYPE(t_hydro_atm_prog),INTENT(IN) :: pt_now  ! prognostic variables at step n

  TYPE(t_hydro_atm_prog),INTENT(INOUT):: pt_new ! input: the estimated value at
                                                ! step n+1 given by leapfrog scheme.
                                                ! output: with the semi-implicit
                                                ! correction added.

#if (defined (__SUNPRO_F95) || defined(__SX__)) && !defined (NOMPI)
  INTEGER  :: jc
#endif

  INTEGER  :: jk, jb, jm
  INTEGER  :: nblks_e, nblks_c, npromz_e, npromz_c, nlen

  REAL(wp) :: z_dt     ! scaled time step in seconds
  REAL(wp) :: z_dtsq   !

  REAL(wp), DIMENSION (nproma,nlev,pt_patch%nblks_c) :: &
            & z_dttdiv,      & ! temporal Laplacian of divergence
            & z_dttpot,      & ! generalized potential
            & z_rhs,         & ! RHS
            & z_rhs_dcpl,    & ! decoupled RHS
            & z_help_c         ! temporary array

  REAL(wp), DIMENSION (nproma,nlev,pt_patch%nblks_e) :: &
            & z_help_e         ! temporary array

  REAL(wp) :: z_dtemp (nproma,nlev)
  REAL(wp) :: z_dlnps (nproma)

  INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx

! for the GMRES solver

  INTEGER, PARAMETER :: nmax_iter = 100
  LOGICAL  :: lmaxiter
  INTEGER  :: niter
  REAL(wp) :: z_residual(nmax_iter)
  REAL(wp) :: z_wrk_cv(nproma,nlev), z_wrk_vc(nproma,nlev)

  CHARACTER(len=MAX_CHAR_LENGTH) :: string

!-----------------------------------------------------------------------

   ! weighted time step (not to mix up with si_offctr)
   z_dt = p_dtime*si_coeff

   nblks_c   = pt_patch%nblks_c
   npromz_c  = pt_patch%npromz_c
   nblks_e   = pt_patch%nblks_e
   npromz_e  = pt_patch%npromz_e

!-----------------------------------------------------------------------
! 1. Calculate the right-hand-side of the 3D Helmholtz equation
!-----------------------------------------------------------------------
! 1.1 second temporal derivative of divergence
!-----------------------------------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1,nblks_e

      IF (jb /= nblks_e) THEN
         nlen = nproma
      ELSE
         nlen = npromz_e
      ENDIF

      z_help_e(1:nlen,1:nlev,jb) = asi_n*pt_new%vn(1:nlen,1:nlev,jb) &
                                 - 2._wp*pt_now%vn(1:nlen,1:nlev,jb) &
                                 + asi_o*pt_old%vn(1:nlen,1:nlev,jb)

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   SELECT CASE(idiv_method)
   CASE (1)

      CALL div(z_help_e, pt_patch, pt_int_state,  &  ! input
               z_dttdiv )  ! output

   CASE (2,3)

      CALL div_avg(z_help_e, pt_patch, pt_int_state, pt_int_state%c_bln_avg,  &
                   z_dttdiv )

   END SELECT

!-----------------------------------------------------------------------
! 1.2 compute second temporal derivative of the generalized potential
!     from temperature and surface pressure.
!-----------------------------------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, nlen, z_dtemp, z_dlnps) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = 1,nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      z_dlnps(1:nlen) = rpr*( - 2._wp*pt_now%pres_sfc(1:nlen,jb) &
                              + asi_n*pt_new%pres_sfc(1:nlen,jb) &
                              + asi_o*pt_old%pres_sfc(1:nlen,jb))

      IF (.NOT. lshallow_water) THEN
         z_dtemp(1:nlen,:) = - 2._wp*pt_now%temp(1:nlen,:,jb) &
                             + asi_n*pt_new%temp(1:nlen,:,jb) &
                             + asi_o*pt_old%temp(1:nlen,:,jb)

         CALL pgrad ( z_dtemp, z_dlnps, nproma, nlen,  & ! input
                      z_dttpot(1:nproma,1:nlev,jb) )     ! output
      ELSE

        DO jk = 1, nlev
          z_dttpot(1:nlen,jk,jb) = z_dlnps(1:nlen) * grav
        ENDDO

      ENDIF

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
!--------------------------------------------------------------
! 1.3 Laplacian of the second temporal derivative
!     of the generalized potential
!--------------------------------------------------------------

   SELECT CASE (idiv_method)
   CASE (1)

     CALL nabla2_scalar( z_dttpot, pt_patch, pt_int_state, z_help_c )

   CASE (2,3)

     CALL nabla2_scalar_avg( z_dttpot, pt_patch, pt_int_state, &
                             pt_int_state%c_bln_avg, z_help_c )

   END SELECT

!$OMP PARALLEL PRIVATE (i_startblk,i_endblk,i_startidx,i_endidx)

   IF (l_limited_area .OR. pt_patch%id > 1) THEN
     ! In case of nested domains, the right-hand-side is set to zero in the
     ! two outermost cell rows
     i_startblk = pt_patch%cells%start_blk(1,1)
     i_endblk   = pt_patch%cells%end_blk(2,1)

     DO jb = i_startblk, i_endblk

       IF (jb == i_startblk) THEN
         i_startidx = pt_patch%cells%start_idx(1,1)
         i_endidx   = nproma
         IF (jb == i_endblk) i_endidx = pt_patch%cells%end_idx(2,1)
       ELSE IF (jb == i_endblk) THEN
         i_startidx = 1
         i_endidx   = pt_patch%cells%end_idx(2,1)
       ELSE
         i_startidx = 1
         i_endidx   = nproma
       ENDIF

!$OMP WORKSHARE
       z_dttdiv(i_startidx:i_endidx,:,jb) = 0._wp
       z_help_c(i_startidx:i_endidx,:,jb) = 0._wp
!$OMP END WORKSHARE

     ENDDO
   ENDIF

!$OMP DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = 1,nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
         z_rhs(nlen+1:nproma,:,jb) = 0.0_wp
      ENDIF
      z_rhs(1:nlen,:,jb)=z_dttdiv(1:nlen,:,jb)-z_dt*asi_n*z_help_c(1:nlen,:,jb)
   ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

   CALL sync_patch_array(SYNC_C,pt_patch,z_rhs)

!-----------------------------------------------------------------------
! 2. decouple the right-hand-side: rhs_dcpl = inverse-eigenvec * rhs
!-----------------------------------------------------------------------

!$OMP PARALLEL
#if (defined (__SUNPRO_F95) || defined(__SX__)) && !defined (NOMPI)
!$OMP DO PRIVATE(jb, nlen, jc, jk, jm) ICON_OMP_DEFAULT_SCHEDULE
#else
!$OMP DO PRIVATE(jb, nlen, z_wrk_cv, z_wrk_vc) ICON_OMP_DEFAULT_SCHEDULE
#endif
   DO jb = 1,nblks_c

     IF (jb /= nblks_c) THEN
       nlen = nproma
     ELSE
       nlen = npromz_c
       z_rhs_dcpl(nlen+1:nproma,:,jb) = 0.0_wp
     ENDIF

#if (defined (__SUNPRO_F95) || defined(__SX__)) && !defined (NOMPI)
! sunf95 fails to OpenMP-parallelize the CALL to dgemm
! On the NEC, results with DGEMM are not invariant against the processor configuration

     z_rhs_dcpl(:,:,jb) = 0._wp

     DO jk = 1, nlev
       DO jm = 1, nlev
         DO jc = 1, nlen
           z_rhs_dcpl(jc,jk,jb) = z_rhs_dcpl(jc,jk,jb) + z_rhs(jc,jm,jb)*reigenvec(jm,jk)
         ENDDO
       ENDDO
     ENDDO

! MATMUL does not vectorize correctly
!     z_rhs_dcpl(1:nlen,:,jb) = MATMUL(z_wrk_cv(1:nlen,:),reigenvec(:,:))
#else

     z_wrk_cv(1:nlen,:) = z_rhs(1:nlen,:,jb)

     CALL dgemm(   & ! BLAS routine: C := alpha*op( A )*op( B ) + beta*C
                 'n',       & ! do not transpose A
                 'n',       & ! do not transpose B
                 nlen,      & ! number of rows of matrix A and C
                 nlev,      & ! number of columns of B' and C
                 nlev,      & ! number of columns of A and rows of B'
                 1._wp,     & ! alpha = 1.
                 z_wrk_cv,  & ! matrix A, unchanged on exit
                 nproma,    & ! leading dimension of A
                 reigenvec, & ! matrix B, unchanged on exit
                 nlev,      & ! leading dimension of B
                 0._wp,     & ! beta
                 z_wrk_vc,  & ! output
                 nproma     ) ! leading dimension of C

     z_rhs_dcpl(1:nlen,:,jb) = z_wrk_vc(1:nlen,:)
#endif

     z_dttdiv (:,:,jb) = 0._wp

   ENDDO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!-----------------------------------------------------------------------
! 3. solve 2d Helmholtz equations
!-----------------------------------------------------------------------
! set initial guess to zero, for modes larger than nmodes, this is already
! the result

! solve linear system

   DO jm = 1, nmodes ! solve the fast modes

      z_dtsq = z_dt*z_dt*asi_n*asi_n*eigenval(jm)

      CALL gmres( z_dttdiv(:,jm:jm,:),   &! x. Input is the first guess
                  lhs_div_eqn_2d,        &! sbr. calculating l.h.s.
                  pt_patch,              &! used for calculating l.h.s.
                  pt_int_state,          &! interpolation state
                  nblks_c,               &! number of blocks
                  npromz_c,              &! length of last block
                  z_dtsq,                &! used for calculating l.h.s.
                  z_rhs_dcpl(:,jm:jm,:), &! right hand side as input
                  si_rtol,               &! relative tolerance
                  .FALSE.,               &! NOT absolute tolerance
                  nmax_iter,             &! max. # of iterations to do
                  lmaxiter,              &! out: .true. = not converged
                  niter,                 &! out: # of iterations done
                  z_residual             &! out: the residual (array)
                  )                       !

      IF (lmaxiter) THEN
         CALL finish('GMRES solver: ','NOT YET CONVERGED !!')
      ENDIF
      IF (msg_level >= 5) THEN
         WRITE(string,'(a,i4,a,e20.10)') 'GMRES solver: iteration ', niter,  &
                                         ', residual = ', ABS(z_residual(niter))
         CALL message('',TRIM(string))
      ENDIF 

   ENDDO !mode loop

!-----------------------------------------------------------------------
! 4. get the updated second temporal derivative of divergence
!-----------------------------------------------------------------------

   i_startblk = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)

!$OMP PARALLEL
#if (defined (__SUNPRO_F95) || defined(__SX__)) && !defined (NOMPI)
!$OMP DO PRIVATE(jb, nlen, i_startidx, i_endidx, z_wrk_cv, jc, jk, jm, &
#else
!$OMP DO PRIVATE(jb, nlen, i_startidx, i_endidx, z_wrk_cv, z_wrk_vc, &
#endif
!$OMP    z_dtemp, z_dlnps) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = 1,nblks_c

     IF (jb /= nblks_c) THEN
       nlen = nproma
     ELSE
       nlen = npromz_c
     ENDIF

     z_wrk_cv(1:nlen,:) = z_dttdiv(1:nlen,:,jb)

#if (defined (__SUNPRO_F95) || defined(__SX__)) && !defined (NOMPI)
! sunf95 fails to OpenMP-parallelize the CALL to dgemm
! On the NEC, results with DGEMM are not invariant against the processor configuration

     z_dttdiv(:,:,jb) = 0._wp

     DO jk = 1, nlev
       DO jm = 1, nlev
         DO jc = 1, nlen
           z_dttdiv(jc,jk,jb) = z_dttdiv(jc,jk,jb) + z_wrk_cv(jc,jm)*eigenvec(jm,jk)
         ENDDO
       ENDDO
     ENDDO

! MATMUL does not vectorize correctly
!     z_dttdiv(1:nlen,:,jb) = MATMUL(z_wrk_cv(1:nlen,:),eigenvec(:,:))
#else

      CALL dgemm(   & ! BLAS routine: C := alpha*op( A )*op( B ) + beta*C
                 'n',       & ! do not transpose A
                 'n',       & ! do not transpose B
                 nlen,      & ! number of rows of matrix A and C
                 nlev,      & ! number of columns of B' and C
                 nlev,      & ! number of columns of A and rows of B'
                 1._wp,     & ! alpha = 1.
                 z_wrk_cv,  & ! matrix A, unchanged on exit
                 nproma,    & ! leading dimension of A
                 eigenvec,  & ! matrix B, unchanged on exit
                 nlev,    & ! leading dimension of B
                 0._wp,     & ! beta
                 z_wrk_vc,  & ! output
                 nproma     ) ! leading dimension of C

      z_dttdiv(1:nlen,:,jb) = z_wrk_vc(1:nlen,:)
#endif

!-----------------------------------------------------------------------
! 5. update the prognostic variables
!-----------------------------------------------------------------------

  IF ( .NOT. lshallow_water) THEN
      CALL conteq( z_dttdiv(1:nproma,1:nlev,jb), nproma, nlen, loperm, &! input
                   z_dtemp,z_dlnps )                                    ! output

!--------------------------------------------------------------
! 5.1 surface pressure and temperature
!--------------------------------------------------------------
!note: plus sign is used in the following 2 expressions
!      because there is a '-' hidden in subroutine conteq.

     ! Restrict update of the prognostic variables to the interior
     ! of nested domains
     IF (jb < i_startblk) THEN
        i_startidx = 1
        i_endidx   = 0
     ELSE IF (jb == i_startblk) THEN
        i_startidx = MAX(1,pt_patch%cells%start_idx(grf_bdywidth_c+1,1))
        i_endidx   = nproma
        IF (jb == nblks_c) i_endidx = npromz_c
      ELSE IF (jb == nblks_c) THEN
        i_startidx = 1
        i_endidx   = npromz_c
      ELSE
        i_startidx = 1
        i_endidx   = nproma
      ENDIF

      pt_new%temp(i_startidx:i_endidx,:,jb) = &
        pt_new%temp(i_startidx:i_endidx,:,jb) + z_dt*z_dtemp(i_startidx:i_endidx,:)
      pt_new%pres_sfc(i_startidx:i_endidx,jb) = &
        pt_new%pres_sfc(i_startidx:i_endidx,jb) + z_dt*apr*z_dlnps(i_startidx:i_endidx)

   ELSE  ! shallow water ...

     IF (jb < i_startblk) THEN
        i_startidx = 1
        i_endidx   = 0
     ELSE IF (jb == i_startblk) THEN
        i_startidx = MAX(1,pt_patch%cells%start_idx(grf_bdywidth_c+1,1))
        i_endidx   = nproma
        IF (jb == nblks_c) i_endidx = npromz_c
      ELSE IF (jb == nblks_c) THEN
        i_startidx = 1
        i_endidx   = npromz_c
      ELSE
        i_startidx = 1
        i_endidx   = nproma
      ENDIF

      pt_new%pres_sfc(i_startidx:i_endidx,jb) =   &
        pt_new%pres_sfc(i_startidx:i_endidx,jb) - &
        z_dt*sw_ref_height*z_dttdiv(i_startidx:i_endidx,1,jb)
   ENDIF

!--------------------------------------------------------------
! 5.2 wind
!--------------------------------------------------------------
!note: zrhs is "free" after step 2. Use it here as a working array

     IF (.NOT. lshallow_water) THEN
        CALL pgrad( z_dtemp, z_dlnps, nproma, nlen, & ! input
                    z_rhs(:,:,jb) )                   ! output
     ENDIF

   ENDDO !block loop
!$OMP END DO

!note: plus sign is used in next expression because zrhs has already got
!      a minus sign from subroutine conteq.

   IF (.NOT. lshallow_water) THEN
!$OMP DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
     DO jb = 1, nblks_c
       IF (jb /= nblks_c) THEN
         nlen = nproma
       ELSE
         nlen = npromz_c
       ENDIF
       z_help_c(1:nlen,:,jb) = z_dttpot(1:nlen,:,jb)  &
         &                   + asi_n*z_dt*z_rhs(1:nlen,:,jb)
     ENDDO
!$OMP END DO NOWAIT
   ELSE
!$OMP DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
     DO jb = 1, nblks_c
       IF (jb /= nblks_c) THEN
         nlen = nproma
       ELSE
         nlen = npromz_c
       ENDIF
       z_help_c(1:nlen,:,jb) = z_dttpot(1:nlen,:,jb)  &
         &  - asi_n*z_dt*grav*sw_ref_height*z_dttdiv(1:nlen,:,jb)
     ENDDO
!$OMP END DO NOWAIT
   ENDIF
!$OMP END PARALLEL

   CALL grad_fd_norm( z_help_c, pt_patch, z_help_e )

   ! Restrict update of vn to the interior of nestd domains
   i_startblk = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk, nblks_e

     IF (jb == i_startblk) THEN
       i_startidx = pt_patch%edges%start_idx(grf_bdywidth_e+1,1)
       i_endidx   = nproma
       IF (jb == nblks_e) i_endidx = npromz_e
     ELSE IF (jb == nblks_e) THEN
       i_startidx = 1
       i_endidx   = npromz_e
     ELSE
       i_startidx = 1
       i_endidx   = nproma
     ENDIF

     pt_new%vn(i_startidx:i_endidx,:,jb) = &
       pt_new%vn(i_startidx:i_endidx,:,jb) - z_dt*z_help_e(i_startidx:i_endidx,:,jb)

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   CALL sync_patch_array(SYNC_E,pt_patch,pt_new%vn)

 END SUBROUTINE add_si_correction_2d
  !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !>
   !!  calculate the lhs of the divergence equation.
   !!
   !!
   !! @par Revision History
   !!  Original version by Hui Wan, MPI-M (2007-12-17)
   !!  Restructuring of code by Almut Gassmann, MPI-M (2008-09-17)
   !!
  FUNCTION lhs_div_eqn_2d( p_div, pt_patch, pt_int_state, p_coeff ) RESULT(p_lhs)
!

  REAL(wp),INTENT(IN)         :: p_div(:,:,:)
  TYPE(t_patch),TARGET,INTENT(INOUT):: pt_patch
  TYPE(t_int_state), INTENT(IN) :: pt_int_state
  REAL(wp),INTENT(IN)         :: p_coeff

  REAL(wp)               :: p_lhs(SIZE(p_div,1),SIZE(p_div,2),SIZE(p_div,3))

! LOCAL VARIABLES:
  REAL(wp)      :: z_lapl    (SIZE(p_div,1),SIZE(p_div,2),SIZE(p_div,3))
  INTEGER       :: i_startblk, i_endblk, i_startidx, i_endidx
  INTEGER       :: jb, jc, nblks_c, npromz_c

!-----------------------------------------------------------------------
   nblks_c   = pt_patch%nblks_c
   npromz_c  = pt_patch%npromz_c

! Laplacian

   SELECT CASE(idiv_method)

   CASE(1)

     CALL nabla2_scalar( p_div, pt_patch, pt_int_state, z_lapl, 1, 1 )

   CASE(2,3)

     CALL nabla2_scalar_avg( p_div, pt_patch, pt_int_state,       &
                             pt_int_state%c_bln_avg, z_lapl, 1, 1 )

   END SELECT

! l.h.s. = ( 1 - coeff* nabla^2 ) D

   i_startblk = pt_patch%cells%start_blk(3,1)

! PARALLEL DO combined for efficiency reasons
!$OMP PARALLEL DO PRIVATE(jb, i_startidx, i_endidx)
   DO jb = i_startblk, nblks_c

     IF (jb == i_startblk) THEN
       i_startidx = pt_patch%cells%start_idx(3,1)
       i_endidx   = nproma
       IF (jb == nblks_c) i_endidx = npromz_c
     ELSE IF (jb == nblks_c) THEN
       i_startidx = 1
       i_endidx   = npromz_c
       ! Fill unused array elements with zero
       p_lhs(npromz_c+1:nproma,1,nblks_c)=0.0_wp
     ELSE
       i_startidx = 1
       i_endidx   = nproma
     ENDIF

     p_lhs(i_startidx:i_endidx,1,jb) = p_div(i_startidx:i_endidx,1,jb) - &
       p_coeff*z_lapl(i_startidx:i_endidx,1,jb)
   ENDDO
!$OMP END PARALLEL DO

   IF (l_limited_area .OR. pt_patch%id > 1) THEN
     ! In case of nested domains, the left-hand-side is set to zero in the
     ! two outermost cell rows
     i_startblk = pt_patch%cells%start_blk(1,1)
     i_endblk   = pt_patch%cells%end_blk(2,1)

! This loop is too short for OpenMP parallelization to be meaningful
     DO jb = i_startblk, i_endblk

       IF (jb == i_startblk) THEN
         i_startidx = pt_patch%cells%start_idx(1,1)
         i_endidx   = nproma
         IF (jb == i_endblk) i_endidx = pt_patch%cells%end_idx(2,1)
       ELSE IF (jb == i_endblk) THEN
         i_startidx = 1
         i_endidx   = pt_patch%cells%end_idx(2,1)
       ELSE
         i_startidx = 1
         i_endidx   = nproma
       ENDIF

       DO jc = i_startidx, i_endidx
         p_lhs(jc,1,jb) = 0._wp
       ENDDO
     ENDDO
   ENDIF

   CALL sync_patch_array(SYNC_C,pt_patch,p_lhs)


 END FUNCTION lhs_div_eqn_2d

!-------------------------------------------------------------------------
!
!

  !>
  !!
  SUBROUTINE add_si_correction_3d( si_coeff, si_rtol, p_dtime,&
                                   pt_patch, pt_int_state,    &
                                   pt_old, pt_now, pt_new  )

! !DESCRIPTION
!  This subroutine added the semi-implicit correction to the
!  estimated values of V, T and ps at time n+1 given by the
!  leapfrog scheme.
!  The central task is to solve a 3D Helmholtz equation of
!  ** the second temporal derivative ** of divergence.
!
! !REVISION HISTORY
!  Original version (solving 2D equations of divergence)
!   by Hui Wan (MPI-M, 2007-08-08)
!  Modifications by Hui Wan, MPI-M (2007-12-15)
!  - implemented the GMRES solver.
!  Modifications by Hui Wan, MPI-M (2008-01-25)
!  - Recoded after GME to solve 2D equations of the second temproal derivative
!    of divergence.
!  Modifications by Hui Wan, MPI-M (2008-01-31)
!  - added the option of directly solving the 3D equation of
!    the second temproal derivative of divergence.
!  Code restructuring by Almut Gassmann, MPI-M (2008-09-19)
!

  CHARACTER(len=*), PARAMETER :: routine = 'mo_si_correction:add_si_correction_3d'

  REAL(wp),INTENT(IN) :: si_coeff, si_rtol, p_dtime
  TYPE(t_patch), TARGET, INTENT(INOUT) :: pt_patch     ! patch
  TYPE(t_int_state),   INTENT(IN) :: pt_int_state ! interpolation state

  TYPE(t_hydro_atm_prog), INTENT(IN)   :: pt_old  ! prognostic variables at step n-1
  TYPE(t_hydro_atm_prog), INTENT(IN)   :: pt_now  ! prognostic variables at step n

  TYPE(t_hydro_atm_prog), INTENT(INOUT):: pt_new  ! input:  the estimated value at
                                                  !         step n+1 given by
                                                  !         leapfrog scheme.
                                                  ! output: with the semi-implicit
                                                  !         correction added.

  INTEGER  :: jk, jb, jc, je
  INTEGER  :: nblks_e, nblks_c, npromz_e, npromz_c, nlen
  INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx

  REAL(wp) :: z_dt     ! scaled time step in seconds
  REAL(wp) :: z_dtsq   !

  REAL(wp), DIMENSION (nproma,nlev,pt_patch%nblks_c) :: &
            & z_dttdiv, & ! temporal Laplacian of divergence
            & z_dttpot, & ! generalized potential
            & z_rhs,    & ! RHS
            & z_help_c    ! help array

  REAL(wp) :: z_dtemp (nproma,nlev)
  REAL(wp) :: z_dlnps (nproma)

  REAL(wp), DIMENSION (nproma,nlev,pt_patch%nblks_e) :: &
            & z_help_e       ! temporary array

! for the GMRES solver
  INTEGER, PARAMETER :: nmax_iter = 100  ! max. number of allowed iterations
  REAL(wp):: z_residual(nmax_iter)       ! residuum of every iteration
  LOGICAL :: lmaxiter !.TRUE. on output if nmax_iter was encountered
  INTEGER :: niter    ! number of actually needed iterations

  CHARACTER(len=MAX_CHAR_LENGTH) :: string

!-----------------------------------------------------------------------

   ! weighted time step (not to mix up with si_offctr)
   z_dt = p_dtime*si_coeff

   nblks_c   = pt_patch%nblks_c
   npromz_c  = pt_patch%npromz_c
   nblks_e   = pt_patch%nblks_e
   npromz_e  = pt_patch%npromz_e

!-----------------------------------------------------------------------
! 1. Calculate the right-hand-side of the 3D Helmholtz equation
!-----------------------------------------------------------------------
! 1.1 second temporal derivative of divergence
!-----------------------------------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = 1,nblks_e

      IF (jb /= nblks_e) THEN
         nlen = nproma
      ELSE
         nlen = npromz_e
      ENDIF

      DO jk = 1,nlev
         z_help_e(1:nlen,jk,jb) =   asi_n*pt_new%vn(1:nlen,jk,jb) &
                                  - 2._wp*pt_now%vn(1:nlen,jk,jb) &
                                  + asi_o*pt_old%vn(1:nlen,jk,jb)
      ENDDO

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   SELECT CASE (idiv_method)
   CASE (1)

      CALL div(z_help_e, pt_patch, pt_int_state, &  ! input
               z_dttdiv )! output

   CASE (2,3)

      CALL div_avg(z_help_e, pt_patch, pt_int_state, pt_int_state%c_bln_avg, &  ! input
                   z_dttdiv )! output

   END SELECT

!-----------------------------------------------------------------------
! 1.2 compute second temporal derivative of the generalized potential
!     from temperature and surface pressure.
!-----------------------------------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, z_dtemp, z_dlnps) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = 1,nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      z_dtemp(1:nlen,:) = - 2._wp*pt_now%temp(1:nlen,:,jb) &
                          + asi_n*pt_new%temp(1:nlen,:,jb) &
                          + asi_o*pt_old%temp(1:nlen,:,jb)

      z_dlnps(1:nlen) =  rpr*(- 2._wp*pt_now%pres_sfc(1:nlen,jb) &
                              + asi_n*pt_new%pres_sfc(1:nlen,jb) &
                              + asi_o*pt_old%pres_sfc(1:nlen,jb))

     CALL pgrad ( z_dtemp,z_dlnps,nproma,nlen,  & ! input
                  z_dttpot(1:nproma,1:nlev,jb))   ! output
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
!--------------------------------------------------------------
! 1.3 Laplacian of the second temporal derivative
!     of the generalized potential
!--------------------------------------------------------------

   SELECT CASE(idiv_method)
   CASE (1)

     CALL nabla2_scalar( z_dttpot, pt_patch, pt_int_state, z_help_c )

   CASE (2,3)

     CALL nabla2_scalar_avg( z_dttpot, pt_patch, pt_int_state, &
                             pt_int_state%c_bln_avg, z_help_c )

   END SELECT


!$OMP PARALLEL PRIVATE (i_startblk,i_endblk,i_startidx,i_endidx)

   IF (l_limited_area .OR. pt_patch%id > 1) THEN
     ! In case of nested domains, the divergence is set to zero in the
     ! two outermost cell rows

     i_startblk = pt_patch%cells%start_blk(1,1)
     i_endblk   = pt_patch%cells%end_blk(2,1)

     DO jb = i_startblk, i_endblk

       IF (jb == i_startblk) THEN
         i_startidx = pt_patch%cells%start_idx(1,1)
         i_endidx   = nproma
         IF (jb == i_endblk) i_endidx = pt_patch%cells%end_idx(2,1)
       ELSE IF (jb == i_endblk) THEN
         i_startidx = 1
         i_endidx   = pt_patch%cells%end_idx(2,1)
       ELSE
         i_startidx = 1
         i_endidx   = nproma
       ENDIF

!$OMP WORKSHARE
       z_dttdiv(i_startidx:i_endidx,:,jb) = 0._wp
       z_help_c(i_startidx:i_endidx,:,jb) = 0._wp
!$OMP END WORKSHARE

     ENDDO
   ENDIF

!$OMP DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
         z_rhs(nlen+1:nproma,:,jb) = 0.0_wp
         z_dttdiv(nlen+1:nproma,:,jb) = 0.0_wp
      ENDIF
      z_rhs(1:nlen,:,jb) = z_dttdiv(1:nlen,:,jb)  &
        &                - z_dt*asi_n*z_help_c(1:nlen,:,jb)
   ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

   CALL sync_patch_array(SYNC_C,pt_patch,z_rhs)

!=======================================================================
!-----------------------------------------------------------------------
! 3. solve 3d Helmholtz equations
!-----------------------------------------------------------------------

! solve linear system

   z_dtsq = z_dt*z_dt*asi_n*asi_n

   CALL gmres( z_dttdiv(:,:,:),     &! x. Input is the first guess
               lhs_div_eqn_3d,      &! sbr. calculating l.h.s.
               pt_patch,            &! used for calculating l.h.s.
               pt_int_state,        &! interpolation state
               nblks_c,             &! number of blocks
               npromz_c,            &! length of last block
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
   ELSE  IF (msg_level >= 5) THEN
      WRITE(string,'(a,i4,a,e20.10)') 'GMRES solver: iteration ', niter,  &
                                      ', residual = ', ABS(z_residual(niter))
      CALL message(TRIM(routine),TRIM(string))
   ENDIF

!======================================================================

!-----------------------------------------------------------------------
! 5. update the prognostic variables
!-----------------------------------------------------------------------

   i_startblk = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, i_startidx, i_endidx, jc, z_dtemp, z_dlnps) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = 1,nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      CALL conteq( z_dttdiv(1:nproma,1:nlev,jb),nproma,nlen,loperm,  & ! input
                   z_dtemp,z_dlnps )                                   ! output

!--------------------------------------------------------------
! 5.1 surface pressure and temperature
!--------------------------------------------------------------
!note: plus sign is used in the following 2 expressions
!      because there is a '-' hidden in subroutine conteq.


     ! Restrict update of the prognostic variables to the interior
     ! of nested domains
     IF (jb < i_startblk) THEN
        i_startidx = 1
        i_endidx   = 0
     ELSE IF (jb == i_startblk) THEN
        i_startidx = pt_patch%cells%start_idx(grf_bdywidth_c+1,1)
        i_endidx   = nproma
        IF (jb == nblks_c) i_endidx = npromz_c
      ELSE IF (jb == nblks_c) THEN
        i_startidx = 1
        i_endidx   = npromz_c
      ELSE
        i_startidx = 1
        i_endidx   = nproma
      ENDIF

      DO jc = i_startidx, i_endidx
        pt_new%temp(jc,:,jb)   = pt_new%temp(jc,:,jb) &
                                   + z_dt*z_dtemp(jc,:)
        pt_new%pres_sfc(jc,jb) = pt_new%pres_sfc(jc,jb) &
                                   + z_dt*apr*z_dlnps(jc)
      ENDDO

!--------------------------------------------------------------
! 5.2 wind
!--------------------------------------------------------------
!note: zrhs is "free" after step 2. Use it here as a working array

      CALL pgrad( z_dtemp,z_dlnps,nproma,nlen, & ! input
                  z_rhs(:,:,jb) )               ! output

  ENDDO !block loop
!$OMP END DO

  !note: plus sign is used in next expression because zrhs has already got
  !      a minus sign from subroutine conteq.

!$OMP DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1, nblks_c
     IF (jb /= nblks_c) THEN
        nlen = nproma
     ELSE
        nlen = npromz_c
     ENDIF
     z_help_c(1:nlen,:,jb) =  &
       &  z_dttpot(1:nlen,:,jb)+asi_n*z_dt*z_rhs(1:nlen,:,jb)
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  CALL grad_fd_norm( z_help_c, pt_patch, z_help_e )

  ! Restrict update of vn to the interior of nestd domains
   i_startblk = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk, nblks_e

     IF (jb == i_startblk) THEN
       i_startidx = pt_patch%edges%start_idx(grf_bdywidth_e+1,1)
       i_endidx   = nproma
       IF (jb == nblks_e) i_endidx = npromz_e
     ELSE IF (jb == nblks_e) THEN
       i_startidx = 1
       i_endidx   = npromz_e
     ELSE
       i_startidx = 1
       i_endidx   = nproma
     ENDIF

     DO je = i_startidx, i_endidx
       pt_new%vn(je,:,jb) = pt_new%vn(je,:,jb) - z_dt*z_help_e(je,:,jb)
     ENDDO

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   CALL sync_patch_array(SYNC_E,pt_patch,pt_new%vn)

 END SUBROUTINE add_si_correction_3d
!-----------------------------------------------------------------------
!
!
!

   !>
   !!  calculate the lhs of the divergence equation.
   !!
   !!
   !! @par Revision History
   !!  Original 2D version by Hui Wan, MPI-M (2007-12-17)
   !!  3D version by Hui Wan, MPI-M (2008-01-31)
   !!  Code restructuring by Almut Gassmann, MPI-M (2008-09-19)
   !!
   FUNCTION lhs_div_eqn_3d( p_x3d, pt_patch, pt_int_state, p_coeff ) RESULT(p_ax3d)
!

!arguments

  TYPE(t_patch),TARGET,INTENT(INOUT) :: pt_patch
  TYPE(t_int_state), INTENT(IN) :: pt_int_state

  REAL(wp),INTENT(IN) :: p_x3d(:,:,:)
  REAL(wp),INTENT(IN) :: p_coeff
  REAL(wp)            :: p_ax3d(SIZE(p_x3d,1),SIZE(p_x3d,2),SIZE(p_x3d,3) )

!local variables

  INTEGER :: jb, nblks_c, npromz_c, nlen
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

  REAL(wp) :: z_help_c(SIZE(p_x3d,1),SIZE(p_x3d,2),SIZE(p_x3d,3) )
  REAL(wp) :: z_dtemp(SIZE(p_x3d,1),SIZE(p_x3d,2))
  REAL(wp) :: z_dlnps(SIZE(p_x3d,1))

!-----------------------------------------------------------------------

   nblks_c   = pt_patch%nblks_c
   npromz_c  = pt_patch%npromz_c

! multiply by the structure matrix
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, z_dtemp, z_dlnps) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = 1,nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      CALL conteq( p_x3d(:,:,jb),nproma,nlen,loperm,  &! input
                   z_dtemp,z_dlnps )                   ! output
      CALL pgrad ( z_dtemp,z_dlnps,nproma,nlen,       &! input
                   p_ax3d(:,:,jb) )

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

! calculate the Laplacian

   SELECT CASE(idiv_method)
   CASE(1)

     CALL nabla2_scalar(p_ax3d, pt_patch, pt_int_state, z_help_c )

   CASE(2,3)

     CALL nabla2_scalar_avg(p_ax3d, pt_patch, pt_int_state, &
                            pt_int_state%c_bln_avg, z_help_c )

   END SELECT

! l.h.s. = ( 1 - coeff* nabla^2 ) D

   !HW: plus sign here because "conteq" has already added a "-".

   i_startblk = pt_patch%cells%start_blk(3,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk, nblks_c

     IF (jb == i_startblk) THEN
       i_startidx = pt_patch%cells%start_idx(3,1)
       i_endidx   = nproma
       IF (jb == nblks_c) i_endidx = npromz_c
     ELSE IF (jb == nblks_c) THEN
       i_startidx = 1
       i_endidx   = npromz_c
       p_ax3d(npromz_c+1:nproma,:,nblks_c) = 0.0_wp
     ELSE
       i_startidx = 1
       i_endidx   = nproma
     ENDIF
      p_ax3d(i_startidx:i_endidx,:,jb) = p_x3d(i_startidx:i_endidx,:,jb) + &
        p_coeff*z_help_c(i_startidx:i_endidx,:,jb)
   ENDDO
!$OMP END DO

   IF (l_limited_area .OR. pt_patch%id > 1) THEN
     ! In case of nested domains, the left-hand-side is set to zero in the
     ! two outermost cell rows
     i_startblk = pt_patch%cells%start_blk(1,1)
     i_endblk   = pt_patch%cells%end_blk(2,1)

     DO jb = i_startblk, i_endblk

       IF (jb == i_startblk) THEN
         i_startidx = pt_patch%cells%start_idx(1,1)
         i_endidx   = nproma
         IF (jb == i_endblk) i_endidx = pt_patch%cells%end_idx(2,1)
       ELSE IF (jb == i_endblk) THEN
         i_startidx = 1
         i_endidx   = pt_patch%cells%end_idx(2,1)
       ELSE
         i_startidx = 1
         i_endidx   = nproma
       ENDIF

!$OMP WORKSHARE
       p_ax3d(i_startidx:i_endidx,:,jb) = 0._wp
!$OMP END WORKSHARE

     ENDDO
   ENDIF
!$OMP END PARALLEL

   CALL sync_patch_array(SYNC_C, pt_patch, p_ax3d)


 END FUNCTION lhs_div_eqn_3d

!-------------------------------------------------------------------------
!
!

  !>
  !!   To calculate linearized increments of the temperature.
  !!
  !!   To calculate linearized increments of the temperature
  !!   and the logarithm of surface pressure for a specified divergence.
  !!   *conteq* is called from *inhysi* in the initial calculation of the
  !!            gravity-wave matrix *bb,* and from subroutines and in the
  !!            calculation of semi-implicit corrections.
  !!
  !! @par Revision History
  !!  A. J. Simmons, ECMWF, November 1981, original source
  !!  L. Kornblueh, MPI, May 1998, f90 rewrite
  !!  U. Schulzweida, MPI, May 1998, f90 rewrite
  !!  H. Wan, MPI, Aug 2007, changed the order of dummy arguments
  !!  to improve readability.
  !!
  SUBROUTINE conteq(p_div,k_bdim,k_len,lpperm,p_dtemp,p_dlnps)

! !INPUT VARIABLES:

  INTEGER, INTENT(IN)  :: k_bdim  ! dimension size of arrays
  INTEGER, INTENT(IN)  :: k_len   ! number of entries for which
                                  !calculation is performed
  LOGICAL, INTENT(IN)  :: lpperm  ! if .true., permute the order of the
                                  !indexes of the 2-dimensional arrays

  REAL(wp),INTENT(IN)  :: p_div(k_bdim,*) ! specified divergence

! !OUTPUT VARIABLES:

  REAL(wp),INTENT(OUT) :: p_dtemp(k_bdim,*) ! computed temperature increment
  REAL(wp),INTENT(OUT) :: p_dlnps(k_bdim)      ! computed increment of ln(ps)


  REAL(wp) :: z_a, z_dp, z_l
  INTEGER  :: ikp, jk, jl

!-------------------------------------------------------------------------

  IF (lpperm) THEN
    DO jl = 1, k_len
      p_dtemp(1,jl) = 0._wp
    END DO
  ELSE
    DO jl = 1, k_len
      p_dtemp(jl,1) = 0._wp
    END DO
  END IF

  DO jk = 1, nlevm1
    ikp  = jk + 1
    z_dp = delpr(jk)
    z_l  = aktlrd(jk)
    z_a  = altrcp(jk)

    IF (lpperm) THEN
      DO jl = 1, k_len
        p_dtemp(ikp,jl) = z_dp*p_div(jk,jl) + p_dtemp(jk,jl)
        p_dtemp(jk,jl)  = z_a*p_div(jk,jl) + z_l*p_dtemp(jk,jl)
      END DO
    ELSE
      DO jl = 1, k_len
        p_dtemp(jl,ikp) = z_dp*p_div(jl,jk) + p_dtemp(jl,jk)
        p_dtemp(jl,jk)  = z_a*p_div(jl,jk) + z_l*p_dtemp(jl,jk)
      END DO

    END IF
  END DO
  z_dp = delpr(nlev)
  z_l = aktlrd(nlev)
  z_a = altrcp(nlev)
  IF (lpperm) THEN
    DO jl = 1, k_len
      p_dlnps(jl)      = -rpr*(z_dp*p_div(nlev,jl)+p_dtemp(nlev,jl))
      p_dtemp(nlev,jl) = z_a*p_div(nlev,jl) + z_l*p_dtemp(nlev,jl)
    END DO
  ELSE
    DO jl = 1, k_len
      p_dlnps(jl)      = -rpr*(z_dp*p_div(jl,nlev)+p_dtemp(jl,nlev))
      p_dtemp(jl,nlev) = z_a*p_div(jl,nlev) + z_l*p_dtemp(jl,nlev)
    END DO
  END IF

  END SUBROUTINE conteq
!-------------------------------------------------------------------------
!
!

  !>
  !!  Calculate linearized geopotential and pressure gradient terms.
  !!
  !!  Calculate linearized geopotential and pressure gradient terms
  !!  for specified temperature and surface pressure fields.
  !!  *pgrad* is called from *inhysi* in the initial calculation
  !!  of the gravity-wave matrix *bb,* and from subroutine for
  !!  the calculation of semi-implicit corrections.
  !!
  SUBROUTINE pgrad(p_t,p_lnps,k_dim,k_len,p_grd)

! !REVISION HISTORY
!  A. J. Simmons, ECMWF, November 1981, original source
!  L. Kornblueh, MPI, May 1998, f90 rewrite
!  U. Schulzweida, MPI, May 1998, f90 rewrite
!
! !INPUT VARIABLES

   INTEGER :: k_dim, k_len

   REAL(wp) :: p_lnps(*)
   REAL(wp) :: p_t(k_dim,*)

! !OUTPUT VARIABLES

   REAL(wp) :: p_grd(k_dim,*) ! computed sum

! !LOCAL VARIABLES

  REAL(wp) :: z_a, z_b
  INTEGER :: ikp, jk, jl

!-------------------------------------------------------------------------

  !  Executable statements

!-- 1. Integrate hydrostatic equation

  z_a = ralphr(nlev)
  DO jl = 1, k_len
    p_grd(jl,nlev) = z_a*p_t(jl,nlev) + rdtr*p_lnps(jl)
  END DO

  DO jk = nlevm1, 1, -1
    ikp = jk + 1
    z_a = ralphr(jk)
    z_b = rlnmar(ikp)

    DO jl = 1, k_len
      p_grd(jl,jk) = z_a*p_t(jl,jk) + z_b*p_t(jl,ikp) + p_grd(jl,ikp)
    END DO

  END DO

  END SUBROUTINE pgrad

!-------------------------------------------------------------------------

END MODULE mo_si_correction

