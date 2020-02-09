!>
!! Flux and reconstruction limiter for vertical tracer transport
!!
!! This module contains various filters and flux limiters for 
!! vertical tracer transport.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision HistoryX1
!! Initial revision by Daniel Reinert, DWD (2010-03-04)
!! Modification by Daniel Reinert, DWD (2010-03-04)
!! - transferred all limiter to this new module
!! - implementation of flux limiter for FFSL-scheme (Miura)
!! Modification by Daniel Reinert, DWD (2010-10-06)
!! - implemented positive definite FCT-limiter for FFSL-scheme
!! Modification by Daniel Reinert, DWD (2012-05-04)
!! - removed slope limiter for horizontal transport, since they were 
!!   no longer used.
!! Modification by Daniel Reinert, DWD (2013-05-23)
!! - removed obsolete slope limiters for vertical MUSCL scheme.
!! Modification by Will Sawyer, CSCS (2015-07-14)
!! - added OpenACC support
!! Modification by Daniel Reinert, DWD (2018-11-29)
!! - Module mo_advection_limiter has been splitted into two parts, 
!!   in order to separate horizontal and vertical limiter.
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
MODULE mo_advection_vlimit

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_parallel_config,     ONLY: nproma, p_test_run
#ifdef _OPENACC
  USE mo_mpi,                 ONLY: i_am_accel_node
#endif

  IMPLICIT NONE

  PRIVATE


  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: vflx_limiter_pd
  PUBLIC :: v_limit_parabola_mo
  PUBLIC :: v_limit_parabola_sm
  PUBLIC :: v_limit_slope_mo
  PUBLIC :: v_limit_slope_sm
  PUBLIC :: v_limit_face_mo
  PUBLIC :: v_limit_face_sm

#if defined( _OPENACC )
#if defined(__ADVECTION_LIMITER_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
  LOGICAL, PARAMETER ::  acc_validate = .FALSE.   ! ONLY SET TO .TRUE. FOR VALIDATION PHASE
#endif

CONTAINS



  !-------------------------------------------------------------------------
  !>
  !! Positive definite flux limiter for vertical advection
  !!
  !! Positive definite Zalesak Flux-Limiter (Flux corrected transport).
  !! for the nonhydrostatic core. Only outward fluxes are re-scaled, in 
  !! order to maintain positive definiteness.
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !! - Harris, L. M. and P. H. Lauritzen (2010): A flux-form version of the
  !!   Conservative Semi-Lagrangian Multi-tracer transport scheme (CSLAM) on
  !!   the cubed sphere grid.  J. Comput. Phys., 230, 1215-1237
  !! - Smolarkiewicz, P. K., 1989: Comment on "A positive definite advection 
  !!   scheme obtained by nonlinear renormalization of the advective fluxes.", 
  !!   Mon. Wea. Rev., 117, 2626-2632
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2011-01-07)
  !!
  SUBROUTINE vflx_limiter_pd( p_dtime, p_cc, p_rhodz_now,             &
    &                         p_mflx_tracer_v, i_startidx, i_endidx,  &
    &                         slev, elev )

    REAL(wp), INTENT(IN) ::     &    !< time step [s]
      &  p_dtime

    REAL(wp), INTENT(IN) ::     &    !< advected cell centered variable at time (n)
      &  p_cc(:,:)                   !< dim: (nproma,nlev)
                                     !< [kg kg^-1]

    REAL(wp), INTENT(IN) ::     &    !< density times cell thickness at timestep n
      &  p_rhodz_now(:,:)            !< dim: (nproma,nlev)

    REAL(wp), INTENT(INOUT) ::  &    !< calculated vertical tracer mass flux
      &  p_mflx_tracer_v(:,:)        !< dim: (nproma,nlevp1)
                                     !< [kg m^-2 s^-1]

    INTEGER, INTENT(IN) ::      &    !< horizontal start and end index of DO loop
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN) ::      &    !< vertical start and end level
      &  slev, elev

    ! local
    !
    REAL(wp) ::                 &    !< fraction with which all outgoing fluxes 
      &  r_m(nproma,SIZE(p_cc,2))    !< of cell jk are multiplied 
                                     !< to guarantee positive definiteness

    REAL(wp) :: p_m(nproma)          !< sum of fluxes out of cell
                                     !< [kg m^-2]

    REAL(wp) :: z_signum             !< sign of mass flux
                                     !< >0: upward; <0: downward

    INTEGER  :: jk, jc               !< index of vert level, cell
    INTEGER  :: jkp1, jkm1

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: r_m,p_m
#endif
  !-------------------------------------------------------------------------


!$ACC DATA CREATE( r_m ), PCOPYIN( p_cc, p_rhodz_now ), PCOPY( p_mflx_tracer_v ), &
!$ACC IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_cc, p_mflx_tracer_v ), WAIT(1) IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

    IF (p_test_run) THEN
!$ACC KERNELS DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
      r_m = 0._wp
!$ACC END KERNELS
    ENDIF

    !
    ! 1. Compute total outward mass (loop over full levels)
    !
!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
    !$ACC LOOP GANG COLLAPSE(2) PRIVATE(p_m)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
        jkp1 = jk+1

        ! Sum of all outgoing fluxes out of cell jk
        p_m(jc) = p_dtime                               &
          &     * (MAX(0._wp,p_mflx_tracer_v(jc,jk))    &  ! upper half level
          &      - MIN(0._wp,p_mflx_tracer_v(jc,jkp1)) )   ! lower half level

        ! fraction with which all the fluxes out of cell jk must be multiplied 
        ! to guarantee no undershoot
        ! Nominator: maximum allowable mass loss \rho^n q^n
        r_m(jc,jk) = MIN(1._wp, (p_cc(jc,jk)*p_rhodz_now(jc,jk)) &
          &         /(p_m(jc) + dbl_eps) )

      ENDDO
    ENDDO
!$ACC END PARALLEL

    !
    ! 2. Limit outward fluxes (loop over half levels)
    !    Choose r_m depending on the sign of p_mflx_tracer_v
    !
!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
    !$ACC LOOP GANG VECTOR PRIVATE(jkm1,z_signum) COLLAPSE(2)
    DO jk = slev+1, elev
      DO jc = i_startidx, i_endidx
        jkm1 = jk-1

        ! NH:
        ! p_mflx_tracer_v(k-1/2) > 0: flux directed from cell k   -> k-1
        ! p_mflx_tracer_v(k-1/2) < 0: flux directed from cell k-1 -> k
        !
        z_signum = SIGN(1._wp,p_mflx_tracer_v(jc,jk))

        p_mflx_tracer_v(jc,jk) =  p_mflx_tracer_v(jc,jk)  * 0.5_wp      &
          &                       * ( (1._wp + z_signum) * r_m(jc,jk)   &
          &                       +   (1._wp - z_signum) * r_m(jc,jkm1) )
  
      ENDDO
    ENDDO
!$ACC END PARALLEL


!$ACC UPDATE HOST( p_mflx_tracer_v ), WAIT(1) IF (acc_validate .AND. i_am_accel_node .AND. acc_on)
!$ACC END DATA

  END SUBROUTINE vflx_limiter_pd



  !-------------------------------------------------------------------------
  !>
  !! Filter for PPM/PSM (3rd order) vertical advection (monotonic version)
  !!
  !! Removes over- and undershoots in first guess parabola by resetting the
  !! upper or lower interface values.
  !! Avoids non-physical over/undershoots in advected fields.
  !!
  !! Note that this limiter was coded assuming a pressure based vertical
  !! coordinate system. Nevertheless this limiter works for a height based
  !! vertical system, too. This is due to a 'wrong' computation of z_delta
  !! in the case of a height based coordinate system (i.e. z_delta is
  !! implicity multiplied by -1)
  !!
  !! In this version we do no longer make use of the monotonized slope 
  !! for detection of local extrema. Therefore this code can be used for 
  !! PSM as well, where the monotonized slope is not available.
  !!
  !! Literature
  !! Lin and Rood (1996), MWR, 124, 2046-2070
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-02-04)
  !! Modification by Daniel Reinert, DWD (2015-11-26)
  !! - do no longer rely on p_slope for detection of extrema. This 
  !!   allows usage of the limiter for PSM as well.
  !! Modification by Daniel Reinert, DWD (2019-01-11)
  !! - ability to preserve smooth extrema (lselective_limit = .TRUE.)
  !!
  SUBROUTINE v_limit_parabola_mo( p_ivlimit_selective, p_cc, p_face, &
    &                           p_face_up, p_face_low, i_startidx,   &
    &                           i_endidx, slev, elev )


    INTEGER, INTENT(IN) ::   &    !< avoids limiting of smooth extrema
      &  p_ivlimit_selective      !< if activated

    REAL(wp), INTENT(IN) ::  &     !< advected cell centered variable
      &  p_cc(:,:)                 !< dim: (nproma,nlev)

    REAL(wp), INTENT(IN) ::  &     !< reconstructed face values of the advected field
      &  p_face(:,:)               !< dim: (nproma,nlevp1)

    REAL(wp), INTENT(INOUT) :: &   !< final face value (upper face, height based)
      &  p_face_up(:,:)            !< dim: (nproma,nlevp)

    REAL(wp), INTENT(INOUT) :: &   !< final face value (lower face, height based)
      &  p_face_low(:,:)           !< dim: (nproma,nlevp)

    INTEGER, INTENT(IN)     :: &   !< horizontal start and end index of DO loop
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN)     :: &   !< vertical start and end index of DO loop
      &  slev, elev

    INTEGER  :: jc, jk             !< index of cell and vertical level
    INTEGER  :: ikp1               !< vertical level plus one

    REAL(wp) :: z_delta            !< lower minus upper face value
    REAL(wp) :: z_a6i              !< curvature of parabola
    LOGICAL  :: is_main_crit       !< is main criterion for limiter activation TRUE
    LOGICAL  :: lselective_limit   !< selective limitation yes/no

    !-----------------------------------------------------------------------

    ! selective limitation yes or no
    lselective_limit = (p_ivlimit_selective == 1)

!$ACC DATA PCOPYIN( p_cc, p_face ), PCOPY( p_face_up, p_face_low ), &
!$ACC      IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_cc, p_face ), WAIT(1) IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
    !$ACC LOOP GANG VECTOR PRIVATE(ikp1, z_delta, z_a6i, is_main_crit) COLLAPSE(2)
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

        ! index of bottom half level
        ikp1 = jk+1 ! WS: put inside inner loop to collapse both loops

        z_delta   = p_face(jc,ikp1) - p_face(jc,jk)
        z_a6i     = 6._wp * (p_cc(jc,jk) - 0.5_wp * (p_face(jc,jk) + p_face(jc,ikp1)))

        ! main criterion upon which it is decided whether an extremum 
        ! is spurious and whether the limiter is activated. 
        is_main_crit = ABS(z_delta) < ABS(z_a6i)

        ! detect spurious extremum
        !
        DETECT:IF (isExtremumSpurious(jc,jk,is_main_crit,z_a6i,p_cc(jc,jk),p_face,              & 
          &                           slev,MIN(elev+1,UBOUND(p_face,2)),lselective_limit) ) THEN

          !
          ! parabola must be modified to remove local extremum
          !

          ! if cell average presents a local extremum, replace parabola 
          ! by piecewise constant function
          IF ( ((p_cc(jc,jk) - p_face(jc,ikp1))*(p_cc(jc,jk)-p_face(jc,jk))) > 0._wp) THEN
            p_face_up(jc,jk)  = p_cc(jc,jk)
            p_face_low(jc,jk) = p_cc(jc,jk)

          ELSE 
            !
            ! monotonize parabola by modifying one of the edge values
            IF (z_delta * z_a6i > z_delta * z_delta) THEN
              p_face_up(jc,jk)  = 3._wp*p_cc(jc,jk) - 2._wp*p_face(jc,ikp1)
              p_face_low(jc,jk) = p_face(jc,ikp1)

            ELSE IF (z_delta * z_a6i < -1._wp * (z_delta * z_delta)) THEN
              p_face_up(jc,jk)  = p_face(jc,jk)
              p_face_low(jc,jk) = 3._wp*p_cc(jc,jk) - 2._wp*p_face(jc,jk)

            ELSE
              ! necessary if z_delta and z_a6i become very tiny.
              p_face_up(jc,jk)  = p_face(jc,jk)
              p_face_low(jc,jk) = p_face(jc,ikp1)
            ENDIF
            !
          ENDIF
          !
        ELSE
          ! no monotonization required
          p_face_up(jc,jk)  = p_face(jc,jk)
          p_face_low(jc,jk) = p_face(jc,ikp1)
        ENDIF DETECT

      END DO  ! jc

    END DO  ! jk
!$ACC END PARALLEL

!$ACC UPDATE HOST( p_face_up, p_face_low ), WAIT(1) IF( acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA

  END SUBROUTINE v_limit_parabola_mo




  !-------------------------------------------------------------------------
  !>
  !! Filter for PPM/PSM (3rd order) vertical advection (semi-monotonic version)
  !!
  !! Removes undershoots in first guess parabola by resetting either the
  !! upper or lower interface value.
  !! Avoids non-physical undershoots in advected fields.
  !!
  !! Note that this limiter was coded assuming a pressure based vertical
  !! coordinate system. Nevertheless this limiter works for a height based
  !! vertical system, too - without any modifications.
  !!
  !! Literature
  !! Lin and Rood (1996), MWR, 124, 2046-2070
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-02-04)
  !!
  SUBROUTINE v_limit_parabola_sm( p_ivlimit_selective, p_cc, p_face, &
    &                           p_face_up, p_face_low, i_startidx,   &
    &                           i_endidx, slev, elev     )


    INTEGER, INTENT(IN) ::   &    !< avoids limiting of smooth extrema
      &  p_ivlimit_selective      !< if activated

    REAL(wp), INTENT(IN) ::  &     !< advected cell centered variable
      &  p_cc(:,:)                 !< dim: (nproma,nlev)

    REAL(wp), INTENT(IN) ::  &     !< reconstructed face values of the advected field
      &  p_face(:,:)               !< dim: (nproma,nlevp1)

    REAL(wp), INTENT(INOUT) :: &   !< final face value (upper face, height based)
      &  p_face_up(:,:)            !< dim: (nproma,nlevp)

    REAL(wp), INTENT(INOUT) :: &   !< final face value (lower face, height based)
      &  p_face_low(:,:)           !< dim: (nproma,nlevp)

    INTEGER, INTENT(IN)     :: &   !< horizontal start and end index of DO loop
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN)     :: &   !< vertical start and end index of DO loop
      &  slev, elev

    REAL(wp) :: z_delta            !< lower minus upper face value
    REAL(wp) :: z_a6i              !< curvature of parabola

    INTEGER  :: jc, jk             !< index of cell and vertical level
    INTEGER  :: ikp1               !< vertical level plus one
    LOGICAL  :: is_main_crit       !< is main criterion for limiter activation TRUE
    LOGICAL  :: lselective_limit   !< selective limitation yes/no

  !-------------------------------------------------------------------------

    ! selective limitation yes or no
    lselective_limit = (p_ivlimit_selective == 1)

!$ACC DATA PCOPYIN( p_cc, p_face ), PCOPYOUT( p_face_up, p_face_low ), &
!$ACC      IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_cc, p_face ), WAIT(1) IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
    !$ACC LOOP GANG VECTOR PRIVATE(ikp1, z_delta, z_a6i, is_main_crit) COLLAPSE(2)
    DO jk = slev, elev

      DO jc = i_startidx, i_endidx

        ! index of bottom half level
        ikp1 = jk+1 ! WS: put inside inner loop to collapse both loops

        z_delta   = p_face(jc,ikp1) - p_face(jc,jk)
        z_a6i     = 6._wp * (p_cc(jc,jk)                      &
          &       - 0.5_wp * (p_face(jc,jk) + p_face(jc,ikp1)))

        ! main criterion upon which it is decided whether an undershoot 
        ! is spurious and whether the limiter is activated. 
        is_main_crit = ABS(z_delta) < -1._wp*z_a6i


        ! detect spurious undershoots
        !
        DETECT:IF (isExtremumSpurious(jc,jk,is_main_crit,z_a6i,p_cc(jc,jk),p_face,              & 
          &                           slev,MIN(elev+1,UBOUND(p_face,2)),lselective_limit) ) THEN

          !
          ! parabola must be modified to remove local undershoots
          !

          ! if cell average presents a local minimum, replace parabola 
          ! by piecewise constant function
          IF (p_cc(jc,jk) < MIN(p_face(jc,jk),p_face(jc,ikp1)) ) THEN
            p_face_up(jc,jk)  = p_cc(jc,jk)
            p_face_low(jc,jk) = p_cc(jc,jk)

          ELSE
            !
            ! monotonize parabola by modifying one of the edge values
            IF (p_face(jc,jk) > p_face(jc,ikp1)) THEN
              p_face_up(jc,jk)  = 3._wp*p_cc(jc,jk) - 2._wp*p_face(jc,ikp1)
              p_face_low(jc,jk) = p_face(jc,ikp1)

            ELSE
              p_face_up(jc,jk)  = p_face(jc,jk)
              p_face_low(jc,jk) = 3._wp*p_cc(jc,jk) - 2._wp*p_face(jc,jk)

            ENDIF
            !
          ENDIF
          !
        ELSE
          p_face_up(jc,jk)  = p_face(jc,jk)
          p_face_low(jc,jk) = p_face(jc,ikp1)
        ENDIF DETECT

      END DO

    END DO
!$ACC END PARALLEL

!$ACC UPDATE HOST( p_face_up, p_face_low ), WAIT(1) IF( acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA

  END SUBROUTINE v_limit_parabola_sm



  !-------------------------------------------------------------------------
  !>
  !! Detect if the sub-grid reconstruction has a spurious extremum
  !!
  !! Detect if the sub-grid reconstruction has a spurious extremum.
  !! In the default case the extremum is deemed spurious, if it occurs 
  !! in the cell interior, i.e. if 0<\zeta_ext<1, where \zeta_ext denotes 
  !! the location of the extremum in dimnesionless coordinates.
  !! This information is provided by the LOGICAL variable is_main_crit.
  !!
  !! If lselective_limit=.TRUE., the detection algorithm is less restrictive 
  !! and avoids clipping of smooth extrema. At least one of 5 additional constraints 
  !! must be true for the extremum to be deemed spurious.  
  !!
  !! Literature
  !! Zerroukat et al. (2005), QJRMS, 131, 2923-2936
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2019-01-11)
  !!
  LOGICAL FUNCTION isExtremumSpurious(jc, jk, is_main_crit, z_a6i, p_cc, p_face, &
    &                                 jk_min, jk_max, lselective_limit)
!$ACC ROUTINE SEQ

    INTEGER , INTENT(IN) :: jc, jk           !< cell at hand
    LOGICAL , INTENT(IN) :: is_main_crit     !< is main criterion for limiter activation TRUE 
    REAL(wp), INTENT(IN) :: z_a6i            !< curvature of parabola
    REAL(wp), INTENT(IN) :: p_cc             !< cell average at given point (jc,jk)
    REAL(wp), INTENT(IN) :: p_face(:,:)      !< face values for vertical column at hand
                                             !  only required for jk+1,jk+2,jk,jk-1,jk-2 at given jc
    LOGICAL , INTENT(IN) :: lselective_limit !< distinguishes between selective 
                                             !  and non-selective limiting
    INTEGER, INTENT(IN)  :: jk_min           !< minimum half level up to which p_face is filled with 
                                             !< meaningful values 
    INTEGER, INTENT(IN)  :: jk_max           !< maximum half level up to which p_face is filled with 
                                             !< meaningful values
                                             !< Accessing indices beyond this range leads to 
                                             !< non-reporducible results.
    ! local
    LOGICAL  :: is_crit(5)          ! additional criteria of selective limiter
    INTEGER  :: jkp1, jkp2, jkp3, jkm1, jkm2  ! neighbour indices
    REAL(wp) :: loc_extremum        ! location of extreme value in dimensionless coordinate
    REAL(wp) :: val_extremum        ! extreme value
    REAL(wp) :: z_a1                ! linear coefficient of parabolic interpolant
    !-----------------------------------------------------------------------

    IF (is_main_crit.AND.lselective_limit) THEN

      jkp1 = MIN(jk+1,jk_max)
      jkp2 = MIN(jk+2,jk_max)
      jkp3 = MIN(jk+3,jk_max)
      jkm1 = MAX(jk-1,jk_min)
      jkm2 = MAX(jk-2,jk_min)

      z_a1  = 6._wp * p_cc - 2._wp*p_face(jc,jk) - 4._wp*p_face(jc,jkp1)
      is_crit(1) = ((p_face(jc,jkp1) - p_face(jc,jkp2)) * (p_face(jc,jkm1) - p_face(jc,jk))  ) >= 0._wp
      is_crit(2) = ((p_face(jc,jkp2) - p_face(jc,jkp3)) * (p_face(jc,jkp1) - p_face(jc,jkp2))) <= 0._wp
      is_crit(3) = ((p_face(jc,jkm1) - p_face(jc,jk)  ) * (p_face(jc,jkm2) - p_face(jc,jkm1))) <= 0._wp
      is_crit(4) = ((p_face(jc,jkp1) - p_face(jc,jkp2)) * z_a1 )                               <= 0._wp
      !
      ! determine location of extremum
      loc_extremum = 0.5_wp*(1._wp + (p_face(jc,jk) - p_face(jc,jkp1))/(SIGN(ABS(z_a6i)+dbl_eps,z_a6i)))
      ! determine extremum
      val_extremum = p_cc - (p_face(jc,jk) - p_face(jc,jkp1))*(0.5_wp - loc_extremum) &
          &         - z_a6i*(1._wp/6._wp - loc_extremum + loc_extremum*loc_extremum)
      is_crit(5) = val_extremum < 0._wp

      ! result
      isExtremumSpurious = ANY(is_crit(:))

    ELSE
      ! result
      isExtremumSpurious = is_main_crit
    ENDIF

  END FUNCTION isExtremumSpurious



  !-------------------------------------------------------------------------
  !>
  !! Monotonicity preserving slope limiter for PPM after Lin et al (1994)
  !!
  !! Monotonicity preserving slope limiter after Lin et al (1994).
  !! Guarantees e.g. that reconstructed PPM edge values lie in 
  !! the range of values defined by neighbouring cells. 
  !!
  !! Literature
  !! Lin et al. (1994), MWR, 122, 1575-1593
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-10-28)
  !! - moved here from mo_advection_vflux
  !!
  SUBROUTINE v_limit_slope_mo( p_cc, i_startidx, i_endidx, slev, elev, slope )


    REAL(wp), INTENT(IN) ::  &     !< advected cell centered variable
      &  p_cc(:,:)                 !< dim: (nproma,nlev)


    INTEGER, INTENT(IN)     :: &   !< horizontal start and end index of DO loop
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN)     :: &   !< vertical start and end index of DO loop
      &  slev, elev

    REAL(wp), INTENT(INOUT) :: &   !< advected cell centered variable
      &  slope(:,:)                !< dim: (nproma,nlev)

    REAL(wp) :: p_cc_min, p_cc_max    !< 3-point max/min values

    INTEGER  :: jc, jk                !< index of cell and vertical level
    INTEGER  :: ikp1, ikm1            !< vertical level plus-minus one

  !-------------------------------------------------------------------------

!$ACC DATA PCOPYIN( p_cc ), PCOPYOUT( slope ), IF( i_am_accel_node .AND. acc_on )

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR PRIVATE(ikm1, ikp1, p_cc_min, p_cc_max) COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! index of top half level
        ikm1    = jk - 1
        ! index of bottom half level
        ikp1    = MIN( jk+1, elev )

        ! equivalent formulation of Colella and Woodward (1984) slope limiter 
        ! following Lin et al (1994).
        p_cc_min = MIN(p_cc(jc,ikm1),p_cc(jc,jk),p_cc(jc,ikp1))
        p_cc_max = MAX(p_cc(jc,ikm1),p_cc(jc,jk),p_cc(jc,ikp1))
        slope(jc,jk) = SIGN(                                                   &
          &            MIN( ABS(slope(jc,jk)), 2._wp*(p_cc(jc,jk)-p_cc_min),   &
          &                                    2._wp*(p_cc_max-p_cc(jc,jk)) ), &
          &            slope(jc,jk) )
      END DO  ! jc

    END DO  ! jk 
!$ACC END PARALLEL

!$ACC END DATA
  END SUBROUTINE v_limit_slope_mo



  !-------------------------------------------------------------------------
  !>
  !! Semi-monotonic slope limiter for PPM after Lin et al (1994)
  !!
  !! Semi-monotonic slope limiter after Lin et al (1994).
  !! Guarantees e.g. that reconstructed PPM edge values do not 
  !! fall below the range of values defined by neighbouring cells. 
  !!
  !! Literature
  !! Lin et al. (1994), MWR, 122, 1575-1593
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-10-28)
  !!
  SUBROUTINE v_limit_slope_sm( p_cc, i_startidx, i_endidx, slev, elev, slope )


    REAL(wp), INTENT(IN) ::  &     !< advected cell centered variable
      &  p_cc(:,:)                 !< dim: (nproma,nlev)


    INTEGER, INTENT(IN)     :: &   !< horizontal start and end index of DO loop
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN)     :: &   !< vertical start and end index of DO loop
      &  slev, elev

    REAL(wp), INTENT(INOUT) :: &   !< advected cell centered variable
      &  slope(:,:)                !< dim: (nproma,nlev)

    REAL(wp) :: p_cc_min           !< 3-point min values

    INTEGER  :: jc, jk                !< index of cell and vertical level
    INTEGER  :: ikp1, ikm1            !< vertical level plus-minus one

  !-------------------------------------------------------------------------

!$ACC DATA PCOPYIN( p_cc ), PCOPYOUT( slope ), IF( i_am_accel_node .AND. acc_on )

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR PRIVATE(ikm1, ikp1, p_cc_min) COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! index of top half level
        ikm1    = jk - 1
        ! index of bottom half level
        ikp1    = MIN( jk+1, elev )

        ! equivalent formulation of Colella and Woodward (1984) slope limiter 
        ! following Lin et al (1994).
        p_cc_min = MIN(p_cc(jc,ikm1),p_cc(jc,jk),p_cc(jc,ikp1))
        slope(jc,jk) = SIGN(                                                   &
          &            MIN( ABS(slope(jc,jk)), 2._wp*(p_cc(jc,jk)-p_cc_min) ), &
          &            slope(jc,jk) )
      END DO  ! jc

    END DO  ! jk 
!$ACC END PARALLEL

!$ACC END DATA
  END SUBROUTINE v_limit_slope_sm



  !-------------------------------------------------------------------------
  !>
  !! Limit reconstructed face values.
  !!
  !! Make sure that face values lie within the range of values 
  !! in the neighbouring cells.
  !!
  !!
  !! Literature
  !! Blossey, P. N. and D. R. Durran (2008), JCP, 227, 5160-5183
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-11-26)
  !!
  SUBROUTINE v_limit_face_mo( p_cc, p_face, i_startidx, i_endidx, slev, elev )


    REAL(wp), INTENT(IN) ::    &   !< advected cell centered variable
      &  p_cc(:,:)                 !< dim: (nproma,nlev)

    REAL(wp), INTENT(INOUT) ::  &  !< reconstructed face values of the advected field
      &  p_face(:,:)               !< dim: (nproma,nlevp1)

    INTEGER, INTENT(IN)     ::  &  !< horizontal start and end index of DO loop
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN)     ::  &  !< vertical start and end index of DO loop
      &  slev, elev

    REAL(wp) :: p_cc_min, p_cc_max    !< 2-point max/min values

    INTEGER  :: jc, jk                !< index of cell and vertical level
    INTEGER  :: ikp1                  !< vertical level plus one

  !-------------------------------------------------------------------------

!$ACC DATA PCOPYIN( p_cc ), PCOPYOUT( p_face ), IF( i_am_accel_node .AND. acc_on )

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR PRIVATE(ikp1, p_cc_min, p_cc_max) COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ikp1 = jk+1

        ! make sure that face values lie within the range of values 
        ! in the neighbouring cells 
        p_cc_min = MIN(p_cc(jc,jk),p_cc(jc,ikp1))
        p_cc_max = MAX(p_cc(jc,jk),p_cc(jc,ikp1))

        p_face(jc,ikp1)= MIN(p_cc_max,MAX(p_cc_min,p_face(jc,ikp1)))

      END DO  ! jc

    END DO  ! jk 
!$ACC END PARALLEL

!$ACC END DATA
  END SUBROUTINE v_limit_face_mo




  !-------------------------------------------------------------------------
  !>
  !! One-sided limitation of reconstructed face values.
  !!
  !! Make sure that face values do not fall below the range of values 
  !! in the neighbouring cells.
  !!
  !!
  !! Literature
  !! Blossey, P. N. and D. R. Durran (2008), JCP, 227, 5160-5183
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-11-26)
  !!
  SUBROUTINE v_limit_face_sm( p_cc, p_face, i_startidx, i_endidx, slev, elev )


    REAL(wp), INTENT(IN) ::    &   !< advected cell centered variable
      &  p_cc(:,:)                 !< dim: (nproma,nlev)

    REAL(wp), INTENT(INOUT) ::  &  !< reconstructed face values of the advected field
      &  p_face(:,:)               !< dim: (nproma,nlevp1)

    INTEGER, INTENT(IN)     ::  &  !< horizontal start and end index of DO loop
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN)     ::  &  !< vertical start and end index of DO loop
      &  slev, elev

    REAL(wp) :: p_cc_min           !< 2-point min values

    INTEGER  :: jc, jk             !< index of cell and vertical level
    INTEGER  :: ikp1               !< vertical level plus one

  !-------------------------------------------------------------------------

!$ACC DATA PCOPYIN( p_cc ), PCOPYOUT( p_face ), IF( i_am_accel_node .AND. acc_on )

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR PRIVATE(ikp1, p_cc_min) COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ikp1 = jk+1

        ! make sure that face values do not fall below the range of values 
        ! in the neighbouring cells 
        p_cc_min = MIN(p_cc(jc,jk),p_cc(jc,ikp1))

        p_face(jc,ikp1)= MAX(p_cc_min,p_face(jc,ikp1))

      END DO  ! jc

    END DO  ! jk 
!$ACC END PARALLEL

!$ACC END DATA
  END SUBROUTINE v_limit_face_sm

END MODULE mo_advection_vlimit

