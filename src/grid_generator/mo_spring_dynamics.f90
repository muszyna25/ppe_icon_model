!>
!! Perform grid optimization using spring dynamics.
!!
!! References:
!! Tomita,H. et al., 2001: Shallow Water Model on a Modified
!! Icosahedral Geodesic Grid by Using Spring Dynamics, JCP 174, 579-613
!! Tomita,H. et al., 2002: An Optimization of the Icosahedral Grid Modified by
!! Spring Dynamics
!!
!! @par Revision History
!! Initial Release by Almut Gassmann (2008-08-28)
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
MODULE mo_spring_dynamics
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  !
  !
  !
  USE mo_kind,           ONLY: wp
  USE mo_math_utilities,  ONLY: arc_length, t_cartesian_coordinates
  USE mo_base_datatypes, ONLY: t_spheres
  USE mo_math_constants, ONLY: pi
  USE mo_exception,      ONLY: message_text, message, finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::  spring_dynamics

CONTAINS

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Spring dynamics adapted from NICAM.
  !!
  !! Abortion and exit criteria are slightly
  !! changed compared to that. Here, it is measured whether the potential energy
  !! might create new maxima, which means that the method is unstable and the
  !! 'beta' must be decreased. If the method is stable, the exit criterion looks
  !! for that state of the system, when kinetic energy has dropped to one tenth
  !! of the maximum value.
  !! Here: d= 2pi/ sqrt(number of triangles) instead of in NICAM: d=2 pi/10*2^(l-1)
  !!
  !! @par Revision History
  !! Initial Release by Almut Gassmann (2008-09-28)
  !! Modification by Almut Gassmann (2009-02-02)
  !! - beta_spring comes via the namelist input
  !!
  SUBROUTINE spring_dynamics (sph, beta_spring, lopt)
    !
    TYPE(t_spheres), INTENT(inout)    :: sph
    REAL(wp), INTENT(IN) :: beta_spring
    LOGICAL,  INTENT(IN) :: lopt

    INTEGER :: i, i_e, i_h, maxit
    REAL (wp) :: len0, lambda, friction, dt, hori, maxekin, maxtest, ekin, test
#ifdef __SX__
    REAL (wp) :: z1, z2, z3
    TYPE(t_cartesian_coordinates), DIMENSION(0:sph%no_vertices-1) :: &
      & ej0, ej1, ej2, ej3, ej4, ej5, feder, rhs
#else
    TYPE(t_cartesian_coordinates) :: ej0, ej1, ej2, ej3, ej4, ej5, rhs
    TYPE(t_cartesian_coordinates), DIMENSION(0:sph%no_vertices-1) :: feder
#endif
    !-----------------------------------------------------------------------


    ! Choice of beta_spring here: choose beta_spring as large as possible for a
    ! good C-grid property, choose beta_spring as small as possible for similarly
    ! shaped (not so strongly deformed) triangles.
    ! An error message will appear, if the potential energy mesured in 'test'
    ! will create new maxima. Then, beta_spring should be chosen further from the
    ! outer limits of the interval [0.9 .. 1.11]
    ! beta_spring     = 1.11_wp ! good C-grid constraint (upper limit)
    ! beta_spring     = 0.90_wp ! good shape constraint
    lambda   = 2.0_wp*pi/(SQRT(REAL(sph%no_triangles,wp)*1.0_wp))
    len0     = beta_spring*lambda
    friction = 1.0_wp
    maxit    = 1000000
    maxekin  = 0.0_wp
    maxtest  = 0.0_wp

    ! initial velocities of vertices is zero
    DO i_h = 0, sph%no_vertices-1
      sph%vs(i_h)%vert_vel%x = 0.0_wp
    ENDDO

    DO i=1, maxit

      ekin = 0._wp
      test = 0._wp

      IF (lopt) THEN  ! Convergence acceleration
        IF (i <= 50) THEN
          dt = 1.6e-2_wp
        ELSE IF (i <= 150) THEN
          dt = 1.6e-2_wp*(1._wp+0.03_wp*(i-50))
        ELSE
          dt = 6.4e-2_wp
        ENDIF
      ELSE ! constant "time step"
        dt = 1.6e-2_wp
      ENDIF

!$OMP PARALLEL
#ifdef __SX__
!$OMP DO PRIVATE(z1, z2, z3)
      DO i_e = 0, sph%no_edges-1
        ! code inlined from arc_length in order to obtain vectorization on the NEC
        z1 = SQRT(DOT_PRODUCT(sph%es(i_e)%vertex0%vertex%x(1:3),&
          &                   sph%es(i_e)%vertex0%vertex%x(1:3)))
        z2 = SQRT(DOT_PRODUCT(sph%es(i_e)%vertex1%vertex%x(1:3),&
          &                   sph%es(i_e)%vertex1%vertex%x(1:3)))
        z3 = DOT_PRODUCT(sph%es(i_e)%vertex0%vertex%x(1:3),      &
          &              sph%es(i_e)%vertex1%vertex%x(1:3))/(z1*z2)
        IF (z3 > 1._wp )  z3 =  1._wp
        IF (z3 < -1._wp ) z3 = -1._wp

        sph%es(i_e)%edge_primal_arc = ACOS(z3)
      ENDDO
!$OMP END DO
#else
!$OMP DO
      DO i_e = 0, sph%no_edges-1
        sph%es(i_e)%edge_primal_arc=arc_length(sph%es(i_e)%vertex0%vertex,&
          & sph%es(i_e)%vertex1%vertex)
      ENDDO
!$OMP END DO
#endif

#ifdef __SX__
      ! Vectorized code
!$OMP DO
      DO i_h = 0, sph%no_vertices-1

        ! Unit vectors P0 --> PJ
        ej0(i_h)%x = sph%vs(i_h)%neighbor0%vertex%x-sph%vs(i_h)%vertex%x
        ej0(i_h)%x = ej0(i_h)%x/SQRT(DOT_PRODUCT(ej0(i_h)%x,ej0(i_h)%x))
        ej1(i_h)%x = sph%vs(i_h)%neighbor1%vertex%x-sph%vs(i_h)%vertex%x
        ej1(i_h)%x = ej1(i_h)%x/SQRT(DOT_PRODUCT(ej1(i_h)%x,ej1(i_h)%x))
        ej2(i_h)%x = sph%vs(i_h)%neighbor2%vertex%x-sph%vs(i_h)%vertex%x
        ej2(i_h)%x = ej2(i_h)%x/SQRT(DOT_PRODUCT(ej2(i_h)%x,ej2(i_h)%x))
        ej3(i_h)%x = sph%vs(i_h)%neighbor3%vertex%x-sph%vs(i_h)%vertex%x
        ej3(i_h)%x = ej3(i_h)%x/SQRT(DOT_PRODUCT(ej3(i_h)%x,ej3(i_h)%x))
        ej4(i_h)%x = sph%vs(i_h)%neighbor4%vertex%x-sph%vs(i_h)%vertex%x
        ej4(i_h)%x = ej4(i_h)%x/SQRT(DOT_PRODUCT(ej4(i_h)%x,ej4(i_h)%x))

        IF (sph%vs(i_h)%edge5%info%idx /= -1 ) THEN
          ej5(i_h)%x = sph%vs(i_h)%neighbor5%vertex%x-sph%vs(i_h)%vertex%x
          ej5(i_h)%x = ej5(i_h)%x/SQRT(DOT_PRODUCT(ej5(i_h)%x,ej5(i_h)%x))
        ELSE
          ej5(i_h)%x = 0._wp
        ENDIF

      ENDDO
!$OMP END DO

!$OMP DO
      DO i_h = 0, sph%no_vertices-1
        ! add all components for the spring force
        IF (sph%vs(i_h)%edge5%info%idx /= -1 ) THEN
          feder(i_h)%x = (sph%vs(i_h)%edge0%edge_primal_arc - len0)*ej0(i_h)%x +&
            &            (sph%vs(i_h)%edge1%edge_primal_arc - len0)*ej1(i_h)%x +&
            &            (sph%vs(i_h)%edge2%edge_primal_arc - len0)*ej2(i_h)%x +&
            &            (sph%vs(i_h)%edge3%edge_primal_arc - len0)*ej3(i_h)%x +&
            &            (sph%vs(i_h)%edge4%edge_primal_arc - len0)*ej4(i_h)%x +&
            &            (sph%vs(i_h)%edge5%edge_primal_arc - len0)*ej5(i_h)%x 
        ELSE
          feder(i_h)%x = (sph%vs(i_h)%edge0%edge_primal_arc - len0)*ej0(i_h)%x +&
            &            (sph%vs(i_h)%edge1%edge_primal_arc - len0)*ej1(i_h)%x +&
            &            (sph%vs(i_h)%edge2%edge_primal_arc - len0)*ej2(i_h)%x +&
            &            (sph%vs(i_h)%edge3%edge_primal_arc - len0)*ej3(i_h)%x +&
            &            (sph%vs(i_h)%edge4%edge_primal_arc - len0)*ej4(i_h)%x
        ENDIF
      ENDDO
!$OMP END DO

!$OMP DO PRIVATE(hori) REDUCTION(+:ekin,test)
      DO i_h = 0, sph%no_vertices-1
        ! RHS of spring equation
        rhs(i_h)%x = feder(i_h)%x-friction*sph%vs(i_h)%vert_vel%x

        ! solve the position equation
        sph%vs(i_h)%vertex%x = sph%vs(i_h)%vertex%x+dt*sph%vs(i_h)%vert_vel%x
        ! - > Normalize
        sph%vs(i_h)%vertex%x = sph%vs(i_h)%vertex%x/&
           SQRT(DOT_PRODUCT(sph%vs(i_h)%vertex%x,sph%vs(i_h)%vertex%x))

        ! solve for spring equation
        sph%vs(i_h)%vert_vel%x = sph%vs(i_h)%vert_vel%x + dt*rhs(i_h)%x
        hori = DOT_PRODUCT(sph%vs(i_h)%vert_vel%x,sph%vs(i_h)%vertex%x)
        ! - > Horizontalize
        sph%vs(i_h)%vert_vel%x = sph%vs(i_h)%vert_vel%x - &
                                 hori * sph%vs(i_h)%vertex%x

        ! kinetic energy for testing
        ekin = ekin + 0.5_wp*DOT_PRODUCT(sph%vs(i_h)%vert_vel%x,&
                                  sph%vs(i_h)%vert_vel%x)

        ! test criterion
        test = test + DOT_PRODUCT(feder(i_h)%x,feder(i_h)%x)
      ENDDO
!$OMP END DO
#else
!$OMP DO PRIVATE(ej0, ej1, ej2, ej3, ej4, ej5)
      DO i_h = 0, sph%no_vertices-1

        ! Unit vectors P0 --> PJ
        ej0%x = sph%vs(i_h)%neighbor0%vertex%x-sph%vs(i_h)%vertex%x
        ej0%x = ej0%x/SQRT(DOT_PRODUCT(ej0%x,ej0%x))
        ej1%x = sph%vs(i_h)%neighbor1%vertex%x-sph%vs(i_h)%vertex%x
        ej1%x = ej1%x/SQRT(DOT_PRODUCT(ej1%x,ej1%x))
        ej2%x = sph%vs(i_h)%neighbor2%vertex%x-sph%vs(i_h)%vertex%x
        ej2%x = ej2%x/SQRT(DOT_PRODUCT(ej2%x,ej2%x))
        ej3%x = sph%vs(i_h)%neighbor3%vertex%x-sph%vs(i_h)%vertex%x
        ej3%x = ej3%x/SQRT(DOT_PRODUCT(ej3%x,ej3%x))
        ej4%x = sph%vs(i_h)%neighbor4%vertex%x-sph%vs(i_h)%vertex%x
        ej4%x = ej4%x/SQRT(DOT_PRODUCT(ej4%x,ej4%x))

        IF (sph%vs(i_h)%edge5%info%idx /= -1 ) THEN
          ej5%x = sph%vs(i_h)%neighbor5%vertex%x-sph%vs(i_h)%vertex%x
          ej5%x = ej5%x/SQRT(DOT_PRODUCT(ej5%x,ej5%x))
          feder(i_h)%x = (sph%vs(i_h)%edge0%edge_primal_arc - len0)*ej0%x +&
            &            (sph%vs(i_h)%edge1%edge_primal_arc - len0)*ej1%x +&
            &            (sph%vs(i_h)%edge2%edge_primal_arc - len0)*ej2%x +&
            &            (sph%vs(i_h)%edge3%edge_primal_arc - len0)*ej3%x +&
            &            (sph%vs(i_h)%edge4%edge_primal_arc - len0)*ej4%x +&
            &            (sph%vs(i_h)%edge5%edge_primal_arc - len0)*ej5%x
        ELSE
          ej5%x = 0._wp
          feder(i_h)%x = (sph%vs(i_h)%edge0%edge_primal_arc - len0)*ej0%x +&
            &            (sph%vs(i_h)%edge1%edge_primal_arc - len0)*ej1%x +&
            &            (sph%vs(i_h)%edge2%edge_primal_arc - len0)*ej2%x +&
            &            (sph%vs(i_h)%edge3%edge_primal_arc - len0)*ej3%x +&
            &            (sph%vs(i_h)%edge4%edge_primal_arc - len0)*ej4%x
        ENDIF

      ENDDO
!$OMP END DO

! GZ: Splitting the loop here fixes a bug that used to be in the non-SX branch until November 2013.
!     Without the loop splitting, a feedback of the vertex motion on the neighbor vertices addressed
!     above occurs, leading to different results when the order of the grid points is changed.
!     The error became apparent because it caused an OpenMP race condition.

!$OMP DO PRIVATE(hori, rhs) REDUCTION(+:ekin,test)
      DO i_h = 0, sph%no_vertices-1
        ! RHS of spring equation
        rhs%x = feder(i_h)%x-friction*sph%vs(i_h)%vert_vel%x

        ! solve the position equation
        sph%vs(i_h)%vertex%x = sph%vs(i_h)%vertex%x+dt*sph%vs(i_h)%vert_vel%x
        ! - > Normalize
        sph%vs(i_h)%vertex%x = sph%vs(i_h)%vertex%x/&
           SQRT(DOT_PRODUCT(sph%vs(i_h)%vertex%x,sph%vs(i_h)%vertex%x))

        ! solve for spring equation
        sph%vs(i_h)%vert_vel%x = sph%vs(i_h)%vert_vel%x + dt*rhs%x
        hori = DOT_PRODUCT(sph%vs(i_h)%vert_vel%x,sph%vs(i_h)%vertex%x)
        ! - > Horizontalize
        sph%vs(i_h)%vert_vel%x = sph%vs(i_h)%vert_vel%x - &
                                 hori * sph%vs(i_h)%vertex%x

        ! kinetic energy for testing
        ekin = ekin + 0.5_wp*DOT_PRODUCT(sph%vs(i_h)%vert_vel%x,&
                                  sph%vs(i_h)%vert_vel%x)

        ! test criterion
        test = test + DOT_PRODUCT(feder(i_h)%x,feder(i_h)%x)
      ENDDO
!$OMP END DO
#endif
!$OMP END PARALLEL

      maxekin = MAX(ekin,maxekin)
      maxtest = MAX(test,maxtest)

      IF (i > 5 .and. test == maxtest) THEN
        CALL finish('mo_spring_dynamics', &
          & 'potential energy increases. Please adjust beta_spring in the range 0.9..1.11')
      ENDIF
      ! kinet energy should decrease to 1/10 of the maxium
      IF (i > 5 .and. ekin < 0.001_wp*maxekin ) THEN
          WRITE(message_text,'(a,i5, a)') ' iterated ', i, " times."
          CALL message ('spring_dynamics', TRIM(message_text))
         EXIT
      ENDIF
      IF (i == maxit) THEN
        CALL finish('mo_spring_dynamics','no convergence in grid iteration')
      ENDIF

    ENDDO

  END SUBROUTINE spring_dynamics


END MODULE mo_spring_dynamics


