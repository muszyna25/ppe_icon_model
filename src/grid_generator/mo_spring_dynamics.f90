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
!! @par Copyright
!! 2002-2008 by DWD and MPI-M
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
  USE mo_base_geometry,  ONLY: arc_length, t_cartesian_coordinates
  USE mo_base_datatypes, ONLY: t_spheres
  USE mo_math_constants, ONLY: pi
  USE mo_exception,      ONLY: message_text, message, finish

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

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
  SUBROUTINE spring_dynamics (sph, beta_spring)
    !
    TYPE(t_spheres), INTENT(inout)    :: sph
    REAL(wp) :: beta_spring

    INTEGER :: i, i_e, i_h, maxit
    REAL (wp) :: len0, lambda, friction, dt, hori, ekin, test, &
      & maxekin, maxtest
#ifdef __SX__
    REAL (wp) :: z1, z2, z3
    TYPE(t_cartesian_coordinates), DIMENSION(0:sph%no_vertices-1) :: &
      & ej0, ej1, ej2, ej3, ej4, ej5, feder, rhs
#else
    TYPE(t_cartesian_coordinates) :: ej0, ej1, ej2, ej3, ej4, ej5, feder, rhs
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
    dt       = 1.6e-2_wp
    maxit    = 1000000
    maxekin  = 0.0_wp
    maxtest  = 0.0_wp

    ! initial velocities of vertices is zero
    DO i_h = 0, sph%no_vertices-1
      sph%vs(i_h)%vert_vel%x = 0.0_wp
    ENDDO

    DO i=1, maxit

      DO i_e = 0, sph%no_edges-1
#ifdef __SX__
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
#else
        sph%es(i_e)%edge_primal_arc=arc_length(sph%es(i_e)%vertex0%vertex,&
          & sph%es(i_e)%vertex1%vertex)
#endif
      ENDDO

      test = 0.0_wp
      ekin = 0.0_wp

#ifdef __SX__
      ! Vectorized code
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
        ekin = 0.5_wp*DOT_PRODUCT(sph%vs(i_h)%vert_vel%x,&
                                  sph%vs(i_h)%vert_vel%x) + ekin

        ! test criterion
        test = DOT_PRODUCT(feder(i_h)%x,feder(i_h)%x) + test
      ENDDO
#else
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
        ENDIF

        ! add all components for the spring force
        feder%x = (sph%vs(i_h)%edge0%edge_primal_arc - len0)*ej0%x +&
          & (sph%vs(i_h)%edge1%edge_primal_arc - len0)*ej1%x +&
          & (sph%vs(i_h)%edge2%edge_primal_arc - len0)*ej2%x +&
          & (sph%vs(i_h)%edge3%edge_primal_arc - len0)*ej3%x +&
          & (sph%vs(i_h)%edge4%edge_primal_arc - len0)*ej4%x
        IF (sph%vs(i_h)%edge5%info%idx /= -1 ) THEN
          feder%x = feder%x + &
            & (sph%vs(i_h)%edge5%edge_primal_arc - len0)*ej5%x
        ENDIF

        ! RHS of spring equation
        rhs%x = feder%x-friction*sph%vs(i_h)%vert_vel%x

        ! solve the position equation
        sph%vs(i_h)%vertex%x = sph%vs(i_h)%vertex%x+dt*sph%vs(i_h)%vert_vel%x
        ! - > Normalize
        sph%vs(i_h)%vertex%x = sph%vs(i_h)%vertex%x/&
          & SQRT(DOT_PRODUCT(sph%vs(i_h)%vertex%x,sph%vs(i_h)%vertex%x))

        ! solve for spring equation
        sph%vs(i_h)%vert_vel%x = sph%vs(i_h)%vert_vel%x + dt*rhs%x
        hori = DOT_PRODUCT(sph%vs(i_h)%vert_vel%x,sph%vs(i_h)%vertex%x)
        ! - > Horizontalize
        sph%vs(i_h)%vert_vel%x = sph%vs(i_h)%vert_vel%x - &
          & hori * sph%vs(i_h)%vertex%x

        ! kinetic energy for testing
        ekin = 0.5_wp*dot_product(sph%vs(i_h)%vert_vel%x,&
          & sph%vs(i_h)%vert_vel%x) + ekin

        ! test criterion
        test = DOT_PRODUCT(feder%x,feder%x) + test
     ENDDO
#endif

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


