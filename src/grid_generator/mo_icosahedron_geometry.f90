!>
!!               The subroutines in this module compute.
!!
!!               The subroutines in this module compute
!!               and store all the geometric information
!!               for the original icosahedron (level -1 of the grid hierarchy).
!!
!! @par Revision History
!! Initial version  by:
!! @par
!! Luis Kornblueh,  MPI-M, Hamburg, January 2004
!! @par
!! Luis Kornblueh, MPI-M, Hamburg, February 2004
!!          - change to compiling and running version
!!          - create level 0 subdivision, connect triangle pointers
!! Luis Kornblueh, MPI-M, Hamburg, March 2004
!!          - include parents
!!          - full working version of tree and graph representation
!! Luca Bonaventura,  MPI-M, Hamburg, October 2004
!!          - include documentation and Protex headers
!! Almut Gassmann, MPI-M (200-01-15)
!! - intruduce pseudo icosahedron on a double periodic plane
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
MODULE mo_icosahedron_geometry
  !
  !
  !

  USE mo_kind,               ONLY: wp
  USE mo_math_constants,     ONLY: pi_5
  USE mo_math_types,         ONLY: t_cartesian_coordinates
  USE mo_math_utilities,     ONLY: rotate_x, rotate_y, rotate_z
  USE mo_base_geometry,      ONLY: x_rot_angle,y_rot_angle,z_rot_angle

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_icosahedron_vertices
  PUBLIC :: init_planar_vertices
  PUBLIC :: get_icosahedron_vertex
  
  TYPE t_coordinates_reference
    TYPE(t_cartesian_coordinates), POINTER :: px
  END TYPE t_coordinates_reference

  TYPE(t_coordinates_reference)          :: icosahedron_vertices(0:2,0:19)
  TYPE(t_cartesian_coordinates), TARGET :: base_vertex(12)



CONTAINS
  !-------------------------------------------------------------------------
  !BOC
  !EOC
  !-------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: init_icosahedron_vertices
  !
  ! !SUBROUTINE INTERFACE:
  SUBROUTINE init_icosahedron_vertices
    !
    ! !DESCRIPTION:
    ! Calculates the Cartesian coordinates of the gridpoints of the
    ! icosahedral grid on the unit sphere. The vertices are saved
    ! in an array from 1 till 12 for comparison with other published
    ! solutions. Computes the  angles associated with the base icosahedron.
    ! Calculate the latitudes of the two times 5 vertices
    ! (the other two important vertices are south and north pole).
    !
    ! !REVISION HISTORY:
    !  Luis Kornblueh, MPI-M, Hamburg, March 2004
    ! !USES:


    REAL(wp) :: z_w, z_cosw, z_sinw, z_rlon
    TYPE(t_cartesian_coordinates) :: temp
    INTEGER :: i_mdist(2:11) ! distribute vertices on latitude ring

    INTEGER :: i_msgn, j

    !EOP
    !-------------------------------------------------
    !BOC

    z_w      = 2.0_wp*acos(1.0_wp/(2.0_wp*sin(pi_5)))
    z_cosw   = COS(z_w)
    z_sinw   = SIN(z_w)

    ! Northern hemisphere are the first 6 elements of base_vertex(1:6)
    ! Southern hemisphere are the other 6 elements of base_vertex(7:12)
    !
    ! set poles first - it is simple

    base_vertex( 1) = t_cartesian_coordinates( (/ 0.0_wp, 0.0_wp,  1.0_wp /) )
    base_vertex(12) = t_cartesian_coordinates( (/ 0.0_wp, 0.0_wp, -1.0_wp /) )

    ! now set the vertices on the two latitude rings

    DO j = 1, 10
      IF (MOD(j,2) == 0) THEN
        i_mdist(j/2+6)     = -1+(j-1)-10*((j-1)/7)
      ELSE
        i_mdist((j+1)/2+1) = -1+(j-1)-10*((j-1)/7)
      END IF
    END DO

    DO j = 2,11

      ! toggle the hemisphere

      IF (j > 6) THEN
        i_msgn = -1    ! southern
      ELSE
        i_msgn =  1    ! northern
      ENDIF

      ! compute the meridian angle for the base vertex.

      z_rlon =  (1.0_wp+REAL(i_mdist(j),wp))*pi_5

      ! now initialize the coordinates

      base_vertex(j) = t_cartesian_coordinates( &
        & (/ z_sinw*cos(z_rlon), z_sinw*sin(z_rlon), z_cosw*REAL(i_msgn,wp) /) &
        & )
    END DO

    !
    !  Perform rotation(s) of original icosahedron
    !

    !
    !  around x axis
    !
    DO j=1,12

      temp=base_vertex(j)
      base_vertex(j)=rotate_x(x_rot_angle,temp)

    ENDDO

    !
    !  around y axis
    !
    DO j=1,12

      temp=base_vertex(j)
      base_vertex(j)=rotate_y(y_rot_angle,temp)

    ENDDO

    !
    !  around z axis
    !
    DO j=1,12

      temp=base_vertex(j)
      base_vertex(j)=rotate_z(z_rot_angle,temp)

    ENDDO

    !LB
    !LB   commented to avoid unused output
    !LB   can be reinserted for debugging purposes
    !LB   or removed entirely
    !LB
    !!$    OPEN (nin,FILE='icosahedron-vertices.dat')
    !!$    DO j = 1, 12
    !!$       position = cc2gc(base_vertex(j))
    !!$       WRITE (nin,'(5x,2f8.1,a,i0)') &
    !!$            rad2deg*position%lon, rad2deg*position%lat, ' 14 0 0 BL ', j
    !!$    ENDDO
    !!$    CLOSE (nin)

    ! upward looking north pole

    icosahedron_vertices(0,0)%px => base_vertex(1)
    icosahedron_vertices(1,0)%px => base_vertex(2)
    icosahedron_vertices(2,0)%px => base_vertex(3)

    icosahedron_vertices(0,1)%px => base_vertex(1)
    icosahedron_vertices(1,1)%px => base_vertex(3)
    icosahedron_vertices(2,1)%px => base_vertex(4)

    icosahedron_vertices(0,2)%px => base_vertex(1)
    icosahedron_vertices(1,2)%px => base_vertex(4)
    icosahedron_vertices(2,2)%px => base_vertex(5)

    icosahedron_vertices(0,3)%px => base_vertex(1)
    icosahedron_vertices(1,3)%px => base_vertex(5)
    icosahedron_vertices(2,3)%px => base_vertex(6)

    icosahedron_vertices(0,4)%px => base_vertex(1)
    icosahedron_vertices(1,4)%px => base_vertex(6)
    icosahedron_vertices(2,4)%px => base_vertex(2)

    ! downward looking main row

    icosahedron_vertices(0,5)%px => base_vertex(7)
    icosahedron_vertices(1,5)%px => base_vertex(3)
    icosahedron_vertices(2,5)%px => base_vertex(2)

    icosahedron_vertices(0,6)%px => base_vertex(8)
    icosahedron_vertices(1,6)%px => base_vertex(4)
    icosahedron_vertices(2,6)%px => base_vertex(3)

    icosahedron_vertices(0,7)%px => base_vertex(9)
    icosahedron_vertices(1,7)%px => base_vertex(5)
    icosahedron_vertices(2,7)%px => base_vertex(4)

    icosahedron_vertices(0,8)%px => base_vertex(10)
    icosahedron_vertices(1,8)%px => base_vertex(6)
    icosahedron_vertices(2,8)%px => base_vertex(5)

    icosahedron_vertices(0,9)%px => base_vertex(11)
    icosahedron_vertices(1,9)%px => base_vertex(2)
    icosahedron_vertices(2,9)%px => base_vertex(6)

    ! upward looking main row

    icosahedron_vertices(0,10)%px => base_vertex(3)
    icosahedron_vertices(1,10)%px => base_vertex(7)
    icosahedron_vertices(2,10)%px => base_vertex(8)

    icosahedron_vertices(0,11)%px => base_vertex(4)
    icosahedron_vertices(1,11)%px => base_vertex(8)
    icosahedron_vertices(2,11)%px => base_vertex(9)

    icosahedron_vertices(0,12)%px => base_vertex(5)
    icosahedron_vertices(1,12)%px => base_vertex(9)
    icosahedron_vertices(2,12)%px => base_vertex(10)

    icosahedron_vertices(0,13)%px => base_vertex(6)
    icosahedron_vertices(1,13)%px => base_vertex(10)
    icosahedron_vertices(2,13)%px => base_vertex(11)

    icosahedron_vertices(0,14)%px => base_vertex(2)
    icosahedron_vertices(1,14)%px => base_vertex(11)
    icosahedron_vertices(2,14)%px => base_vertex(7)

    ! downward looking south pole

    icosahedron_vertices(0,15)%px => base_vertex(12)
    icosahedron_vertices(1,15)%px => base_vertex(8)
    icosahedron_vertices(2,15)%px => base_vertex(7)

    icosahedron_vertices(0,16)%px => base_vertex(12)
    icosahedron_vertices(1,16)%px => base_vertex(9)
    icosahedron_vertices(2,16)%px => base_vertex(8)

    icosahedron_vertices(0,17)%px => base_vertex(12)
    icosahedron_vertices(1,17)%px => base_vertex(10)
    icosahedron_vertices(2,17)%px => base_vertex(9)

    icosahedron_vertices(0,18)%px => base_vertex(12)
    icosahedron_vertices(1,18)%px => base_vertex(11)
    icosahedron_vertices(2,18)%px => base_vertex(10)

    icosahedron_vertices(0,19)%px => base_vertex(12)
    icosahedron_vertices(1,19)%px => base_vertex(7)
    icosahedron_vertices(2,19)%px => base_vertex(11)


  END SUBROUTINE init_icosahedron_vertices
  !---------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Calculates the Cartesian coordinates of the gridpoints of the planar grid.
  !!
  !!
  !! @par Revision History
  !! Almut Gassmann, MPI-M (2009-01-09)
  !!
  SUBROUTINE init_planar_vertices
    !

    REAL(wp) :: sqrt3o2

    !-------------------------------------------------

    sqrt3o2 = SQRT(3.0_wp)/2.0_wp

    !There are 9 base vertices, the center with coord. center: number 5
    !      1__2____3
    !     /\  /\  /
    !   4/__5/__\6
    !   /\  /\  /
    ! 7/__8/__\9
    !
    base_vertex(1) = t_cartesian_coordinates( (/ -0.5_wp,  sqrt3o2, 0.0_wp /) )
    base_vertex(2) = t_cartesian_coordinates( (/  0.5_wp,  sqrt3o2, 0.0_wp /) )
    base_vertex(3) = t_cartesian_coordinates( (/  1.5_wp,  sqrt3o2, 0.0_wp /) )
    base_vertex(4) = t_cartesian_coordinates( (/ -1.0_wp,   0.0_wp, 0.0_wp /) )
    base_vertex(5) = t_cartesian_coordinates( (/  0.0_wp,   0.0_wp, 0.0_wp /) )
    base_vertex(6) = t_cartesian_coordinates( (/  1.0_wp,   0.0_wp, 0.0_wp /) )
    base_vertex(7) = t_cartesian_coordinates( (/ -1.5_wp, -sqrt3o2, 0.0_wp /) )
    base_vertex(8) = t_cartesian_coordinates( (/ -0.5_wp, -sqrt3o2, 0.0_wp /) )
    base_vertex(9) = t_cartesian_coordinates( (/  0.5_wp, -sqrt3o2, 0.0_wp /) )

    icosahedron_vertices(0,0)%px => base_vertex(5)
    icosahedron_vertices(1,0)%px => base_vertex(2)
    icosahedron_vertices(2,0)%px => base_vertex(1)

    icosahedron_vertices(0,1)%px => base_vertex(6)
    icosahedron_vertices(1,1)%px => base_vertex(3)
    icosahedron_vertices(2,1)%px => base_vertex(2)

    icosahedron_vertices(0,2)%px => base_vertex(1)
    icosahedron_vertices(1,2)%px => base_vertex(4)
    icosahedron_vertices(2,2)%px => base_vertex(5)

    icosahedron_vertices(0,3)%px => base_vertex(2)
    icosahedron_vertices(1,3)%px => base_vertex(5)
    icosahedron_vertices(2,3)%px => base_vertex(6)

    icosahedron_vertices(0,4)%px => base_vertex(8)
    icosahedron_vertices(1,4)%px => base_vertex(5)
    icosahedron_vertices(2,4)%px => base_vertex(4)

    icosahedron_vertices(0,5)%px => base_vertex(9)
    icosahedron_vertices(1,5)%px => base_vertex(6)
    icosahedron_vertices(2,5)%px => base_vertex(5)

    icosahedron_vertices(0,6)%px => base_vertex(4)
    icosahedron_vertices(1,6)%px => base_vertex(7)
    icosahedron_vertices(2,6)%px => base_vertex(8)

    icosahedron_vertices(0,7)%px => base_vertex(5)
    icosahedron_vertices(1,7)%px => base_vertex(8)
    icosahedron_vertices(2,7)%px => base_vertex(9)

  END SUBROUTINE init_planar_vertices

  !-------------------------------------------------------------------------

  FUNCTION get_icosahedron_vertex (k_vertex_number, k_triangle_number) &
      & result (vertex)

    INTEGER, INTENT(in) :: k_vertex_number, k_triangle_number
    TYPE(t_cartesian_coordinates) :: vertex

    vertex = icosahedron_vertices(k_vertex_number,k_triangle_number)%px

  END FUNCTION get_icosahedron_vertex

END MODULE mo_icosahedron_geometry

!----------------------------------------------










