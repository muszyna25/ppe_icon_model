!>
!! Some utilities which are specific to the transport algorithm.
!!
!! Module contains some functions and procedures which are specifically related
!! to the transport schemes. These subroutines or functions are needed at
!! various places within the transport scheme. Therefore outsourcing these
!! routines protects from possible circular dependencies.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2010-03-04)
!! Modification by Daniel Reinert, DWD (2010-04-23)
!! - implemented generalized Lax-Friedrich flux function
!!   laxfr_upflux_v, which allows to use the same transport
!!   code for pressure and height based vertical coordinate
!!   systems.
!! Modification by Daniel Reinert, DWD (2010-05-17)
!! - added subroutines back_traj_dreg_o1, prep_gauss_quadrature and function
!!   jac which are part of the Gauss-Legendre quadrature applied in the
!!   Miura-scheme.
!! Modification by Daniel Reinert, DWD (2010-10-14)
!! - added subroutine prep_gauss_quadrature_c for integrating a cubic polynomial.
!!   Renamed old prep_gauss_quadrature to prep_gauss_quadrature_q
!! Modification by Daniel Reinert, DWD (2011-04-21)
!! - moved setup_transport to mo_advection_nml
!! Modification by Daniel Reinert, DWD (2013-10-30)
!! - moved divide_flux_area to mo_advection_geometry
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
MODULE mo_advection_traj

  USE mo_kind,                ONLY: wp, vp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_e
  USE mo_impl_constants,      ONLY: min_rledge_int
  USE mo_timer,               ONLY: timer_start, timer_stop, timers_level, timer_back_traj
  USE mo_advection_utils,     ONLY: t_list2D
!!$  USE mo_math_constants,      ONLY: dbl_eps


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: &
    &  version = '$Id$'


  PUBLIC :: btraj
  PUBLIC :: btraj_dreg
  PUBLIC :: btraj_o2


CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! Computation of first order backward trajectories for FFSL transport scheme
  !!
  !! Computes backward trajectories in order to determine the barycenter of the
  !! departure region. Here, a simple first order method is used. Computations are
  !! performed on a plane tangent to the edge midpoint. Coordinate axes point into
  !! the local normal and tangential direction.
  !! Once the barycenter of the departure region is known, the distance vector
  !! between the circumcenter of the upstream cell and the barycenter is computed.
  !! In a final step, this vector is transformed into a rotated coordinate system
  !! which has its origin at the circumcenter. The coordinate axes point to the local
  !! east and local north. Note that this subroutine has specifically been designed
  !! for the MIURA scheme with second order (linear) reconstruction of the subgrid 
  !! distribution.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-03-17)
  !!
  SUBROUTINE btraj( ptr_p, ptr_int, p_vn, p_vt, p_dthalf, p_cell_idx, &
    &               p_cell_blk, p_distv_bary, opt_rlstart, opt_rlend, &
    &               opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(in) ::      &  !< patch on which computation is performed
         &  ptr_p

    TYPE(t_int_state), TARGET, INTENT(in) ::  &  !< pointer to data structure for interpolation
         &  ptr_int

    REAL(wp), INTENT(IN)    ::  &  !< normal component of velocity vector at edge midpoints
         &  p_vn(:,:,:)

    REAL(wp), INTENT(IN)    ::  &  !< tangential component of velocity vector at
         &  p_vt(:,:,:)               !< edge midpoints

    REAL(wp), INTENT(IN)    ::  &  !< $0.5 \Delta t$
         &  p_dthalf

    REAL(vp), INTENT(OUT)   ::  &  !< distance vectors cell center --> barycenter of
         &  p_distv_bary(:,:,:,:)     !< departure region. (geographical coordinates)
                                      !< dim: (nproma,nlev,ptr_p%nblks_e,2)

    INTEGER, INTENT(OUT)  ::    &  !< line and block indices of upwind cell
         &  p_cell_idx(:,:,:)         !< dim: (nproma,nlev,ptr_p%nblks_e)
    INTEGER, INTENT(OUT)  ::    &  !< line and block indices of upwind cell
         &  p_cell_blk(:,:,:)         !< dim: (nproma,nlev,ptr_p%nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
         &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
         &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
         &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
         &  opt_elev

    REAL(wp) :: z_ntdistv_bary(2)      !< cell center --> barycenter in 'normal' and
                                       !< 'tangential' coordinates.

    INTEGER :: je, jk, jb        !< index of edge, vert level, block
    INTEGER :: nlev              !< number of full levels
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart, i_rlend, i_nchdom
    INTEGER :: slev, elev        !< vertical start and end level
    LOGICAL :: lvn_pos

    !-------------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = ptr_p%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    ! number of vertical levels
    nlev   = ptr_p%nlev

    ! number of child domains
    i_nchdom   = MAX(1,ptr_p%n_childdom)

    i_startblk = ptr_p%edges%start_blk(i_rlstart,1)
    i_endblk   = ptr_p%edges%end_blk(i_rlend,i_nchdom)

    !-------------------------------------------------------------------------
    IF (timers_level > 5) THEN
      CALL timer_start(timer_back_traj)
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,z_ntdistv_bary,lvn_pos) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_p, jb, i_startblk, i_endblk,        &
           i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = slev, elev
        DO je = i_startidx, i_endidx

          !
          ! Calculate backward trajectories
          !

          ! position of barycenter in normal direction
          ! pos_barycenter(1) = - p_vn(je,jk,jb) * p_dthalf

          ! position of barycenter in tangential direction
          ! pos_barycenter(2) = - p_vt(je,jk,jb) * p_dthalf

          ! logical auxiliary for MERGE operations: .TRUE. for vn >= 0
          lvn_pos = p_vn(je,jk,jb) * p_dthalf >= 0._wp

          ! If vn > 0 (vn < 0), the upwind cell is cell 1 (cell 2)

          ! line and block indices of neighbor cell with barycenter
          p_cell_idx(je,jk,jb) = &
             &   MERGE(ptr_p%edges%cell_idx(je,jb,1),ptr_p%edges%cell_idx(je,jb,2),lvn_pos)

          p_cell_blk(je,jk,jb) = &
             &   MERGE(ptr_p%edges%cell_blk(je,jb,1),ptr_p%edges%cell_blk(je,jb,2),lvn_pos)


          ! Calculate the distance cell center --> barycenter for the cell,
          ! in which the barycenter is located. The distance vector points
          ! from the cell center to the barycenter.
          z_ntdistv_bary(1) =  - ( p_vn(je,jk,jb) * p_dthalf     &
               & + MERGE(ptr_int%pos_on_tplane_e(je,jb,1,1),        &
               &         ptr_int%pos_on_tplane_e(je,jb,2,1),lvn_pos))

          z_ntdistv_bary(2) =  - ( p_vt(je,jk,jb) * p_dthalf     &
               & + MERGE(ptr_int%pos_on_tplane_e(je,jb,1,2),        &
               &         ptr_int%pos_on_tplane_e(je,jb,2,2),lvn_pos))

          ! In a last step, transform this distance vector into a rotated
          ! geographical coordinate system with its origin at the circumcenter
          ! of the upstream cell. Coordinate axes point to local East and local
          ! North.

          ! component in longitudinal direction
          p_distv_bary(je,jk,jb,1) =                                                        &
               &   z_ntdistv_bary(1)*MERGE(ptr_p%edges%primal_normal_cell(je,jb,1)%v1,         &
               &                           ptr_p%edges%primal_normal_cell(je,jb,2)%v1,lvn_pos) &
               & + z_ntdistv_bary(2)*MERGE(ptr_p%edges%dual_normal_cell(je,jb,1)%v1,           &
               &                           ptr_p%edges%dual_normal_cell(je,jb,2)%v1,lvn_pos)

          ! component in latitudinal direction
          p_distv_bary(je,jk,jb,2) =                                                        &
               &   z_ntdistv_bary(1)*MERGE(ptr_p%edges%primal_normal_cell(je,jb,1)%v2,         &
               &                           ptr_p%edges%primal_normal_cell(je,jb,2)%v2,lvn_pos) &
               & + z_ntdistv_bary(2)*MERGE(ptr_p%edges%dual_normal_cell(je,jb,1)%v2,           &
               &                           ptr_p%edges%dual_normal_cell(je,jb,2)%v2,lvn_pos)

        ENDDO ! loop over edges
      ENDDO   ! loop over vertical levels
    END DO    ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (timers_level > 5) CALL timer_stop(timer_back_traj)

  END SUBROUTINE btraj



  !-------------------------------------------------------------------------
  !>
  !! Computation of first order backward trajectories for FFSL transport scheme
  !!
  !! Computes backward trajectories in order to determine an approximation to the
  !! departure region. Here, the departure region is approximated by a rhomboid,
  !! using a simple first order accurate backward trajectory. Computations
  !! are performed on a plane tangent to the edge midpoint. Coordinate axes point into
  !! the local normal and tangential direction.
  !! Once the 4 vertices of the departure region are known, the distance vector
  !! between the circumcenter of the upstream cell and the vertices is computed.
  !! In a final step, these vectors are transformed into a rotated coordinate system
  !! which has its origin at the circumcenter. The coordinate axes point to the local
  !! east and local north. This subroutine may be combined with any reconstruction method
  !! for the subgrid distribution.
  !! Note: So far, we take care that the vertices of the departure region are stored in
  !! counterclockwise order. This ensures that the gaussian quadrature is positive definite.
  !!
  !! NOTE_1: Since we are only interested in the departure region average rather than 
  !!       the departure region integral, counterclockwise numbering is not strictly 
  !!       necessary. Maybe we should remove the computational overhead of counterclockwise 
  !!       numbering at some time. However, the points must not be numbered in random order.
  !!       Care must be taken that the points are numbered either clockwise or counterclockwise. 
  !!       
  !! Note_2: The coordinates for 2 of the 4 vertices do not change with time. However, 
  !!       tests indicated that re-computing these coordinates is faster than fetching 
  !!       precomputed ones from memory. 
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-05-12)
  !! Modification by Daniel Reinert, DWD (2012-04-24)
  !! - bug fix: counterclockwise numbering is now ensured independent of the 
  !!   edge system-orientation.
  !! Modification by Daniel Reinert, DWD (2013-11-01)
  !! - optionally derive list of edges for which the standard Miura scheme is 
  !!   potentially insufficient
  !!
  SUBROUTINE btraj_dreg( ptr_p, ptr_int, p_vn, p_vt, p_dt, lcounterclock, &
       &                   p_cell_idx, p_cell_blk, p_coords_dreg_v,       &
       &                   opt_rlstart, opt_rlend, opt_slev, opt_elev,    &
       &                   opt_falist )

    TYPE(t_patch), TARGET, INTENT(IN) ::     &  !< patch on which computation is performed
         &  ptr_p

    TYPE(t_int_state), TARGET, INTENT(IN) :: &  !< pointer to data structure for interpolation
         &  ptr_int

    REAL(wp), INTENT(IN)    ::  &  !< normal component of velocity vector at edge midpoints
         &  p_vn(:,:,:)

    REAL(wp), INTENT(IN)    ::  &  !< tangential component of velocity vector at
         &  p_vt(:,:,:)               !< edge midpoints

    REAL(wp), INTENT(IN)    ::  &  !< time step $\Delta t$
         &  p_dt

    LOGICAL, INTENT(IN)     ::  &  !< if TRUE, flux area vertices are ordered 
         &  lcounterclock            !< counterclockwise. If FALSE, some are ordered
                                     !< counterclockwise, some clockwise

    REAL(vp), INTENT(OUT) ::    &  !< coordinates of departure region vertices. The origin
         &  p_coords_dreg_v(:,:,:,:,:)!< of the coordinate system is at the circumcenter of
                                      !< the upwind cell. Unit vectors point to local East
                                      !< and North. (geographical coordinates)
                                      !< dim: (nproma,4,2,nlev,ptr_p%nblks_e)

    INTEGER, INTENT(OUT)  ::    &  !< line and block indices of upwind cell
         &  p_cell_idx(:,:,:)   !< dim: (nproma,nlev,ptr_p%nblks_e)
    INTEGER, INTENT(OUT)  ::    &  !< line and block indices of upwind cell
         &  p_cell_blk(:,:,:)   !< dim: (nproma,nlev,ptr_p%nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
         &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
         &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
         &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
         &  opt_elev

    TYPE(t_list2D), INTENT(INOUT), OPTIONAL :: & !< list with points for which a local
         &  opt_falist                      !< polynomial approximation is insufficient 
                                            !< and a piecewise approximation is needed, 
                                            !< instead

    REAL(wp) ::            &       !< coordinates of departure points 
         &  depart_pts(nproma,2,2) !< in edge-based coordinate system

    REAL(wp) ::            &       !< coordinates of departure region vertices
         &  pos_dreg_vert_c(nproma,4,2)      !< as seen from translated coordinate system.
                                      !< origin at circumcenter of upwind cell

    REAL(wp) ::            &       !< position on tangential plane depending
         &  pos_on_tplane_e(nproma,2)        !< on the sign of vn

    REAL(wp) ::            &       !< primal and dual normals of cell lying
         &  pn_cell(2), dn_cell(2)    !< in the direction of vn

    REAL(wp) ::            &       !< edge vertices
         &  edge_verts(nproma,2,2)

    INTEGER :: je, jk, jb          !< index of edge, vert level, block
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart, i_rlend, i_nchdom
    INTEGER :: slev, elev          !< vertical start and end level
    LOGICAL :: lvn_pos             !< vn > 0: TRUE/FALSE
    LOGICAL :: lvn_sys_pos(nproma,ptr_p%nlev)   !< vn*system_orient > 0

    ! for index list generation
    LOGICAL :: llist_gen           !< if TRUE, generate index list
    INTEGER :: ie                  !< counter
    REAL(wp):: traj_length         !< backward trajectory length [m]
    REAL(wp):: e2c_length          !< edge-upwind cell circumcenter length [m]
    !-------------------------------------------------------------------------

    IF (timers_level > 5) THEN
      CALL timer_start(timer_back_traj)
    ENDIF

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = ptr_p%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    IF (PRESENT(opt_falist)) THEN
      llist_gen = .TRUE.
    ELSE
      llist_gen = .FALSE.
    ENDIF

    ! number of child domains
    i_nchdom   = MAX(1,ptr_p%n_childdom)

    i_startblk = ptr_p%edges%start_blk(i_rlstart,1)
    i_endblk   = ptr_p%edges%end_blk(i_rlend,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,ie,i_startidx,i_endidx,traj_length,e2c_length, &
!$OMP depart_pts,pos_dreg_vert_c,pos_on_tplane_e,pn_cell,dn_cell,lvn_pos,&
!$OMP lvn_sys_pos,edge_verts) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_p, jb, i_startblk, i_endblk, &
           i_startidx, i_endidx, i_rlstart, i_rlend)


      ! get local copy of edge vertices
      DO je = i_startidx, i_endidx
        edge_verts(je,1:2,1:2) = ptr_int%pos_on_tplane_e(je,jb,7:8,1:2)
      ENDDO

      ! logical switch for merge options regarding the counterclockwise numbering
      IF (lcounterclock) THEN
        DO jk = slev, elev
          DO je = i_startidx, i_endidx
            lvn_sys_pos(je,jk) = p_vn(je,jk,jb)*ptr_p%edges%system_orientation(je,jb) >= 0._wp
          ENDDO
        ENDDO
      ELSE
        lvn_sys_pos(i_startidx:i_endidx,slev:elev) = .FALSE.
      ENDIF

      ! generate list of points that require special treatment
      !
      IF (llist_gen) THEN
        ie = 0
        DO jk = slev, elev
          DO je = i_startidx, i_endidx
            ! logical switch for MERGE operations: .TRUE. for p_vn >= 0
            lvn_pos     = p_vn(je,jk,jb) * p_dt >= 0._wp

            ! compute length of backward trajectory
            traj_length = SQRT(p_vn(je,jk,jb)**2 + p_vt(je,jk,jb)**2) * p_dt
!!$            traj_length = ABS(p_vn(je,jk,jb)) * p_dt

            ! distance from edge midpoint to upwind cell circumcenter [m]
            e2c_length  = MERGE(ptr_p%edges%edge_cell_length(je,jb,1),       &
              &                 ptr_p%edges%edge_cell_length(je,jb,2),lvn_pos)

            IF (traj_length > 1.25_wp*e2c_length) THEN   ! add point to index list
!!$            IF (traj_length > ((1.4_wp - MIN(0.5_wp,(0.1_wp*ABS(p_vt(je,jk,jb)/MAX(dbl_eps,p_vn(je,jk,jb)))))) &
!!$                               *e2c_length)) THEN   ! add point to index list

              ie = ie + 1
              opt_falist%eidx(ie,jb) = je
              opt_falist%elev(ie,jb) = jk
            ENDIF
          ENDDO ! loop over edges
        ENDDO   ! loop over vertical levels

        ! store list dimension
        opt_falist%len(jb) = ie
      ENDIF

      DO jk = slev, elev
!DIR$ IVDEP
        DO je = i_startidx, i_endidx



          ! departure region and correct counterclockwise numbering of vertices
          !--------------------------------------------------------------------
          !
          !        -1                          +1            : system orientation
          !
          !  3\--------------\2       4\--------------\3
          !    \              \         \              \         [vn < 0]
          !     \      N       \         \      N       \
          !      \      \       \         \      \       \      <- edge normal
          !   4,2 \------\-------\1      1 \------\-------\2,4  <- triangle edge
          !        \              \         \              \
          !         \              \         \              \
          !          \              \         \              \   [vn > 0]
          !         3 \--------------\4      2 \--------------\3


          ! Determine the upwind cell
          ! Cell indices are chosen such that the direction from cell 1 to cell 2
          ! is the positive direction of the normal vector N.
          ! Edge indices are chosen such that the direction from edge 1 to edge 2
          ! is the positive directions of the tangential vector T.
          !
          ! If (T,N,Z) form a right-hand system, the system orientation is 1.
          ! If (T,N,Z) form a left-hand system, the system orientation is -1.

          !
          ! Calculate backward trajectories, starting at the two edge vertices
          ! (arrival points). It is assumed that the velocity vector is constant 
          ! along the edge.
          !

          ! logical switch for MERGE operations: .TRUE. for p_vn >= 0
          lvn_pos     = p_vn(je,jk,jb) * p_dt >= 0._wp


          ! get line and block indices of upwind cell
          p_cell_idx(je,jk,jb) = MERGE(ptr_p%edges%cell_idx(je,jb,1),       &
            &                          ptr_p%edges%cell_idx(je,jb,2),lvn_pos)
          p_cell_blk(je,jk,jb) = MERGE(ptr_p%edges%cell_blk(je,jb,1),       &
            &                          ptr_p%edges%cell_blk(je,jb,2),lvn_pos)


          ! departure points of the departure cell. Point 1 belongs to edge-vertex 1, 
          ! point 2 belongs to edge_vertex 2.
          !
          ! position of vertex 4 (vn > 0) / vertex 2(vn < 0) in normal direction
          depart_pts(je,1,1)      = edge_verts(je,1,1) - p_vn(je,jk,jb) * p_dt

          ! position of vertex 4 (vn > 0) / vertex 2(vn < 0) in tangential direction
          depart_pts(je,1,2)      = edge_verts(je,1,2) - p_vt(je,jk,jb) * p_dt

          ! position of vertex 3 in normal direction
          depart_pts(je,2,1)      = edge_verts(je,2,1) - p_vn(je,jk,jb) * p_dt

          ! position of vertex 3 in tangential direction
          depart_pts(je,2,2)      = edge_verts(je,2,2) - p_vt(je,jk,jb) * p_dt



          ! determine correct position on tangential plane
          pos_on_tplane_e(je,1) = MERGE(ptr_int%pos_on_tplane_e(je,jb,1,1), &
               &                        ptr_int%pos_on_tplane_e(je,jb,2,1),lvn_pos)

          pos_on_tplane_e(je,2) = MERGE(ptr_int%pos_on_tplane_e(je,jb,1,2), &
               &                        ptr_int%pos_on_tplane_e(je,jb,2,2),lvn_pos)

          ! Calculate position of departure region vertices in a translated
          ! coordinate system. The origin is located at the circumcenter
          ! of the upwind cell. The distance vectors point from the cell center
          ! to the vertices.
          !
          ! Take care of correct counterclockwise numbering below
          !
          pos_dreg_vert_c(je,1,1:2) = edge_verts(je,1,1:2) - pos_on_tplane_e(je,1:2)

          pos_dreg_vert_c(je,2,1)   = MERGE(depart_pts(je,1,1),edge_verts(je,2,1),lvn_sys_pos(je,jk)) &
               &                    - pos_on_tplane_e(je,1)
          pos_dreg_vert_c(je,2,2)   = MERGE(depart_pts(je,1,2),edge_verts(je,2,2),lvn_sys_pos(je,jk)) &
               &                    - pos_on_tplane_e(je,2)

          pos_dreg_vert_c(je,3,1:2) = depart_pts(je,2,1:2) - pos_on_tplane_e(je,1:2)

          pos_dreg_vert_c(je,4,1)   = MERGE(edge_verts(je,2,1),depart_pts(je,1,1),lvn_sys_pos(je,jk)) &
               &                    - pos_on_tplane_e(je,1)
          pos_dreg_vert_c(je,4,2)   = MERGE(edge_verts(je,2,2),depart_pts(je,1,2),lvn_sys_pos(je,jk)) &
               &                    - pos_on_tplane_e(je,2)



          ! In a last step, these distance vectors are transformed into a rotated
          ! geographical coordinate system, which still has its origin at the circumcenter
          ! of the upwind cell. Now the coordinate axes point to local East and local
          ! North.
          !
          ! Determine primal and dual normals of the cell lying in the direction of vn
          pn_cell(1) = MERGE(ptr_p%edges%primal_normal_cell(je,jb,1)%v1,       &
               &             ptr_p%edges%primal_normal_cell(je,jb,2)%v1,lvn_pos)

          pn_cell(2) = MERGE(ptr_p%edges%primal_normal_cell(je,jb,1)%v2,       &
               &             ptr_p%edges%primal_normal_cell(je,jb,2)%v2,lvn_pos)

          dn_cell(1) = MERGE(ptr_p%edges%dual_normal_cell(je,jb,1)%v1,       &
               &             ptr_p%edges%dual_normal_cell(je,jb,2)%v1,lvn_pos)

          dn_cell(2) = MERGE(ptr_p%edges%dual_normal_cell(je,jb,1)%v2,       &
               &             ptr_p%edges%dual_normal_cell(je,jb,2)%v2,lvn_pos)

          ! components in longitudinal direction
          p_coords_dreg_v(je,1:4,1,jk,jb) =                                         &
               & pos_dreg_vert_c(je,1:4,1) * pn_cell(1) + pos_dreg_vert_c(je,1:4,2) * dn_cell(1)

          ! components in latitudinal direction
          p_coords_dreg_v(je,1:4,2,jk,jb) =                                         &
               & pos_dreg_vert_c(je,1:4,1) * pn_cell(2) + pos_dreg_vert_c(je,1:4,2) * dn_cell(2)

        ENDDO ! loop over edges
      ENDDO   ! loop over vertical levels
    END DO    ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (timers_level > 5) CALL timer_stop(timer_back_traj)

  END SUBROUTINE btraj_dreg



  !-------------------------------------------------------------------------
  !>
  !! Computation of second order backward trajectories for FFSL transport scheme
  !!
  !! Computes backward trajectories in order to determine the barycenter of the
  !! departure region. Here, an iterative second order method is used. Computations
  !! are performed on a plane tangent to the edge midpoint. Coordinate axes point
  !! into the local normal and tangential direction. A bilinear interpolation is
  !! applied using the velocity vectors at the three edge midpoints of the upwind cell
  !! in order to derive an improved estimate of the velocity.
  !!
  !! Once the barycenter of the departure region is known, the distance vector
  !! between the circumcenter of the upstream cell and the barycenter is computed.
  !! In a final step, this vector is transformed into a different coordinate system
  !! which has its origin at the circumcenter. The coordinate axes point to the local
  !! east and local north.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-03-24)
  !!
  SUBROUTINE btraj_o2( ptr_p, ptr_int, p_vn, p_vt, p_dthalf, p_cell_idx, &
       &               p_cell_blk, p_distv_bary, opt_rlstart, opt_rlend, &
       &               opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::      &  !< patch on which computation is performed
         &  ptr_p

    TYPE(t_int_state), TARGET, INTENT(IN) ::  &  !< pointer to data structure for interpolation
         &  ptr_int

    REAL(wp), INTENT(IN)    ::  &  !< normal component of velocity vector at edge midpoints
         &  p_vn(:,:,:)

    REAL(wp), INTENT(IN)    ::  &  !< tangential component of velocity vector at
         &  p_vt(:,:,:)               !< edge midpoints

    REAL(wp), INTENT(IN)  ::    &  !< $0.5 \Delta t$
         &  p_dthalf

    REAL(vp), INTENT(OUT) ::    &  !< distance vectors cell center --> barycenter of
         &  p_distv_bary(:,:,:,:)     !< departure region. (geographical coordinates)
                                      !< dim: (nproma,nlev,ptr_p%nblks_e,2)

    INTEGER, INTENT(OUT)  ::    &  !< line indices of cell centers in which the
         &  p_cell_idx(:,:,:)         !< computed barycenters are located.
                                      !< dim: (nproma,nlev,ptr_p%nblks_e)

    INTEGER, INTENT(OUT)  ::    &  !< block indices of cell centers in which the
         &  p_cell_blk(:,:,:)         !< computed barycenters are located.
                                      !< dim: (nproma,nlev,ptr_p%nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
         &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
         &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
         &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
         &  opt_elev

    REAL(wp) :: pos_barycenter(2),  &   !< position of barycenter and distance vector
         &         z_ntdistv_bary(2)       !< cell center --> barycenter in 'normal' and
    !< 'tangential' coordinates.

    REAL(wp), POINTER  ::  &     !< pointer to positions of quadrilateral edge midpoints
         &  ptr_em(:,:,:,:)         !< on tangent plane

    REAL(wp)           ::  &     ! normal and tangential components of the velocity vectors
         &  z_vn_plane(nproma,ptr_p%nlev,ptr_p%nblks_e,4), &!< projected into local edge-based
         &  z_vt_plane(nproma,ptr_p%nlev,ptr_p%nblks_e,4)   !< coordinate system


    REAL(wp)  :: w1, w2, w3      !< weights for bilinear interpolation

    REAL(wp)  ::           &     !< resulting normal and tangential velocity component
         &  vn_new, vt_new          !< after bilinear interpolation

    INTEGER :: je, jk, jb        !< index of edge, vert level, block
    INTEGER :: nlev              !< number of full levels
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart, i_rlend, i_nchdom
    INTEGER :: slev, elev        !< vertical start and end level
    INTEGER :: zcell             !< determines whether the barycenter is located
    !< in cell 1 or 2
    INTEGER, POINTER ::    &     !< pointer for line and block indices of edge
         & iidx(:,:,:), iblk(:,:,:) !< midpoints for quadrilateral cell

    !DR    REAL(wp) :: z_vabs_orig, z_vabs_new

    !-------------------------------------------------------------------------


    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = ptr_p%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 1
    ENDIF

    ! number of vertical levels
    nlev = ptr_p%nlev

    ! number of child domains
    i_nchdom   = MAX(1,ptr_p%n_childdom)

    i_startblk = ptr_p%edges%start_blk(i_rlstart,1)
    i_endblk   = ptr_p%edges%end_blk(i_rlend,i_nchdom)


    ! For each edge the velocity vectors at the quadrilateral edges are
    ! transformed into the new edge-based coordinate system (coordinate axes
    ! normal and tangential to the inner edge).
    ! - projection of normal component into local normal direction
    ! - projection of tangential component into local normal direction
    ! - projection of normal component into local tangential direction
    ! - projection of tangential component into local tangential direction

    ! pointer to the quadrilateral edge midpoints
    ! Note that ptr_em(:,:,1:2,:) belong to cell number 1 and
    ! ptr_em(:,:,3:4,:) belong to cell number 2 as seen from the edge.
    ptr_em => ptr_int%pos_on_tplane_e(:,:,3:6,:)

    ! pointer to line and block indices of outer quadrilateral edge midpoints
    iidx => ptr_p%edges%quad_idx
    iblk => ptr_p%edges%quad_blk

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_p, jb, i_startblk, i_endblk,        &
           i_startidx, i_endidx, i_rlstart, i_rlend)

#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
          DO jk = slev, elev
            DO je = i_startidx, i_endidx
#endif

              z_vn_plane(je,jk,jb,1) =                                                         &
                   &   p_vn(iidx(je,jb,1),jk,iblk(je,jb,1)) * ptr_int%tplane_e_dotprod(je,jb,1,1) &
                   & + p_vt(iidx(je,jb,1),jk,iblk(je,jb,1)) * ptr_int%tplane_e_dotprod(je,jb,1,2)

              z_vn_plane(je,jk,jb,2) =                                                         &
                   &   p_vn(iidx(je,jb,2),jk,iblk(je,jb,2)) * ptr_int%tplane_e_dotprod(je,jb,2,1) &
                   & + p_vt(iidx(je,jb,2),jk,iblk(je,jb,2)) * ptr_int%tplane_e_dotprod(je,jb,2,2)

              z_vn_plane(je,jk,jb,3) =                                                         &
                   &   p_vn(iidx(je,jb,3),jk,iblk(je,jb,3)) * ptr_int%tplane_e_dotprod(je,jb,3,1) &
                   & + p_vt(iidx(je,jb,3),jk,iblk(je,jb,3)) * ptr_int%tplane_e_dotprod(je,jb,3,2)

              z_vn_plane(je,jk,jb,4) =                                                         &
                   &   p_vn(iidx(je,jb,4),jk,iblk(je,jb,4)) * ptr_int%tplane_e_dotprod(je,jb,4,1) &
                   & + p_vt(iidx(je,jb,4),jk,iblk(je,jb,4)) * ptr_int%tplane_e_dotprod(je,jb,4,2)


              z_vt_plane(je,jk,jb,1) =                                                         &
                   &   p_vn(iidx(je,jb,1),jk,iblk(je,jb,1)) * ptr_int%tplane_e_dotprod(je,jb,1,3) &
                   & + p_vt(iidx(je,jb,1),jk,iblk(je,jb,1)) * ptr_int%tplane_e_dotprod(je,jb,1,4)

              z_vt_plane(je,jk,jb,2) =                                                         &
                   &   p_vn(iidx(je,jb,2),jk,iblk(je,jb,2)) * ptr_int%tplane_e_dotprod(je,jb,2,3) &
                   & + p_vt(iidx(je,jb,2),jk,iblk(je,jb,2)) * ptr_int%tplane_e_dotprod(je,jb,2,4)

              z_vt_plane(je,jk,jb,3) =                                                         &
                   &   p_vn(iidx(je,jb,3),jk,iblk(je,jb,3)) * ptr_int%tplane_e_dotprod(je,jb,3,3) &
                   & + p_vt(iidx(je,jb,3),jk,iblk(je,jb,3)) * ptr_int%tplane_e_dotprod(je,jb,3,4)

              z_vt_plane(je,jk,jb,4) =                                                         &
                   &   p_vn(iidx(je,jb,4),jk,iblk(je,jb,4)) * ptr_int%tplane_e_dotprod(je,jb,4,3) &
                   & + p_vt(iidx(je,jb,4),jk,iblk(je,jb,4)) * ptr_int%tplane_e_dotprod(je,jb,4,4)


              !            ! Re-normalize the projected velocity vector
              !            z_vabs_orig = SQRT( p_vn(ilq,jk,ibq)*p_vn(ilq,jk,ibq)                &
              !              &         + p_vt(ilq,jk,ibq)*p_vt(ilq,jk,ibq) )
              !            z_vabs_new  = SQRT( z_vn_plane(je,jk,jb,ne)*z_vn_plane(je,jk,jb,ne)  &
              !              &         + z_vt_plane(je,jk,jb,ne)*z_vt_plane(je,jk,jb,ne) )

              !            z_vn_plane(je,jk,jb,ne) = z_vn_plane(je,jk,jb,ne) * (z_vabs_orig/z_vabs_new)
              !            z_vt_plane(je,jk,jb,ne) = z_vt_plane(je,jk,jb,ne) * (z_vabs_orig/z_vabs_new)


            ENDDO ! loop over edges
          ENDDO   ! loop over vertical levels
        END DO    ! loop over blocks
!$OMP END DO



!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,pos_barycenter,zcell,w1,w2,w3, &
!$OMP            vn_new,vt_new,z_ntdistv_bary) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(ptr_p, jb, i_startblk, i_endblk,        &
               i_startidx, i_endidx, i_rlstart, i_rlend)

          DO jk = slev, elev
            DO je = i_startidx, i_endidx


              !
              ! Calculate backward trajectories
              !

              ! First guess
              ! position of barycenter in normal direction
              pos_barycenter(1) = - p_vn(je,jk,jb) * p_dthalf

              ! position of barycenter in tangential direction
              pos_barycenter(2) = - p_vt(je,jk,jb) * p_dthalf



              IF (p_vn(je,jk,jb) * p_dthalf >= 0._wp) THEN

                !! we are in cell 1 !!

                ! line and block indices of neighboring cell with barycenter
                p_cell_idx(je,jk,jb) = ptr_p%edges%cell_idx(je,jb,1)
                p_cell_blk(je,jk,jb) = ptr_p%edges%cell_blk(je,jb,1)

                zcell = 1

                ! calculate weights for bilinear interpolation of velocities onto the
                ! barycenter
                w3 = ( pos_barycenter(2)                                          &
                     &   - pos_barycenter(1) * ptr_em(je,jb,1,2)/ptr_em(je,jb,1,1) ) &
                     &   / ( ptr_em(je,jb,2,2)                                       &
                     &   - ptr_em(je,jb,2,1) * ptr_em(je,jb,1,2)/ptr_em(je,jb,1,1) )

                w2 = (pos_barycenter(1) - w3 * ptr_em(je,jb,2,1))       &
                     &  / ptr_em(je,jb,1,1)

                w1 = 1._wp - w2 - w3


                ! calculate improved normal and tangential velocity components using
                ! the bilinear interpolation weights.
                vn_new = w1 * p_vn(je,jk,jb) + w2 * z_vn_plane(je,jk,jb,1) &
                     &    + w3 * z_vn_plane(je,jk,jb,2)
                vt_new = w1 * p_vt(je,jk,jb) + w2 * z_vt_plane(je,jk,jb,1) &
                     &    + w3 * z_vt_plane(je,jk,jb,2)

              ELSE

                !! we are in cell 2 !!

                ! line and block indices of neighboring cell with barycenter
                p_cell_idx(je,jk,jb) = ptr_p%edges%cell_idx(je,jb,2)
                p_cell_blk(je,jk,jb) = ptr_p%edges%cell_blk(je,jb,2)

                zcell = 2


                ! calculate weights for bilinear interpolation of velocities onto the
                ! barycenter
                w3 = ( pos_barycenter(2)                                          &
                     &   - pos_barycenter(1) * ptr_em(je,jb,3,2)/ptr_em(je,jb,3,1) ) &
                     &   / ( ptr_em(je,jb,4,2)                                       &
                     &   - ptr_em(je,jb,4,1) * ptr_em(je,jb,3,2)/ptr_em(je,jb,3,1) )

                w2 = (pos_barycenter(1) - w3 * ptr_em(je,jb,4,1))       &
                     &  / ptr_em(je,jb,3,1)

                w1 = 1._wp - w2 - w3


                ! calculate improved normal and tangential velocity components using
                ! the bilinear interpolation weights.
                vn_new = w1 * p_vn(je,jk,jb) + w2 * z_vn_plane(je,jk,jb,3) &
                     &    + w3 * z_vn_plane(je,jk,jb,4)
                vt_new = w1 * p_vt(je,jk,jb) + w2 * z_vt_plane(je,jk,jb,3) &
                     &    + w3 * z_vt_plane(je,jk,jb,4)

              ENDIF


              ! Improved guess
              ! position of barycenter in normal direction
              pos_barycenter(1) = - vn_new * p_dthalf

              ! position of barycenter in tangential direction
              pos_barycenter(2) = - vt_new * p_dthalf


              ! Calculate the distance cell center --> barycenter for the cell,
              ! in which the barycenter is located. The distance vector points
              ! from the cell center to the barycenter.
              z_ntdistv_bary(1:2) = pos_barycenter(1:2)                      &
                   &                 - ptr_int%pos_on_tplane_e(je,jb,zcell,1:2)


              ! In a last step, transform this distance vector into a rotated
              ! geographical coordinate system with its origin at the circumcenter
              ! of the upstream cell. Coordinate axes point to local East and local
              ! North.

              ! component in longitudinal direction
              p_distv_bary(je,jk,jb,1) =                                                 &
                   &    z_ntdistv_bary(1) * ptr_p%edges%primal_normal_cell(je,jb,zcell)%v1  &
                   &  + z_ntdistv_bary(2) * ptr_p%edges%dual_normal_cell(je,jb,zcell)%v1

              ! component in latitudinal direction
              p_distv_bary(je,jk,jb,2) =                                                 &
                   &    z_ntdistv_bary(1) * ptr_p%edges%primal_normal_cell(je,jb,zcell)%v2  &
                   &  + z_ntdistv_bary(2) * ptr_p%edges%dual_normal_cell(je,jb,zcell)%v2



            ENDDO ! loop over edges
          ENDDO   ! loop over vertical levels
        END DO    ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      END SUBROUTINE btraj_o2

END MODULE mo_advection_traj

