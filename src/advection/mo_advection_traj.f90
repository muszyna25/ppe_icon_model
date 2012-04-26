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
!!   jac which are part of the Gauss-Legendre quadrature apllied in the
!!   Miura-scheme.
!! Modification by Daniel Reinert, DWD (2010-10-14)
!! - added subroutine prep_gauss_quadrature_c for integrating a cubic polynomial.
!!   Renamed old prep_gauss_quadrature to prep_gauss_quadrature_q
!! Modification by Daniel Reinert, DWD (2011-04-21)
!! - moved setup_transport to mo_advection_nml
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

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_e
  USE mo_impl_constants,      ONLY: min_rledge_int, max_char_length
  USE mo_timer,               ONLY: timer_start, timer_stop, timers_level, new_timer
  USE mo_math_utilities,      ONLY: ccw, lintersect, line_intersect, t_line, &
    &                               t_geographical_coordinates
  USE mo_exception,           ONLY: finish, message, message_text


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: btraj
  PUBLIC :: btraj_dreg
  PUBLIC :: btraj_o2
  PUBLIC :: btraj_dreg_nosort
  PUBLIC :: divide_flux_area

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  INTEGER :: timer_back_traj_o1   = 0

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
  !! In a final step, this vector is transformed into a different coordinate system
  !! which has its origin at the circumcenter. The coordinate axes point to the local
  !! east and local north. Note that this subroutine has specifically been designed
  !! for the MIURA scheme with second order reconstruction of the subgrid distribution.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-03-17)
  !!
  SUBROUTINE btraj( ptr_p, ptr_int, p_vn, p_vt, p_dthalf, p_cell_indices, &
    &               p_distv_bary, opt_rlstart, opt_rlend, opt_slev,       &
    &               opt_elev )

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

    REAL(wp), INTENT(OUT)   ::  &  !< distance vectors cell center --> barycenter of
      &  p_distv_bary(:,:,:,:)     !< departure region. (geographical coordinates)
                                   !< dim: (nproma,nlev,ptr_p%nblks_e,2)

    INTEGER, INTENT(OUT)    ::  &  !< line and block indices of cell centers in which the
      &  p_cell_indices(:,:,:,:)   !< calculated barycenters are located.
                                   !< dim: (nproma,nlev,ptr_p%nblks_e,2)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    REAL(wp) :: z_ntdistv_bary(2)       !< cell center --> barycenter in 'normal' and
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
      timer_back_traj_o1 = new_timer("back_traj_o1", timer_back_traj_o1)
      CALL timer_start(timer_back_traj_o1)
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
          lvn_pos = p_vn(je,jk,jb) >= 0._wp

          ! If vn > 0 (vn < 0), the upwind cell is cell 1 (cell 2)

          ! line and block indices of neighbor cell with barycenter
          p_cell_indices(je,jk,jb,1) = &
            & MERGE(ptr_p%edges%cell_idx(je,jb,1),ptr_p%edges%cell_idx(je,jb,2),lvn_pos)

          p_cell_indices(je,jk,jb,2) = &
            & MERGE(ptr_p%edges%cell_blk(je,jb,1),ptr_p%edges%cell_blk(je,jb,2),lvn_pos)


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
    
    IF (timers_level > 5) CALL timer_stop(timer_back_traj_o1)

  END SUBROUTINE btraj


  !-------------------------------------------------------------------------
  !>
  !! Computation of first order backward trajectories for FFSL transport scheme
  !!
  !! Computes backward trajectories in order to determine an approximation to the
  !! departure region. Here, the departure region is approximated by a rhomboidal
  !! region, using a simple first order accurate backward trajectory. Computations
  !! are performed on a plane tangent to the edge midpoint. Coordinate axes point into
  !! the local normal and tangential direction.
  !! Once the 4 vertices of the departure region are known, the distance vector
  !! between the circumcenter of the upstream cell and the vertices is computed.
  !! In a final step, these vectors are transformed into a different coordinate system
  !! which has its origin at the circumcenter. The coordinate axes point to the local
  !! east and local north. This subroutine may be combined with any reconstruction method
  !! for the subgrid distribution.
  !! Note: So far, we take care that the vertices of the departure region are stored in
  !! counterclockwise order. This ensures that the gaussian quadrature is positive definite.
  !!
  !! NOTE: Since we are only interested in the departure region average rather than 
  !!       the departure region integral, counterclockwise numbering is not strictly 
  !!       necessary. Maybe we should remove the computational overhead of counterclockwise 
  !!       numbering at some time. However, the points must not be numbered in random order.
  !!       Care must be taken that the points are numbered either clockwise or counterclockwise. 
  !!       
  !! Note: The coordinates for 2 of the 4 vertices do not change with time. However, 
  !!       tests indicated that re-computing these coordinates is faster than fetching 
  !!       precomputed ones from memory. 
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-05-12)
  !! Modification by Daniel Reinert, DWD (2112-04-24)
  !! - bug fix: counterclockwise numbering is now ensured independent of the 
  !!   edge system-orientation.
  !!
  SUBROUTINE btraj_dreg( ptr_p, ptr_int, p_vn, p_vt, p_dt, p_cell_indices,  &
    &                    p_coords_dreg_v, opt_rlstart, opt_rlend, opt_slev, &
    &                    opt_elev )

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

    REAL(wp), INTENT(OUT) ::    &  !< coordinates of departure region vertices. The origin
      &  p_coords_dreg_v(:,:,:,:,:)!< of the coordinate system is at the circumcenter of
                                   !< the upwind cell. Unit vectors point to local East
                                   !< and North. (geographical coordinates)
                                   !< dim: (nproma,4,2,nlev,ptr_p%nblks_e)

    INTEGER, INTENT(OUT)  ::    &  !< line and block indices of upwind cell
      &  p_cell_indices(:,:,:,:)   !< dim: (nproma,nlev,ptr_p%nblks_e,2)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    REAL(wp), POINTER  ::  &       !< pointer to coordinates of edge vertices
      &  ptr_ve(:,:,:,:)           !< on tangent plane

    REAL(wp) ::            &       !< coordinates of departure points 
      &  depart_pts(2,2)           !< in edge-based coordinate system

    REAL(wp) ::            &       !< coordinates of departure region vertices
      &  pos_dreg_vert_c(4,2)      !< as seen from translated coordinate system.
                                   !< origin at circumcenter of upwind cell

    REAL(wp) ::            &       !< position on tangential plane depending
      &  pos_on_tplane_e(2)        !< on the sign of vn

    REAL(wp) ::            &       !< primal and dual normals of cell lying
      &  pn_cell(2), dn_cell(2)    !< in the direction of vn

    INTEGER :: je, jk, jb          !< index of edge, vert level, block
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart, i_rlend, i_nchdom
    INTEGER :: slev, elev          !< vertical start and end level
    LOGICAL :: lvn_pos             !< vn > 0: TRUE/FALSE
    LOGICAL :: lvn_sys_pos         !< vn*system_orient > 0
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


    ! number of child domains
    i_nchdom   = MAX(1,ptr_p%n_childdom)

    i_startblk = ptr_p%edges%start_blk(i_rlstart,1)
    i_endblk   = ptr_p%edges%end_blk(i_rlend,i_nchdom)

    ! pointer to the edge vertices
    ptr_ve => ptr_int%pos_on_tplane_e(:,:,7:8,:)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,depart_pts,pos_dreg_vert_c, &
!$OMP            pos_on_tplane_e,pn_cell,dn_cell,lvn_pos,lvn_sys_pos) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

     CALL get_indices_e(ptr_p, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = slev, elev
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

          !
          ! Calculate backward trajectories, starting at the two edge vertices
          ! (arrival points). It is assumed that the velocity vector is constant 
          ! along the edge.
          !

          ! logical switch for MERGE operations: .TRUE. for p_vn >= 0
          lvn_pos     = p_vn(je,jk,jb) >= 0._wp

          ! logical switch for merge options regarding the counterclockwise numbering
          lvn_sys_pos = (p_vn(je,jk,jb)*ptr_p%edges%system_orientation(je,jb)) >= 0._wp

          ! get line and block indices of upwind cell
          p_cell_indices(je,jk,jb,1) = MERGE(ptr_p%edges%cell_idx(je,jb,1),       &
            &                                ptr_p%edges%cell_idx(je,jb,2),lvn_pos)
          p_cell_indices(je,jk,jb,2) = MERGE(ptr_p%edges%cell_blk(je,jb,1),       &
            &                                ptr_p%edges%cell_blk(je,jb,2),lvn_pos)


          ! departure points of the departure cell. Point 1 belongs to edge-vertex 1, 
          ! point 2 belongs to edge_vertex 2.
          !
          ! position of vertex 4 (vn > 0) / vertex 2(vn < 0) in normal direction
          depart_pts(1,1)      = ptr_ve(je,jb,1,1) - p_vn(je,jk,jb) * p_dt

          ! position of vertex 4 (vn > 0) / vertex 2(vn < 0) in tangential direction
          depart_pts(1,2)      = ptr_ve(je,jb,1,2) - p_vt(je,jk,jb) * p_dt

          ! position of vertex 3 in normal direction
          depart_pts(2,1)      = ptr_ve(je,jb,2,1) - p_vn(je,jk,jb) * p_dt

          ! position of vertex 3 in tangential direction
          depart_pts(2,2)      = ptr_ve(je,jb,2,2) - p_vt(je,jk,jb) * p_dt



          ! determine correct position on tangential plane
          pos_on_tplane_e(1) = MERGE(ptr_int%pos_on_tplane_e(je,jb,1,1), &
            &                        ptr_int%pos_on_tplane_e(je,jb,2,1),lvn_pos)

          pos_on_tplane_e(2) = MERGE(ptr_int%pos_on_tplane_e(je,jb,1,2), &
            &                        ptr_int%pos_on_tplane_e(je,jb,2,2),lvn_pos)

          ! Calculate position of departure region vertices in a translated
          ! coordinate system. The origin is located at the circumcenter
          ! of the upwind cell. The distance vectors point from the cell center
          ! to the vertices.
          !
          ! Take care of correct counterclockwise numbering below
          !
          pos_dreg_vert_c(1,1:2) = ptr_ve(je,jb,1,1:2) - pos_on_tplane_e(1:2)

          pos_dreg_vert_c(2,1)   = MERGE(depart_pts(1,1),ptr_ve(je,jb,2,1),lvn_sys_pos) &
            &                    - pos_on_tplane_e(1)
          pos_dreg_vert_c(2,2)   = MERGE(depart_pts(1,2),ptr_ve(je,jb,2,2),lvn_sys_pos) &
            &                    - pos_on_tplane_e(2)

          pos_dreg_vert_c(3,1:2) = depart_pts(2,1:2) - pos_on_tplane_e(1:2)

          pos_dreg_vert_c(4,1)   = MERGE(ptr_ve(je,jb,2,1),depart_pts(1,1),lvn_sys_pos) &
            &                    - pos_on_tplane_e(1)
          pos_dreg_vert_c(4,2)   = MERGE(ptr_ve(je,jb,2,2),depart_pts(1,2),lvn_sys_pos) &
            &                    - pos_on_tplane_e(2)



          ! In a last step, these distance vectors are transformed into a rotated
          ! geographical coordinate system, which still has its origin at the circumcenter
          ! of the upwind cell. Now the coordinate axes point to local East and local
          ! North.
          !
          ! Determine primal and dual normals of the cell lying in the direction of vn
          pn_cell(1) = MERGE(ptr_p%edges%primal_normal_cell(je,jb,1)%v1,       &
            &                ptr_p%edges%primal_normal_cell(je,jb,2)%v1,lvn_pos)

          pn_cell(2) = MERGE(ptr_p%edges%primal_normal_cell(je,jb,1)%v2,       &
            &                ptr_p%edges%primal_normal_cell(je,jb,2)%v2,lvn_pos)

          dn_cell(1) = MERGE(ptr_p%edges%dual_normal_cell(je,jb,1)%v1,       &
            &                ptr_p%edges%dual_normal_cell(je,jb,2)%v1,lvn_pos)

          dn_cell(2) = MERGE(ptr_p%edges%dual_normal_cell(je,jb,1)%v2,       &
            &                ptr_p%edges%dual_normal_cell(je,jb,2)%v2,lvn_pos)

          ! components in longitudinal direction
          p_coords_dreg_v(je,1,1,jk,jb) =                                         &
            & pos_dreg_vert_c(1,1) * pn_cell(1) + pos_dreg_vert_c(1,2) * dn_cell(1)

          p_coords_dreg_v(je,2,1,jk,jb) =                                         &
            & pos_dreg_vert_c(2,1) * pn_cell(1) + pos_dreg_vert_c(2,2) * dn_cell(1)

          p_coords_dreg_v(je,3,1,jk,jb) =                                        &
            & pos_dreg_vert_c(3,1) * pn_cell(1) + pos_dreg_vert_c(3,2) * dn_cell(1)

          p_coords_dreg_v(je,4,1,jk,jb) =                                         &
            & pos_dreg_vert_c(4,1) * pn_cell(1) + pos_dreg_vert_c(4,2) * dn_cell(1)


          ! components in latitudinal direction
          p_coords_dreg_v(je,1,2,jk,jb) =                                         &
            & pos_dreg_vert_c(1,1) * pn_cell(2) + pos_dreg_vert_c(1,2) * dn_cell(2)

          p_coords_dreg_v(je,2,2,jk,jb) =                                         &
            & pos_dreg_vert_c(2,1) * pn_cell(2) + pos_dreg_vert_c(2,2) * dn_cell(2)

          p_coords_dreg_v(je,3,2,jk,jb) =                                         &
            & pos_dreg_vert_c(3,1) * pn_cell(2) + pos_dreg_vert_c(3,2) * dn_cell(2)

          p_coords_dreg_v(je,4,2,jk,jb) =                                         &
            & pos_dreg_vert_c(4,1) * pn_cell(2) + pos_dreg_vert_c(4,2) * dn_cell(2)


        ENDDO ! loop over edges
      ENDDO   ! loop over vertical levels
    END DO    ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SUBROUTINE btraj_dreg




  !-------------------------------------------------------------------------
  !>
  !! Computation of first order backward trajectories for FFSL transport scheme
  !!
  !! Computes backward trajectories in order to determine an approximation to the
  !! departure region. Here, the departure region is approximated by a rhomboidal
  !! region, using a simple first order accurate backward trajectory. Computations
  !! are performed on a plane tangent to the edge midpoint. Coordinate axes point into
  !! the local normal and tangential direction.
  !! Once the 4 vertices of the departure region are known, the distance vector
  !! between the circumcenter of the upstream cell and the vertices is computed.
  !! In a final step, these vectors are transformed into a rotated coordinate system
  !! which has its origin at the circumcenter. The coordinate axes point to the local
  !! east and local north. This subroutine may be combined with any reconstruction method
  !! for the subgrid distribution.
  !! Note: The coordinates for 2 of the 4 vertices do not change with time. However, 
  !!       tests indicated that re-computing these coordinates is faster than fetching 
  !!       precomputed ones from memory. 
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-05-12)
  !!
  SUBROUTINE btraj_dreg_nosort( ptr_p, ptr_int, p_vn, p_vt, p_dt,                &
    &                     upwind_cell_idx, upwind_cell_blk, arrival_pts,         &
    &                     depart_pts, opt_rlstart, opt_rlend, opt_slev, opt_elev )

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

    REAL(wp), INTENT(OUT) ::    &  !< coordinates of arrival points. The origin
      &  arrival_pts(:,:,:,:,:)    !< of the coordinate system is at the circumcenter of
                                   !< the upwind cell. Unit vectors point to local East
                                   !< and North. (geographical coordinates)
                                   !< dim: (nproma,2,2,nlev,ptr_p%nblks_e)

    REAL(wp), INTENT(OUT) ::    &  !< coordinates of departure points. The origin
      &  depart_pts(:,:,:,:,:)     !< of the coordinate system is at the circumcenter of
                                   !< the upwind cell. Unit vectors point to local East
                                   !< and North. (geographical coordinates)
                                   !< dim: (nproma,2,2,nlev,ptr_p%nblks_e)

    INTEGER, INTENT(OUT)  ::    &  !< line and block index of upwind cell
      &  upwind_cell_idx(:,:,:),&  !< dim: (nproma,nlev,ptr_p%nblks_e)
      &  upwind_cell_blk(:,:,:)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    REAL(wp), POINTER  ::  &       !< pointer to coordinates of edge vertices
      &  ptr_ve(:,:,:,:)           !< on tangent plane

    REAL(wp) ::            &       !< coordinates of departure points 
      &  dep_pts(2,2)              !< in edge-based coordinate system

    REAL(wp) ::            &       !< coordinates of arrival points 
      &  arr_pts(2,2)              !< in edge-based coordinate system

    REAL(wp) ::            &       !< position on tangential plane depending
      &  pos_on_tplane_e(2)        !< on the sign of vn

    REAL(wp) ::            &       !< primal and dual normals of cell lying
      &  pn_cell(2), dn_cell(2)    !< in the direction of vn

    INTEGER :: je, jk, jb          !< index of edge, vert level, block
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart, i_rlend, i_nchdom
    INTEGER :: slev, elev          !< vertical start and end level
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


    ! number of child domains
    i_nchdom   = MAX(1,ptr_p%n_childdom)

    i_startblk = ptr_p%edges%start_blk(i_rlstart,1)
    i_endblk   = ptr_p%edges%end_blk(i_rlend,i_nchdom)

    ! pointer to the edge vertices
    ptr_ve => ptr_int%pos_on_tplane_e(:,:,7:8,:)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,dep_pts,arr_pts, &
!$OMP            pos_on_tplane_e,pn_cell,dn_cell,lvn_pos)
    DO jb = i_startblk, i_endblk

     CALL get_indices_e(ptr_p, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = slev, elev
        DO je = i_startidx, i_endidx



          ! Departure and arrival points. 
          ! Departure point 1 belongs to arrival point 1. The numbering of the 
          ! arrival points is adapted from the numbering of the edge vertices.
          !---------------------------------------------------------------------------
          !
          !  system orientation: -1 in this figure
          !
          ! D2\--------------\D1 
          !    \              \
          !     \      N       \       <- edge normal        [vn < 0]
          !      \      \       \
          !    A2 \------\-------\A1   <- triangle edge
          !        \              \
          !         \              \                         [vn > 0]
          !          \              \
          !        D2 \--------------\D1


          ! Determine the upwind cell
          ! Cell indices are chosen such that the direction from cell 1 to cell 2
          ! is the positive direction of the normal vector N.

          !
          ! Calculate backward trajectories, starting at the two edge vertices
          ! (arrival points). It is assumed that the velocity vector is constant 
          ! along the edge.
          !

          ! logical switch for MERGE operations: .TRUE. for p_vn >= 0
          lvn_pos = p_vn(je,jk,jb) >= 0._wp

          ! get line and block indices of upwind cell
          upwind_cell_idx(je,jk,jb) = MERGE(ptr_p%edges%cell_idx(je,jb,1),       &
            &                               ptr_p%edges%cell_idx(je,jb,2),lvn_pos)
          upwind_cell_blk(je,jk,jb) = MERGE(ptr_p%edges%cell_blk(je,jb,1),       &
            &                               ptr_p%edges%cell_blk(je,jb,2),lvn_pos)


          ! determine correct position on tangential plane (coordinates of 
          ! upwind cell center
          pos_on_tplane_e(1) = MERGE(ptr_int%pos_on_tplane_e(je,jb,1,1), &
            &                        ptr_int%pos_on_tplane_e(je,jb,2,1),lvn_pos)

          pos_on_tplane_e(2) = MERGE(ptr_int%pos_on_tplane_e(je,jb,1,2), &
            &                        ptr_int%pos_on_tplane_e(je,jb,2,2),lvn_pos)


          ! departure points of the departure cell. Point 1 belongs to edge-vertex 1, 
          ! point 2 belongs to edge_vertex 2.
          !
          ! position of departure point 1 in normal direction
          dep_pts(1,1) = ptr_ve(je,jb,1,1) - p_vn(je,jk,jb) * p_dt

          ! position of departure point 1 in tangential direction
          dep_pts(1,2) = ptr_ve(je,jb,1,2) - p_vt(je,jk,jb) * p_dt

          ! position of departure point 2 in normal direction
          dep_pts(2,1) = ptr_ve(je,jb,2,1) - p_vn(je,jk,jb) * p_dt

          ! position of departure point 2 in tangential direction
          dep_pts(2,2) = ptr_ve(je,jb,2,2) - p_vt(je,jk,jb) * p_dt



          ! Calculate position of departure region vertices in a translated
          ! coordinate system. The origin is located at the circumcenter
          ! of the upwind cell. The distance vectors point from the cell center
          ! to the vertices.
          !
          ! No sorting!
          !
          arr_pts(1,1:2) = ptr_ve(je,jb,1,1:2) - pos_on_tplane_e(1:2)

          arr_pts(2,1:2) = ptr_ve(je,jb,2,1:2) - pos_on_tplane_e(1:2)

          dep_pts(1,1:2) = dep_pts(1,1:2)      - pos_on_tplane_e(1:2)

          dep_pts(2,1:2) = dep_pts(2,1:2)      - pos_on_tplane_e(1:2)




          ! In a last step, these distance vectors are transformed into a rotated
          ! geographical coordinate system, which still has its origin at the circumcenter
          ! of the upwind cell. Now the coordinate axes point to local East and local
          ! North.
          !
          ! Determine primal and dual normals of the cell lying in the direction of vn
          pn_cell(1) = MERGE(ptr_p%edges%primal_normal_cell(je,jb,1)%v1,       &
            &                ptr_p%edges%primal_normal_cell(je,jb,2)%v1,lvn_pos)

          pn_cell(2) = MERGE(ptr_p%edges%primal_normal_cell(je,jb,1)%v2,       &
            &                ptr_p%edges%primal_normal_cell(je,jb,2)%v2,lvn_pos)

          dn_cell(1) = MERGE(ptr_p%edges%dual_normal_cell(je,jb,1)%v1,         &
            &                ptr_p%edges%dual_normal_cell(je,jb,2)%v1,lvn_pos)

          dn_cell(2) = MERGE(ptr_p%edges%dual_normal_cell(je,jb,1)%v2,         &
            &                ptr_p%edges%dual_normal_cell(je,jb,2)%v2,lvn_pos)


          ! components in longitudinal direction
          arrival_pts(je,1,1,jk,jb) = arr_pts(1,1) * pn_cell(1) &
            &                       + arr_pts(1,2) * dn_cell(1)

          arrival_pts(je,2,1,jk,jb) = arr_pts(2,1) * pn_cell(1) &
            &                       + arr_pts(2,2) * dn_cell(1)

          depart_pts(je,1,1,jk,jb)  = dep_pts(1,1) * pn_cell(1) &
            &                       + dep_pts(1,2) * dn_cell(1)

          depart_pts(je,2,1,jk,jb)  = dep_pts(2,1) * pn_cell(1) &
            &                       + dep_pts(2,2) * dn_cell(1)


          ! components in latitudinal direction
          arrival_pts(je,1,2,jk,jb) = arr_pts(1,1) * pn_cell(2) &
            &                       + arr_pts(1,2) * dn_cell(2)

          arrival_pts(je,2,2,jk,jb) = arr_pts(2,1) * pn_cell(2) &
            &                       + arr_pts(2,2) * dn_cell(2)

          depart_pts(je,1,2,jk,jb)  = dep_pts(1,1) * pn_cell(2) &
            &                       + dep_pts(1,2) * dn_cell(2)

          depart_pts(je,2,2,jk,jb)  = dep_pts(2,1) * pn_cell(2) &
            &                       + dep_pts(2,2) * dn_cell(2)


        ENDDO ! loop over edges
      ENDDO   ! loop over vertical levels
    END DO    ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE btraj_dreg_nosort


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
  SUBROUTINE btraj_o2( ptr_p, ptr_int, p_vn, p_vt, p_dthalf, p_cell_indices,  &
    &                  p_distv_bary, opt_rlstart, opt_rlend, opt_slev,        &
    &                  opt_elev )

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

    REAL(wp), INTENT(OUT) ::    &  !< distance vectors cell center --> barycenter of
      &  p_distv_bary(:,:,:,:)     !< departure region. (geographical coordinates)
                                   !< dim: (nproma,nlev,ptr_p%nblks_e,2)

    INTEGER, INTENT(OUT)  ::    &  !< line and block indices of cell centers in which the
      &  p_cell_indices(:,:,:,:)   !< computed barycenters are located.
                                   !< dim: (nproma,nlev,ptr_p%nblks_e,2)

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



          IF (p_vn(je,jk,jb) >= 0._wp) THEN

            !! we are in cell 1 !!

            ! line and block indices of neighboring cell with barycenter
            p_cell_indices(je,jk,jb,1) = ptr_p%edges%cell_idx(je,jb,1)
            p_cell_indices(je,jk,jb,2) = ptr_p%edges%cell_blk(je,jb,1)

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
            p_cell_indices(je,jk,jb,1) = ptr_p%edges%cell_idx(je,jb,2)
            p_cell_indices(je,jk,jb,2) = ptr_p%edges%cell_blk(je,jb,2)

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



  !-------------------------------------------------------------------------
  !>
  !! Divide flux area
  !!
  !! Flux area (aka. departure region) is subdivided according to its overlap 
  !! with the underlying grid.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2012-04-03)
  !!
  SUBROUTINE divide_flux_area(p_patch, p_int, p_vn, p_vt, depart_pts,        &
    &                         arrival_pts, dreg_patch0,                      &
    &                         dreg_patch1, dreg_patch2,                      &
    &                         patch1_cell_idx, patch1_cell_blk,              &
    &                         patch2_cell_idx, patch2_cell_blk,              &
    &                         opt_rlstart, opt_rlend, opt_slev,              &
    &                         opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::     &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) :: &  !< pointer to data structure for interpolation
      &  p_int

    REAL(wp), INTENT(IN)    ::  &  !< normal component of velocity vector at edge midpoints
      &  p_vn(:,:,:)

    REAL(wp), INTENT(IN)    ::  &  !< tangential component of velocity vector at
      &  p_vt(:,:,:)               !< edge midpoints

    REAL(wp), INTENT(IN) ::    &   !< coordinates of arrival points. The origin
      &  arrival_pts(:,:,:,:,:)    !< of the coordinate system is at the circumcenter of
                                   !< the upwind cell. Unit vectors point to local East
                                   !< and North. (geographical coordinates)
                                   !< dim: (nproma,2,2,nlev,ptr_p%nblks_e)

    REAL(wp), INTENT(IN) ::    &   !< coordinates of departure points. The origin
      &  depart_pts(:,:,:,:,:)     !< of the coordinate system is at the circumcenter of
                                   !< the upwind cell. Unit vectors point to local East
                                   !< and North. (geographical coordinates)
                                   !< dim: (nproma,2,2,nlev,ptr_p%nblks_e)

    REAL(wp), INTENT(OUT) ::    &  !< patch 0,1,2 of subdivided departure region
      & dreg_patch0(:,:,:,:,:), &  !< coordinates
      & dreg_patch1(:,:,:,:,:), &  !< dim: (nproma,4,2,nlev,ptr_p%nblks_e)
      & dreg_patch2(:,:,:,:,:)

    INTEGER, INTENT(OUT)  ::    &  !< line and block indices of underlying cell
      & patch1_cell_idx(:,:,:), &  !< dim: (nproma,nlev,ptr_p%nblks_e)
      & patch1_cell_blk(:,:,:), &
      & patch2_cell_idx(:,:,:), &
      & patch2_cell_blk(:,:,:)


    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    INTEGER :: je, jk, jb              !< loop index of edge, vert level, block
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart, i_rlend, i_nchdom
    INTEGER :: slev, elev              !< vertical start and end level

    TYPE(t_line) ::                  & !< departure-line segment
      &  fl_line(nproma)

    TYPE(t_line) ::                  & !< departure area edges
      &  fl_e1,                      & !< edge 1
      &  fl_e2                         !< edge 2

    TYPE(t_line) ::                  & !< triangle edge
      &  tri_line1(nproma),          &
      &  tri_line2(nproma)

    TYPE(t_geographical_coordinates), POINTER :: & !< pointer to coordinates of vertex3
      &  ptr_v3(:,:,:)
    TYPE(t_geographical_coordinates), POINTER :: & !< pointer to coordinates of 
      &  ptr_bfcc(:,:,:,:)                         !< of butterfly cell centers


    REAL(wp) :: ps1(nproma,2),          & !< coordinates of intersection 
      &         ps2(nproma,2)             !< points S1, S2

    REAL(wp) :: pi1(nproma,2),          & !< coordinates of intersection 
      &         pi2(nproma,2)             !< points I1, I2

    REAL(wp) :: bf_cc(2,2)                !< coordinates of butterfly cell centers

    LOGICAL :: lintersect_line1, lintersect_line2
    LOGICAL :: lintersect_e2_line1, lintersect_e1_line2
    LOGICAL :: lvn_pos, lvn_sys_pos
    INTEGER :: cfl_stable(nproma)


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_traj: divide_flux_area'
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
      elev = p_patch%nlev
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


    ! number of child domains
    i_nchdom   = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)


    ! pointer to coordinates of vertex3 (i.e. vertex which belongs to 
    ! the upwind cell but does not belong to the current edge.)
    !
    ptr_v3   => p_int%pos_on_tplane_c_edge(:,:,:,3)

    ! pointer to coordinates of butterfly cell centers
    !
    ptr_bfcc => p_int%pos_on_tplane_c_edge(:,:,:,4:5)



!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,lvn_pos,fl_line,tri_line1,  &
!$OMP            tri_line2,cfl_stable,lvn_sys_pos,lintersect_line1,       &
!$OMP            lintersect_line2,fl_e1,fl_e2,lintersect_e2_line1,        &
!$OMP            lintersect_e1_line2,ps1,ps2,pi1,pi2,bf_cc)
    DO jb = i_startblk, i_endblk

     CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = slev, elev

         DO je = i_startidx, i_endidx

          lvn_pos = p_vn(je,jk,jb) >= 0._wp

          !
          ! I) check whether departure-line segment is located 
          !    within the upwind row
          !
          ! get flux area departure-line segment
          !
          fl_line(je)%p1%lon = depart_pts(je,1,1,jk,jb)
          fl_line(je)%p1%lat = depart_pts(je,1,2,jk,jb)
          fl_line(je)%p2%lon = depart_pts(je,2,1,jk,jb)
          fl_line(je)%p2%lat = depart_pts(je,2,2,jk,jb)

          ! get triangle edge 1 (A1V3)
          !
          tri_line1(je)%p1%lon = arrival_pts(je,1,1,jk,jb)
          tri_line1(je)%p1%lat = arrival_pts(je,1,2,jk,jb)
          tri_line1(je)%p2%lon = MERGE(ptr_v3(je,jb,1)%lon,ptr_v3(je,jb,2)%lon,lvn_pos)
          tri_line1(je)%p2%lat = MERGE(ptr_v3(je,jb,1)%lat,ptr_v3(je,jb,2)%lat,lvn_pos)

          ! get triangle edge 2 (A2V3)
          !
          tri_line2(je)%p1%lon = arrival_pts(je,2,1,jk,jb)
          tri_line2(je)%p1%lat = arrival_pts(je,2,2,jk,jb)
          tri_line2(je)%p2%lon = MERGE(ptr_v3(je,jb,1)%lon,ptr_v3(je,jb,2)%lon,lvn_pos)
          tri_line2(je)%p2%lat = MERGE(ptr_v3(je,jb,1)%lat,ptr_v3(je,jb,2)%lat,lvn_pos)

 
          ! fl_line%p1 -> fl_line%p2 -> vert3 and
          ! fl_line%p1 -> fl_line%p2 -> arr1 must have different orientations
          ! ccw() * ccw() = -1
          cfl_stable(je) = (ccw(fl_line(je)%p1,fl_line(je)%p2,tri_line1(je)%p2)  &
            &            *  ccw(fl_line(je)%p1,fl_line(je)%p2,tri_line1(je)%p1))

        ENDDO ! loop over edges

        !
        ! I) check whether departure-line segment is still located within the upwind row
        !
        IF ( ANY(cfl_stable(i_startidx:i_endidx) == 1) ) THEN
          WRITE(message_text,'(a,2i4)') 'horizontal CFL number exceeded at jk,jb = ',&
            &  jk,jb
          CALL message(TRIM(routine), TRIM(message_text))
        ENDIF




        DO je = i_startidx, i_endidx

          lvn_pos = p_vn(je,jk,jb) >= 0._wp

          lvn_sys_pos = (p_vn(je,jk,jb) * p_patch%edges%system_orientation(je,jb)) >= 0._wp


          ! II) Check whether departure line intersects with edge A1V3 and/or A2V3


          ! does departure-line segment intersect with A1V3?
          !
!CDIR NEXPAND(lintersect)
          lintersect_line1 = lintersect(fl_line(je), tri_line1(je))

          ! does departure-line segment intersect with A2V3?
          !
!CDIR NEXPAND(lintersect)
          lintersect_line2 = lintersect(fl_line(je), tri_line2(je))



          IF ( (p_patch%edges%system_orientation(je,jb) * p_vt(je,jk,jb)) >= 0._wp ) THEN


            ! get flux area edge 2
            !
            fl_e2%p1%lon = arrival_pts(je,2,1,jk,jb)
            fl_e2%p1%lat = arrival_pts(je,2,2,jk,jb)
            fl_e2%p2%lon = depart_pts (je,2,1,jk,jb)
            fl_e2%p2%lat = depart_pts (je,2,2,jk,jb)


            ! get intersection point of flux area edge 2 with triangle edge 1
            !
!CDIR NEXPAND(lintersect)
            lintersect_e2_line1 = lintersect(fl_e2, tri_line1(je))


            IF ( lintersect_line1 .AND. lintersect_line2 ) THEN
!!$WRITE(0,*) "CASE I: je,jk,jb", je,jk,jb
              ! CASE I
              ! Compute intersection point of fl_line with tri_line1
              ! Compute intersection point of fl_line with tri_line2
              !
              ps1(je,1:2) = line_intersect(fl_line(je), tri_line1(je))
              ps2(je,1:2) = line_intersect(fl_line(je), tri_line2(je))

              ! store corners of flux area patches (counterclockwise)
              ! patch 0
              ! vn > 0: A1 A2 S2 S1
              ! vn < 0: A1 S1 S2 A2
              !
              dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
              dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(je,2,1,jk,jb), &
                &                           ps1(je,1),lvn_sys_pos)
              dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(je,2,2,jk,jb), &
                &                           ps1(je,2),lvn_sys_pos)
              dreg_patch0(je,3,1:2,jk,jb) = ps2(je,1:2)
              dreg_patch0(je,4,1,jk,jb)   = MERGE(ps1(je,1),arrival_pts(je,2,1,jk,jb), &
                &                           lvn_sys_pos)
              dreg_patch0(je,4,2,jk,jb)   = MERGE(ps1(je,2),arrival_pts(je,2,2,jk,jb), &
                &                           lvn_sys_pos)

              ! patch 1
              ! vn > 0: A1 S1 D1 A1 (degenerated)
              ! vn < 0: A1 D1 S1 A1 (degenerated)
              !
              dreg_patch1(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
              dreg_patch1(je,2,1,jk,jb)   = MERGE(ps1(je,1),                    &
                &                           depart_pts(je,1,1,jk,jb), lvn_sys_pos)
              dreg_patch1(je,2,2,jk,jb)   = MERGE(ps1(je,2),                    &
                &                           depart_pts(je,1,2,jk,jb), lvn_sys_pos)
              dreg_patch1(je,3,1,jk,jb)   = MERGE(depart_pts(je,1,1,jk,jb),  &
                &                           ps1(je,1),lvn_sys_pos)
              dreg_patch1(je,3,2,jk,jb)   = MERGE(depart_pts(je,1,2,jk,jb),  &
                &                           ps1(je,2), lvn_sys_pos)
              dreg_patch1(je,4,1:2,jk,jb) = dreg_patch1(je,1,1:2,jk,jb)


              ! patch 2
              ! vn > 0: A2 D2 S2 A2 (degenerated)
              ! vn < 0: A2 S2 D2 A2 (degenerated)
              !
              dreg_patch2(je,1,1:2,jk,jb) = arrival_pts(je,2,1:2,jk,jb)
              dreg_patch2(je,2,1,jk,jb)   = MERGE(depart_pts(je,2,1,jk,jb),  &
                &                           ps2(je,1), lvn_sys_pos)
              dreg_patch2(je,2,2,jk,jb)   = MERGE(depart_pts(je,2,2,jk,jb),  &
                &                           ps2(je,2), lvn_sys_pos)
              dreg_patch2(je,3,1,jk,jb)   = MERGE(ps2(je,1),                    &
                &                           depart_pts(je,2,1,jk,jb), lvn_sys_pos)
              dreg_patch2(je,3,2,jk,jb)   = MERGE(ps2(je,2),                    &
                &                           depart_pts(je,2,2,jk,jb), lvn_sys_pos)
              dreg_patch2(je,4,1:2,jk,jb) = dreg_patch2(je,1,1:2,jk,jb)



            ELSE IF ( lintersect_line1 .AND. (.NOT. lintersect_line2) ) THEN
!!$WRITE(0,*) "CASE IIa: je,jk,jb", je,jk,jb
              ! CASE 2a
              ! Compute intersection point of fl_line with tri_line1
              !
              ps1(je,1:2) = line_intersect(fl_line(je), tri_line1(je))

              ! store corners of flux area patches (counterclockwise)
              ! patch 0
              ! vn > 0: A1 A2 D2 S1
              ! vn < 0: A1 S1 D2 A2
              !
              dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
              dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(je,2,1,jk,jb), &
                &                           ps1(je,1),lvn_sys_pos)
              dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(je,2,2,jk,jb), &
                &                           ps1(je,2),lvn_sys_pos)
              dreg_patch0(je,3,1:2,jk,jb) = depart_pts(je,2,1:2,jk,jb)
              dreg_patch0(je,4,1,jk,jb)   = MERGE(ps1(je,1),arrival_pts(je,2,1,jk,jb), &
                &                           lvn_sys_pos)
              dreg_patch0(je,4,2,jk,jb)   = MERGE(ps1(je,2),arrival_pts(je,2,2,jk,jb), &
                &                           lvn_sys_pos)

              ! patch 1
              ! vn > 0: A1 S1 D1 A1 (degenerated)
              ! vn < 0: A1 D1 S1 A1 (degenerated)
              !
              dreg_patch1(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
              dreg_patch1(je,2,1,jk,jb)   = MERGE(ps1(je,1),                    &
                &                           depart_pts(je,1,1,jk,jb), lvn_sys_pos)
              dreg_patch1(je,2,2,jk,jb)   = MERGE(ps1(je,2),                    &
                &                           depart_pts(je,1,2,jk,jb), lvn_sys_pos)
              dreg_patch1(je,3,1,jk,jb)   = MERGE(depart_pts(je,1,1,jk,jb),  &
                &                           ps1(je,1),lvn_sys_pos)
              dreg_patch1(je,3,2,jk,jb)   = MERGE(depart_pts(je,1,2,jk,jb),  &
                &                           ps1(je,2), lvn_sys_pos)
              dreg_patch1(je,4,1:2,jk,jb) = dreg_patch1(je,1,1:2,jk,jb)


              ! patch 2 (non-existing)
              !
              dreg_patch2(je,1,1:2,jk,jb) = 0._wp
              dreg_patch2(je,2,1:2,jk,jb) = 0._wp
              dreg_patch2(je,3,1:2,jk,jb) = 0._wp
              dreg_patch2(je,4,1:2,jk,jb) = 0._wp


            ELSE IF ( lintersect_e2_line1 ) THEN
!!$WRITE(0,*) "CASE IIIa: je,jk,jb", je,jk,jb
              ! CASE 3a
              ! Compute intersection point of fl_e2 with tri_line1
              !
              pi1(je,1:2) = line_intersect(fl_e2, tri_line1(je))

              ! store corners of flux area patches (counterclockwise)
              ! patch 0
              ! vn > 0: A1 A2 I1 A1 (degenerated)
              ! vn < 0: A1 I1 A2 A1 (degenerated)
              !
              dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
              dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(je,2,1,jk,jb), &
                &                           pi1(je,1),lvn_sys_pos)
              dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(je,2,2,jk,jb), &
                &                           pi1(je,2),lvn_sys_pos)
              dreg_patch0(je,3,1,jk,jb)   = MERGE(pi1(je,1),                    &
                &                           arrival_pts(je,2,1,jk,jb),lvn_sys_pos)
              dreg_patch0(je,3,2,jk,jb)   = MERGE(pi1(je,2),                    &
                &                           arrival_pts(je,2,2,jk,jb),lvn_sys_pos) 
              dreg_patch0(je,4,1:2,jk,jb) = dreg_patch0(je,1,1:2,jk,jb)


              ! patch 1
              ! vn > 0: A1 I1 D2 D1
              ! vn < 0: A1 D1 D2 I1
              !
              dreg_patch1(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
              dreg_patch1(je,2,1,jk,jb)   = MERGE(pi1(je,1),                    &
                &                           depart_pts(je,1,1,jk,jb), lvn_sys_pos)
              dreg_patch1(je,2,2,jk,jb)   = MERGE(pi1(je,2),                    &
                &                           depart_pts(je,1,2,jk,jb), lvn_sys_pos)
              dreg_patch1(je,3,1:2,jk,jb) = depart_pts(je,2,1:2,jk,jb)
              dreg_patch1(je,4,1,jk,jb)   = MERGE(depart_pts(je,1,1,jk,jb),  &
                &                           pi1(je,1),lvn_sys_pos)
              dreg_patch1(je,4,2,jk,jb)   = MERGE(depart_pts(je,1,2,jk,jb),  &
                &                           pi1(je,2), lvn_sys_pos)


              ! patch 2 (non-existing)
              !
              dreg_patch2(je,1,1:2,jk,jb) = 0._wp
              dreg_patch2(je,2,1:2,jk,jb) = 0._wp
              dreg_patch2(je,3,1:2,jk,jb) = 0._wp
              dreg_patch2(je,4,1:2,jk,jb) = 0._wp

             ENDIF



          ELSE IF ( (p_patch%edges%system_orientation(je,jb) * p_vt(je,jk,jb)) < 0._wp ) THEN


            ! get flux area edge 1
            !
            fl_e1%p1%lon = arrival_pts(je,1,1,jk,jb)
            fl_e1%p1%lat = arrival_pts(je,1,2,jk,jb)
            fl_e1%p2%lon = depart_pts (je,1,1,jk,jb)
            fl_e1%p2%lat = depart_pts (je,1,2,jk,jb)

            ! get intersection point of flux area edge 1 with triangle edge 2
            !
!CDIR NEXPAND(lintersect)
            lintersect_e1_line2 = lintersect(fl_e1, tri_line2(je))


            IF ( lintersect_line1 .AND. lintersect_line2 ) THEN
!!$WRITE(0,*) "CASE Iagain: je,jk,jb", je,jk,jb
              ! CASE I
              ! Compute intersection point of fl_line with tri_line1
              ! Compute intersection point of fl_line with tri_line2
              !
              ps1(je,1:2) = line_intersect(fl_line(je), tri_line1(je))
              ps2(je,1:2) = line_intersect(fl_line(je), tri_line2(je))

              ! store corners of flux area patches (counterclockwise)
              ! patch 0
              ! vn > 0: A1 A2 S2 S1
              ! vn < 0: A1 S1 S2 A2
              !
              dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
              dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(je,2,1,jk,jb), &
                &                           ps1(je,1),lvn_sys_pos)
              dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(je,2,2,jk,jb), &
                &                           ps1(je,2),lvn_sys_pos)
              dreg_patch0(je,3,1:2,jk,jb) = ps2(je,1:2)
              dreg_patch0(je,4,1,jk,jb)   = MERGE(ps1(je,1),arrival_pts(je,2,1,jk,jb), &
                &                           lvn_sys_pos)
              dreg_patch0(je,4,2,jk,jb)   = MERGE(ps1(je,2),arrival_pts(je,2,2,jk,jb), &
                &                           lvn_sys_pos)

              ! patch 1
              ! vn > 0: A1 S1 D1 A1 (degenerated)
              ! vn < 0: A1 D1 S1 A1 (degenerated)
              !
              dreg_patch1(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
              dreg_patch1(je,2,1,jk,jb)   = MERGE(ps1(je,1),                    &
                &                           depart_pts(je,1,1,jk,jb), lvn_sys_pos)
              dreg_patch1(je,2,2,jk,jb)   = MERGE(ps1(je,2),                    &
                &                           depart_pts(je,1,2,jk,jb), lvn_sys_pos)
              dreg_patch1(je,3,1,jk,jb)   = MERGE(depart_pts(je,1,1,jk,jb),  &
                &                           ps1(je,1),lvn_sys_pos)
              dreg_patch1(je,3,2,jk,jb)   = MERGE(depart_pts(je,1,2,jk,jb),  &
                &                           ps1(je,2), lvn_sys_pos)
              dreg_patch1(je,4,1:2,jk,jb) = dreg_patch1(je,1,1:2,jk,jb)

              ! patch 2
              ! vn > 0: A2 D2 S2 A2 (degenerated)
              ! vn < 0: A2 S2 D2 A2 (degenerated)
              !
              dreg_patch2(je,1,1:2,jk,jb) = arrival_pts(je,2,1:2,jk,jb)
              dreg_patch2(je,2,1,jk,jb)   = MERGE(depart_pts(je,2,1,jk,jb),  &
                &                           ps2(je,1), lvn_sys_pos)
              dreg_patch2(je,2,2,jk,jb)   = MERGE(depart_pts(je,2,2,jk,jb),  &
                &                           ps2(je,2), lvn_sys_pos)
              dreg_patch2(je,3,1,jk,jb)   = MERGE(ps2(je,1),                    &
                &                           depart_pts(je,2,1,jk,jb), lvn_sys_pos)
              dreg_patch2(je,3,2,jk,jb)   = MERGE(ps2(je,2),                    &
                &                           depart_pts(je,2,2,jk,jb), lvn_sys_pos)
              dreg_patch2(je,4,1:2,jk,jb) = dreg_patch2(je,1,1:2,jk,jb)




            ELSE IF ( lintersect_line2 .AND. (.NOT. lintersect_line1) ) THEN
!!$WRITE(0,*) "CASE IIb: je,jk,jb", je,jk,jb
              ! CASE 2b
              ! Compute intersection point of fl_line with tri_line2
              !

              ps2(je,1:2) = line_intersect(fl_line(je), tri_line2(je))

              ! store corners of flux area patches (counterclockwise)
              ! patch 0
              ! vn > 0: A1 A2 S2 D1
              ! vn < 0: A1 D1 S2 A2
              !
              dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
              dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(je,2,1,jk,jb), &
                &                           depart_pts(je,1,1,jk,jb),lvn_sys_pos)
              dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(je,2,2,jk,jb), &
                &                           depart_pts(je,1,2,jk,jb),lvn_sys_pos)
              dreg_patch0(je,3,1:2,jk,jb) = ps2(je,1:2)
              dreg_patch0(je,4,1,jk,jb)   = MERGE(depart_pts(je,1,1,jk,jb),  &
                &                           arrival_pts(je,2,1,jk,jb), lvn_sys_pos)
              dreg_patch0(je,4,2,jk,jb)   = MERGE(depart_pts(je,1,2,jk,jb),  &
                &                           arrival_pts(je,2,2,jk,jb), lvn_sys_pos)


              ! patch 1 (non-existing)
              !
              dreg_patch1(je,1,1:2,jk,jb) = 0._wp
              dreg_patch1(je,2,1:2,jk,jb) = 0._wp
              dreg_patch1(je,3,1:2,jk,jb) = 0._wp
              dreg_patch1(je,4,1:2,jk,jb) = 0._wp


              ! patch 2
              ! vn > 0: A2 D2 S2 A2 (degenerated)
              ! vn < 0: A2 S2 D2 A2 (degenerated)
              !
              dreg_patch2(je,1,1:2,jk,jb) = arrival_pts(je,2,1:2,jk,jb)
              dreg_patch2(je,2,1,jk,jb)   = MERGE(depart_pts(je,2,1,jk,jb),  &
                &                           ps2(je,1), lvn_sys_pos)
              dreg_patch2(je,2,2,jk,jb)   = MERGE(depart_pts(je,2,2,jk,jb),  &
                &                           ps2(je,2), lvn_sys_pos)
              dreg_patch2(je,3,1,jk,jb)   = MERGE(ps2(je,1),                    &
                &                           depart_pts(je,2,1,jk,jb), lvn_sys_pos)
              dreg_patch2(je,3,2,jk,jb)   = MERGE(ps2(je,2),                    &
                &                           depart_pts(je,2,2,jk,jb), lvn_sys_pos)
              dreg_patch2(je,4,1:2,jk,jb) = dreg_patch2(je,1,1:2,jk,jb)



            ELSE IF ( lintersect_e1_line2 ) THEN
!!$WRITE(0,*) "CASE IIIb: je,jk,jb", je,jk,jb
              ! CASE 3b
              ! Compute intersection point of fl_e1 with tri_line2
              pi2(je,1:2) = line_intersect(fl_e1, tri_line2(je))

              ! store corners of flux area patches (counterclockwise)
              ! patch 0
              ! vn > 0: A1 A2 I2 A1 (degenerated)
              ! vn < 0: A1 I2 A2 A1 (degenerated)
              !
              dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
              dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(je,2,1,jk,jb), &
                &                           pi2(je,1),lvn_sys_pos)
              dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(je,2,2,jk,jb), &
                &                           pi2(je,2),lvn_sys_pos)
              dreg_patch0(je,3,1,jk,jb)   = MERGE(pi2(je,1),                    &
                &                           arrival_pts(je,2,1,jk,jb),lvn_sys_pos)
              dreg_patch0(je,3,2,jk,jb)   = MERGE(pi2(je,2),                    &
                &                           arrival_pts(je,2,2,jk,jb),lvn_sys_pos) 
              dreg_patch0(je,4,1:2,jk,jb) = dreg_patch0(je,1,1:2,jk,jb)


              ! patch 1 (non-existing)
              !
              dreg_patch1(je,1,1:2,jk,jb) = 0._wp
              dreg_patch1(je,2,1:2,jk,jb) = 0._wp
              dreg_patch1(je,3,1:2,jk,jb) = 0._wp
              dreg_patch1(je,4,1:2,jk,jb) = 0._wp



              ! patch 2
              ! vn > 0: A2 D2 D1 I2
              ! vn < 0: A2 I2 D1 D2
              !
              dreg_patch2(je,1,1:2,jk,jb) = arrival_pts(je,2,1:2,jk,jb)
              dreg_patch2(je,2,1,jk,jb)   = MERGE(depart_pts(je,2,1,jk,jb),  &
                &                           pi2(je,1), lvn_sys_pos)
              dreg_patch2(je,2,2,jk,jb)   = MERGE(depart_pts(je,2,2,jk,jb),  &
                &                           pi2(je,2), lvn_sys_pos)
              dreg_patch2(je,3,1:2,jk,jb) = depart_pts(je,1,1:2,jk,jb)
              dreg_patch2(je,4,1,jk,jb)   = MERGE(pi2(je,1),                    &
                &                           depart_pts(je,2,1,jk,jb), lvn_sys_pos)
              dreg_patch2(je,4,2,jk,jb)   = MERGE(pi2(je,2),                    &
                &                           depart_pts(je,2,2,jk,jb), lvn_sys_pos)

            ENDIF  ! lintersect
          ENDIF  ! p_vt



!!$WRITE(0,*)
!!$WRITE(0,*) "x: dreg_patch0(je,1:4,1,jk,jb): ", dreg_patch0(je,1:4,1,jk,jb)
!!$WRITE(0,*) "y: dreg_patch0(je,1:4,2,jk,jb): ", dreg_patch0(je,1:4,2,jk,jb)
!!$WRITE(0,*)
!!$WRITE(0,*) "x: dreg_patch1(je,1:4,1,jk,jb): ", dreg_patch1(je,1:4,1,jk,jb)
!!$WRITE(0,*) "y: dreg_patch1(je,1:4,2,jk,jb): ", dreg_patch1(je,1:4,2,jk,jb)
!!$WRITE(0,*) "cell center of underlying cell (patch1) (x,y): ", bf_cc(1,1:2)
!!$WRITE(0,*)
!!$WRITE(0,*) "x: dreg_patch2(je,1:4,1,jk,jb): ", dreg_patch2(je,1:4,1,jk,jb)
!!$WRITE(0,*) "y: dreg_patch2(je,1:4,2,jk,jb): ", dreg_patch2(je,1:4,2,jk,jb)
!!$WRITE(0,*) "cell center of underlying cell (patch2) (x,y): ", bf_cc(2,1:2)
!!$WRITE(0,*)
!!$WRITE(0,*) "up_vert1 (x,y): ",MERGE(p_int%pos_on_tplane_c_edge(je,jb,1,1)%lon,          &
!!$  &                                 p_int%pos_on_tplane_c_edge(je,jb,2,1)%lon,lvn_pos), &
!!$  &                           MERGE(p_int%pos_on_tplane_c_edge(je,jb,1,1)%lat,          &
!!$  &                                 p_int%pos_on_tplane_c_edge(je,jb,2,1)%lat,lvn_pos)
!!$WRITE(0,*) "up_vert2 (x,y): ",MERGE(p_int%pos_on_tplane_c_edge(je,jb,1,2)%lon,          &
!!$  &                                 p_int%pos_on_tplane_c_edge(je,jb,2,2)%lon,lvn_pos), &
!!$  &                           MERGE(p_int%pos_on_tplane_c_edge(je,jb,1,2)%lat,          &
!!$  &                                 p_int%pos_on_tplane_c_edge(je,jb,2,2)%lat,lvn_pos)
!!$WRITE(0,*) "up_vert3 (x,y): ",MERGE(p_int%pos_on_tplane_c_edge(je,jb,1,3)%lon,          &
!!$  &                                 p_int%pos_on_tplane_c_edge(je,jb,2,3)%lon,lvn_pos), &
!!$  &                           MERGE(p_int%pos_on_tplane_c_edge(je,jb,1,3)%lat,          &
!!$  &                                 p_int%pos_on_tplane_c_edge(je,jb,2,3)%lat,lvn_pos)
!!$WRITE(0,*)
!!$WRITE(0,*) "system orientation, v_t: ", p_patch%edges%system_orientation(je,jb), p_vt(je,jk,jb)
!!$WRITE(0,*)


          ! get coordinates of butterfly cell centers (from upwind cell)
          bf_cc(1,1) = MERGE(ptr_bfcc(je,jb,1,1)%lon,         &
            &                ptr_bfcc(je,jb,2,1)%lon, lvn_pos ) 
          bf_cc(1,2) = MERGE(ptr_bfcc(je,jb,1,1)%lat,         &
            &                ptr_bfcc(je,jb,2,1)%lat, lvn_pos )
          bf_cc(2,1) = MERGE(ptr_bfcc(je,jb,1,2)%lon,         &
            &                ptr_bfcc(je,jb,2,2)%lon, lvn_pos ) 
          bf_cc(2,2) = MERGE(ptr_bfcc(je,jb,1,2)%lat,         &
            &                ptr_bfcc(je,jb,2,2)%lat, lvn_pos ) 

          ! patch 1 in translated system
          !
          dreg_patch1(je,1,1:2,jk,jb) = dreg_patch1(je,1,1:2,jk,jb) - bf_cc(1,1:2)
          dreg_patch1(je,2,1:2,jk,jb) = dreg_patch1(je,2,1:2,jk,jb) - bf_cc(1,1:2)
          dreg_patch1(je,3,1:2,jk,jb) = dreg_patch1(je,3,1:2,jk,jb) - bf_cc(1,1:2)
          dreg_patch1(je,4,1:2,jk,jb) = dreg_patch1(je,4,1:2,jk,jb) - bf_cc(1,1:2)


          ! patch 2 in translated system
          !
          dreg_patch2(je,1,1:2,jk,jb) = dreg_patch2(je,1,1:2,jk,jb) - bf_cc(2,1:2)
          dreg_patch2(je,2,1:2,jk,jb) = dreg_patch2(je,2,1:2,jk,jb) - bf_cc(2,1:2)
          dreg_patch2(je,3,1:2,jk,jb) = dreg_patch2(je,3,1:2,jk,jb) - bf_cc(2,1:2)
          dreg_patch2(je,4,1:2,jk,jb) = dreg_patch2(je,4,1:2,jk,jb) - bf_cc(2,1:2)



          ! store global index of the underlying grid cell
          !
          patch1_cell_idx(je,jk,jb) = MERGE(p_patch%edges%butterfly_idx(je,jb,1,1), &
            &                               p_patch%edges%butterfly_idx(je,jb,2,1), &
            &                               lvn_pos) 
          patch2_cell_idx(je,jk,jb) = MERGE(p_patch%edges%butterfly_idx(je,jb,1,2), &
            &                               p_patch%edges%butterfly_idx(je,jb,2,2), &
            &                               lvn_pos)

          patch1_cell_blk(je,jk,jb) = MERGE(p_patch%edges%butterfly_blk(je,jb,1,1), &
            &                               p_patch%edges%butterfly_blk(je,jb,2,1), &
            &                               lvn_pos)
          patch2_cell_blk(je,jk,jb) = MERGE(p_patch%edges%butterfly_blk(je,jb,1,2), &
            &                               p_patch%edges%butterfly_blk(je,jb,2,2), &
            &                               lvn_pos)


        ENDDO ! loop over edges


      ENDDO   ! loop over vertical levels
    END DO    ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE divide_flux_area

END MODULE mo_advection_traj

