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
MODULE mo_advection_utils

  USE mo_advection_config,    ONLY: shape_func, zeta, eta, wgt_zeta, wgt_eta
  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_interpolation,       ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants,      ONLY: min_rlcell_int, min_rledge_int
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_timer,               ONLY: timer_start, timer_stop, timers_level, new_timer

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: laxfr_upflux, laxfr_upflux_v, ptr_delp_mc_now, ptr_delp_mc_new,  &
    &       back_traj_o1, back_traj_dreg_o1, back_traj_o2,                   &
    &       prep_gauss_quadrature_q, prep_gauss_quadrature_cpoor,            &
    &       prep_gauss_quadrature_c
  
  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  ! In order to avoid circular dependencies these two pointers
  ! have been moved from mo_advection_stepping to this module.
  REAL(wp), POINTER ::  &
    &  ptr_delp_mc_now(:,:,:) => NULL() !< pointer to old layer thickness
                                        !< at cell center
  REAL(wp), POINTER ::  &
    &  ptr_delp_mc_new(:,:,:) => NULL() !< pointer to new layer thickness
                                        !< at cell center

  INTEGER :: timer_back_traj_o1   = 0

CONTAINS

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Integrates tracer continuity equation from old time step to new time step
  !!
  !! This subroutine integrates the tracer continuity equation using a simple
  !! forward time step.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-02-24)
  !!
  SUBROUTINE tupdate_tracer( p_patch, p_dtime, p_tracer_now, p_density_c_now, &
    &                        p_density_c_new, p_fluxdiv_c, p_tracer_new,      &
    &                        opt_rlstart, opt_rlend )

    TYPE(t_patch), TARGET, INTENT(IN) ::   & !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) :: p_dtime      !< time step

    REAL(wp), INTENT(IN) ::     &        !< tracer field at current time
      &  p_tracer_now(:,:,:)             !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::     &        !< density (or layer thickness) at current time
      &  p_density_c_now(:,:,:)          !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::     &        !< density (or layer thickness) at new time
      &  p_density_c_new(:,:,:)          !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::     &        !< flux divergence at current time
      &  p_fluxdiv_c(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(INOUT) ::  &        !< tracer field at current time
      &  p_tracer_new(:,:,:)             !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL :: &   !< optional: refinement control start level
     &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: &   !< optional: refinement control end level
     &  opt_rlend                        !< (to avoid calculation of halo points)

    INTEGER :: nlev                      !< number of full levels
    INTEGER :: jb, jk, jc                !< loop indices
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_rlstart, i_rlend, i_nchdom !< start and end values of refined grid
    INTEGER :: i_startidx, i_endidx
   !-----------------------------------------------------------------------

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell_int
    ENDIF

    ! number of vertical levels
    nlev = p_patch%nlev

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
          DO jb = i_startblk, i_endblk

           CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

           DO jk = 1, nlev

             DO jc = i_startidx, i_endidx

               p_tracer_new(jc,jk,jb) =                                     &
                 &   ( p_tracer_now(jc,jk,jb) * p_density_c_now(jc,jk,jb)   &
                 &    - p_dtime * p_fluxdiv_c(jc,jk,jb) )                   &
                 &    / p_density_c_new(jc,jk,jb)

             ENDDO
           ENDDO
         ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE tupdate_tracer



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
  SUBROUTINE back_traj_o1( ptr_p, ptr_int, p_vn, p_vt, p_dthalf, p_cell_indices, &
    &                     p_distv_bary, opt_rlstart, opt_rlend, opt_slev,        &
    &                     opt_elev )

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
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,z_ntdistv_bary,lvn_pos)
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
!$OMP END DO
!$OMP END PARALLEL
    
    IF (timers_level > 5) CALL timer_stop(timer_back_traj_o1)

  END SUBROUTINE back_traj_o1



  !-------------------------------------------------------------------------
  !>
  !! Computation of first order backward trajectories for FFSL transport scheme
  !!
  !! Computes backward trajectories in order to determine an approximation to the
  !! departure region. Here, the departure region is approximated by a rhomboidal
  !! region, using a simple first order accurate trajectory computation. Computations
  !! are performed on a plane tangent to the edge midpoint. Coordinate axes point into
  !! the local normal and tangential direction.
  !! Once the 4 vertices of the departure region are known, the distance vector
  !! between the circumcenter of the upstream cell and the vertices is computed.
  !! In a final step, these vectors are transformed into a different coordinate system
  !! which has its origin at the circumcenter. The coordinate axes point to the local
  !! east and local north. This subroutine may be combined with any reconstruction method
  !! for the subgrid distribution.
  !! Note: Care has to be taken that the vertices of the departure region are stored in
  !! counterclockwise order. This is a requirement of the gaussian quadrature which
  !! follows lateron.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-05-12)
  !!
  SUBROUTINE back_traj_dreg_o1( ptr_p, ptr_int, p_vn, p_vt, p_dt, p_cell_indices, &
    &                     p_coords_dreg_v, opt_rlstart, opt_rlend, opt_slev,      &
    &                     opt_elev )

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
                                   !< dim: (nproma,nlev,ptr_p%nblks_e,4,2)

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

    REAL(wp) ::            &       !< coordinates of departure region vertices
      &  pos_dreg_vert_e(3:4,2)    !< in edge-based coordinate system

    REAL(wp) ::            &       !< coordinates of departure region vertices
      &  pos_dreg_vert_c(4,2)      !< as seen from translated coordinate system.
                                   !< origin at circumcenter of upwind cell

    REAL(wp) ::            &       !< position on tangential plane depending
      &  pos_on_tplane_e(2)        !< on the sign of vn

    REAL(wp) ::            &       !< primal and dual normals of cell lying
      &  pn_cell(2), dn_cell(2)    !< in the direction of vn

    INTEGER :: je, jk, jb          !< index of edge, vert level, block
    INTEGER :: nlev                !< number of full levels
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

    ! number of vertical levels
    nlev = ptr_p%nlev

    ! number of child domains
    i_nchdom   = MAX(1,ptr_p%n_childdom)

    i_startblk = ptr_p%edges%start_blk(i_rlstart,1)
    i_endblk   = ptr_p%edges%end_blk(i_rlend,i_nchdom)

    ! pointer to the edge vertices
    ptr_ve => ptr_int%pos_on_tplane_e(:,:,7:8,:)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,pos_dreg_vert_e,          &
!$OMP            pos_dreg_vert_c,pos_on_tplane_e,pn_cell,dn_cell,lvn_pos)
    DO jb = i_startblk, i_endblk

     CALL get_indices_e(ptr_p, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = slev, elev
        DO je = i_startidx, i_endidx



          ! departure region and correct counterclockwise numbering of vertices
          !--------------------------------------------------------------------
          !
          !  3\--------------\2       <- correct counterclockwise numbering
          !    \              \
          !     \      N       \      <- edge normal
          !      \      \       \
          !     4 \------\-------\1   <- triangle edge


          ! Determine the upwind cell
          ! Cell indices are chosen such that the direction from cell 1 to cell 2
          ! is the positive direction of the normal vector N.

          !
          ! Calculate backward trajectories, starting at the two edge vertices
          ! It is assumed the velocity vector is constant along the edge.
          !

          ! logical switch for MERGE operations: .TRUE. for p_vn >= 0
          lvn_pos = p_vn(je,jk,jb) >= 0._wp

          ! line and block indices of upwind cell
          p_cell_indices(je,jk,jb,1) = MERGE(ptr_p%edges%cell_idx(je,jb,1),       &
            &                                ptr_p%edges%cell_idx(je,jb,2),lvn_pos)
          p_cell_indices(je,jk,jb,2) = MERGE(ptr_p%edges%cell_blk(je,jb,1),       &
            &                                ptr_p%edges%cell_blk(je,jb,2),lvn_pos)


          ! Vertices 3 and 4 of the departure cell
          ! Take care of correct counterclockwise numbering below

          ! position of vertex 3 in normal direction
          pos_dreg_vert_e(3,1) = ptr_ve(je,jb,2,1) - p_vn(je,jk,jb) * p_dt

          ! position of vertex 3 in tangential direction
          pos_dreg_vert_e(3,2) = ptr_ve(je,jb,2,2) - p_vt(je,jk,jb) * p_dt

          ! position of vertex 4 (vn > 0) / vertex 2(vn < 0) in normal direction
          pos_dreg_vert_e(4,1) = ptr_ve(je,jb,1,1) - p_vn(je,jk,jb) * p_dt

          ! position of vertex 4 (vn > 0) / vertex 2(vn < 0) in tangential direction
          pos_dreg_vert_e(4,2) = ptr_ve(je,jb,1,2) - p_vt(je,jk,jb) * p_dt


          ! determine correct position on tangential plane
          pos_on_tplane_e(1:2) = MERGE(ptr_int%pos_on_tplane_e(je,jb,1,1:2), &
            &                          ptr_int%pos_on_tplane_e(je,jb,2,1:2),lvn_pos)


          ! Calculate position of departure region vertices in a translated
          ! coordinate system. The origin is located at the circumcenter
          ! of the upwind cell. The distance vectors point from the cell center
          ! to the vertices.
          pos_dreg_vert_c(1,1:2) = ptr_ve(je,jb,1,1:2) - pos_on_tplane_e(1:2)

          pos_dreg_vert_c(2,1:2) =                                       &
            & MERGE(ptr_ve(je,jb,2,1:2),pos_dreg_vert_e(4,1:2),lvn_pos)  &
            &     - pos_on_tplane_e(1:2)

          pos_dreg_vert_c(3,1:2) = pos_dreg_vert_e(3,1:2) - pos_on_tplane_e(1:2)

          pos_dreg_vert_c(4,1:2) =                                       &
            & MERGE(pos_dreg_vert_e(4,1:2),ptr_ve(je,jb,2,1:2),lvn_pos)  &
            &     - pos_on_tplane_e(1:2)


          ! In a last step, these distance vectors are transformed into a rotated
          ! geographical coordinate system, which still has its origin at the circumcenter
          ! of the upwind cell. Now the coordinate axes point to local East and local
          ! North.

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
          p_coords_dreg_v(je,jk,jb,1,1) =                                         &
            & pos_dreg_vert_c(1,1) * pn_cell(1) + pos_dreg_vert_c(1,2) * dn_cell(1)

          p_coords_dreg_v(je,jk,jb,2,1) =                                         &
            & pos_dreg_vert_c(2,1) * pn_cell(1) + pos_dreg_vert_c(2,2) * dn_cell(1)

           p_coords_dreg_v(je,jk,jb,3,1) =                                        &
            & pos_dreg_vert_c(3,1) * pn_cell(1) + pos_dreg_vert_c(3,2) * dn_cell(1)

          p_coords_dreg_v(je,jk,jb,4,1) =                                         &
            & pos_dreg_vert_c(4,1) * pn_cell(1) + pos_dreg_vert_c(4,2) * dn_cell(1)


          ! components in latitudinal direction
          p_coords_dreg_v(je,jk,jb,1,2) =                                         &
            & pos_dreg_vert_c(1,1) * pn_cell(2) + pos_dreg_vert_c(1,2) * dn_cell(2)

          p_coords_dreg_v(je,jk,jb,2,2) =                                         &
            & pos_dreg_vert_c(2,1) * pn_cell(2) + pos_dreg_vert_c(2,2) * dn_cell(2)

          p_coords_dreg_v(je,jk,jb,3,2) =                                         &
            & pos_dreg_vert_c(3,1) * pn_cell(2) + pos_dreg_vert_c(3,2) * dn_cell(2)

          p_coords_dreg_v(je,jk,jb,4,2) =                                         &
            & pos_dreg_vert_c(4,1) * pn_cell(2) + pos_dreg_vert_c(4,2) * dn_cell(2)


        ENDDO ! loop over edges
      ENDDO   ! loop over vertical levels
    END DO    ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE back_traj_dreg_o1



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
  SUBROUTINE back_traj_o2( ptr_p, ptr_int, p_vn, p_vt, p_dthalf, p_cell_indices,  &
    &                     p_distv_bary, opt_rlstart, opt_rlend, opt_slev,         &
    &                     opt_elev )

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
!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx)
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
!$OMP            vn_new,vt_new,z_ntdistv_bary)
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
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE back_traj_o2




  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Lax Friedrichs first order upwind flux,.
  !!
  !! Lax Friedrichs first order upwind flux,
  !! used in conservative advection routines.
  !! For passive advection, equivalent to
  !! any other first order upwind flux.
  !!
  !! @par Revision History
  !! Developed  by L.Bonaventura  (2004).
  !!
  FUNCTION laxfr_upflux( p_vn, p_psi1, p_psi2 )  RESULT(p_upflux)
    !

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: p_vn
    REAL(wp), INTENT(in) :: p_psi1, p_psi2

    REAL(wp) :: p_upflux

    !-----------------------------------------------------------------------
    p_upflux = 0.5_wp * (        p_vn  *( p_psi1 + p_psi2 )    &
      &                   - ABS( p_vn )*( p_psi2 - p_psi1 ) )

  END FUNCTION laxfr_upflux



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Lax Friedrichs first order upwind flux for vertical advection,.
  !!
  !! Generalized Lax Friedrichs first order upwind flux,
  !! used in conservative vertical advection routines.
  !! For passive advection, equivalent to any other first
  !! order upwind flux.
  !! Applicable to both pressure based and height based vertical
  !! coordinate systems. Depending on the coordinate system chosen,
  !! the sign of the second term in the flux equation changes.
  !! - (-) for pressure based vertical coordinate systems
  !! - (+) for height based coordinate systems
  !! In order to get the correct sign, the variable p_coeff_grid
  !! has been introduced which is =1 for pressure based and =-1
  !! for height based coordinate systems.
  !!
  !! @par Revision History
  !! Developed  by L.Bonaventura  (2004).
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized for p- and z-based vertical coordinate systems
  !!
  FUNCTION laxfr_upflux_v( p_vn, p_psi1, p_psi2, p_coeff_grid )  RESULT(p_upflux)
    !

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: p_vn
    REAL(wp), INTENT(in) :: p_psi1, p_psi2
    REAL(wp), INTENT(in) :: p_coeff_grid

    REAL(wp) :: p_upflux

    !-----------------------------------------------------------------------
    p_upflux = 0.5_wp * (                       p_vn  *( p_psi1 + p_psi2 )    &
      &                   - p_coeff_grid * ABS( p_vn )*( p_psi2 - p_psi1 ) )

  END FUNCTION laxfr_upflux_v




  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of quadratic tracer subgrid distribution
  !!
  !! Provides tracer-independent parts for a gauss quadrature of order 2.
  !! I.e. the 4 quadrature points in physical space and the product of weights
  !! and the determinant of the Jacobian for each quadrature point.
  !! This subroutine is specific to a reconstruction based on a quadratic
  !! polynomial. It needs to be called only once per time step, independent
  !! of the number of advected fields.
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-05-14)
  !!
  !!
  SUBROUTINE prep_gauss_quadrature_q( p_patch, p_coords_dreg_v,         &
    &                                 p_quad_vector_sum, p_rdreg_area,  &
    &                                 opt_rlstart, opt_rlend, opt_slev, &
    &                                 opt_elev )

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                           !< performed

    REAL(wp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:,:) !< in 2D cartesian coordinates
                                    !< dim: (nproma,nlev,nblks_e,4,2)

    REAL(wp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:) !< dim: (nproma,6,nlev,nblks_e)

    REAL(wp), INTENT(OUT) :: &      !< reciprocal area of departure region
      &  p_rdreg_area(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

   ! local variables
    REAL(wp) ::                &    !< coordinates of gaussian quadrature points
      &  z_gauss_pts(4,2)           !< in physical space

    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac(4)            !< each gaussian quadrature point.

    REAL(wp) ::                &    !< quadrature vector for single integration point
      &  z_quad_vector(4,6)

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: nlev                !< number of full levels
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: slev, elev          !< vertical start and end level

  !-----------------------------------------------------------------------

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

    ! number of vertical levels
    nlev = p_patch%nlev

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,jg,i_startidx,i_endidx,z_gauss_pts,wgt_t_detjac,z_quad_vector)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = slev, elev

        DO je = i_startidx, i_endidx

          ! get coordinates of the quadrature points in physical space (mapping)
          z_gauss_pts(1,1) = DOT_PRODUCT(shape_func(1:4,1),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(1,2) = DOT_PRODUCT(shape_func(1:4,1),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(2,1) = DOT_PRODUCT(shape_func(1:4,2),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(2,2) = DOT_PRODUCT(shape_func(1:4,2),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(3,1) = DOT_PRODUCT(shape_func(1:4,3),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(3,2) = DOT_PRODUCT(shape_func(1:4,3),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(4,1) = DOT_PRODUCT(shape_func(1:4,4),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(4,2) = DOT_PRODUCT(shape_func(1:4,4),p_coords_dreg_v(je,jk,jb,1:4,2))


          ! get Jacobian determinant for each quadrature point and multiply with
          ! corresponding weights
          ! Note: dbl_eps is added, in order to have a meaningful 'edge value' 
          ! (better: area-average) even when the integration-area tends to zero.
          wgt_t_detjac(1) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(1),eta(1)) &
            &                      * wgt_zeta(1) *  wgt_eta(1) ) + dbl_eps
          wgt_t_detjac(2) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(2),eta(2)) &
            &                      * wgt_zeta(2) *  wgt_eta(2) ) + dbl_eps
          wgt_t_detjac(3) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(3),eta(3)) &
            &                      * wgt_zeta(3) *  wgt_eta(3) ) + dbl_eps
          wgt_t_detjac(4) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(4),eta(4)) &
            &                      * wgt_zeta(4) *  wgt_eta(4) ) + dbl_eps


          ! Get quadrature vector for each integration point and multiply by
          ! corresponding wgt_t_detjac
          DO jg=1, 4
            z_quad_vector(jg,1) = wgt_t_detjac(jg)
            z_quad_vector(jg,2) = wgt_t_detjac(jg) * z_gauss_pts(jg,1)
            z_quad_vector(jg,3) = wgt_t_detjac(jg) * z_gauss_pts(jg,2)
            z_quad_vector(jg,4) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1)**2)
            z_quad_vector(jg,5) = wgt_t_detjac(jg) * (z_gauss_pts(jg,2)**2)
            z_quad_vector(jg,6) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1) * z_gauss_pts(jg,2))
          ENDDO


          ! Sum quadrature vectors over all integration points
          p_quad_vector_sum(je,1,jk,jb) = SUM(z_quad_vector(:,1))
          p_quad_vector_sum(je,2,jk,jb) = SUM(z_quad_vector(:,2))
          p_quad_vector_sum(je,3,jk,jb) = SUM(z_quad_vector(:,3))
          p_quad_vector_sum(je,4,jk,jb) = SUM(z_quad_vector(:,4))
          p_quad_vector_sum(je,5,jk,jb) = SUM(z_quad_vector(:,5))
          p_quad_vector_sum(je,6,jk,jb) = SUM(z_quad_vector(:,6))


          ! reciprocal area of departure region
          p_rdreg_area(je,jk,jb) = 1._wp/SUM(wgt_t_detjac(1:4))

        ENDDO ! loop over edges

      ENDDO  ! loop over levels

    ENDDO  ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE prep_gauss_quadrature_q


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of cubic tracer subgrid distribution
  !!
  !! Provides tracer-independent parts for a gauss quadrature of order 2.
  !! I.e. the 4 quadrature points in physical space and the product of weights
  !! and the determinant of the Jacobian for each quadrature point.
  !! This subroutine is specific to a cubic reconstruction without third order
  !! cross derivatives.
  !! It needs to be called only once per time step, independent of the number
  !! of advected fields.
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-10-13)
  !!
  !!
  SUBROUTINE prep_gauss_quadrature_cpoor( p_patch, p_coords_dreg_v,         &
    &                                     p_quad_vector_sum, p_rdreg_area,  &
    &                                     opt_rlstart, opt_rlend, opt_slev, &
    &                                     opt_elev                          )

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(IN) :: & !< patch on which computation is
      &  p_patch                           !< performed

    REAL(wp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:,:) !< in 2D cartesian coordinates
                                    !< dim: (nproma,nlev,nblks_e,4,2)

    REAL(wp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:) !< dim: (nproma,8,nlev,nblks_e)

    REAL(wp), INTENT(OUT) :: &      !< reciprocal area of departure region
      &  p_rdreg_area(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

   ! local variables
    REAL(wp) ::                &    !< coordinates of gaussian quadrature points
      &  z_gauss_pts(4,2)           !< in physical space

    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac(4)            !< each gaussian quadrature point.

    REAL(wp) ::                &    !< quadrature vector for single integration point
      &  z_quad_vector(4,8)

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: nlev                !< number of full levels
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: slev, elev          !< vertical start and end level

  !-----------------------------------------------------------------------

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

    ! number of vertical levels
    nlev = p_patch%nlev

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,jg,i_startidx,i_endidx,z_gauss_pts,wgt_t_detjac,z_quad_vector)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = slev, elev

        DO je = i_startidx, i_endidx

          ! get coordinates of the quadrature points in physical space (mapping)
          z_gauss_pts(1,1) = DOT_PRODUCT(shape_func(1:4,1),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(1,2) = DOT_PRODUCT(shape_func(1:4,1),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(2,1) = DOT_PRODUCT(shape_func(1:4,2),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(2,2) = DOT_PRODUCT(shape_func(1:4,2),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(3,1) = DOT_PRODUCT(shape_func(1:4,3),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(3,2) = DOT_PRODUCT(shape_func(1:4,3),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(4,1) = DOT_PRODUCT(shape_func(1:4,4),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(4,2) = DOT_PRODUCT(shape_func(1:4,4),p_coords_dreg_v(je,jk,jb,1:4,2))


          ! get Jacobian determinant for each quadrature point and multiply with
          ! corresponding weights
          ! Note: dbl_eps is added, in order to have a meaningful 'edge value' 
          ! (better: area-average) even when the integration-area tends to zero.
          wgt_t_detjac(1) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(1),eta(1)) &
            &                      * wgt_zeta(1) *  wgt_eta(1) ) + dbl_eps
          wgt_t_detjac(2) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(2),eta(2)) &
            &                      * wgt_zeta(1) *  wgt_eta(2) ) + dbl_eps
          wgt_t_detjac(3) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(3),eta(3)) &
            &                      * wgt_zeta(2) *  wgt_eta(1) ) + dbl_eps
          wgt_t_detjac(4) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(4),eta(4)) &
            &                      * wgt_zeta(2) *  wgt_eta(2) ) + dbl_eps


          ! Get quadrature vector for each integration point and multiply by
          ! corresponding wgt_t_detjac
          DO jg=1, 4
            z_quad_vector(jg,1) = wgt_t_detjac(jg)
            z_quad_vector(jg,2) = wgt_t_detjac(jg) * z_gauss_pts(jg,1)
            z_quad_vector(jg,3) = wgt_t_detjac(jg) * z_gauss_pts(jg,2)
            z_quad_vector(jg,4) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1)**2)
            z_quad_vector(jg,5) = wgt_t_detjac(jg) * (z_gauss_pts(jg,2)**2)
            z_quad_vector(jg,6) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1) * z_gauss_pts(jg,2))
            z_quad_vector(jg,7) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1)**3)
            z_quad_vector(jg,8) = wgt_t_detjac(jg) * (z_gauss_pts(jg,2)**3)
          ENDDO


          ! Sum quadrature vectors over all integration points
          p_quad_vector_sum(je,1,jk,jb) = SUM(z_quad_vector(:,1))
          p_quad_vector_sum(je,2,jk,jb) = SUM(z_quad_vector(:,2))
          p_quad_vector_sum(je,3,jk,jb) = SUM(z_quad_vector(:,3))
          p_quad_vector_sum(je,4,jk,jb) = SUM(z_quad_vector(:,4))
          p_quad_vector_sum(je,5,jk,jb) = SUM(z_quad_vector(:,5))
          p_quad_vector_sum(je,6,jk,jb) = SUM(z_quad_vector(:,6))
          p_quad_vector_sum(je,7,jk,jb) = SUM(z_quad_vector(:,7))
          p_quad_vector_sum(je,8,jk,jb) = SUM(z_quad_vector(:,8))


          ! reciprocal area of departure region
          p_rdreg_area(je,jk,jb) = 1._wp/SUM(wgt_t_detjac(1:4))

        ENDDO ! loop over edges

      ENDDO  ! loop over levels

    ENDDO  ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE prep_gauss_quadrature_cpoor


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of cubic tracer subgrid distribution
  !!
  !! Provides tracer-independent parts for a gauss quadrature of order 2.
  !! I.e. the 4 quadrature points in physical space and the product of weights
  !! and the determinant of the Jacobian for each quadrature point.
  !! This subroutine is specific to a cubic reconstruction.
  !! It needs to be called only once per time step, independent of the number
  !! of advected fields.
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-10-13)
  !!
  !!
  SUBROUTINE prep_gauss_quadrature_c( p_patch, p_coords_dreg_v,         &
    &                                 p_quad_vector_sum, p_rdreg_area,  &
    &                                 opt_rlstart, opt_rlend, opt_slev, &
    &                                 opt_elev                          )

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                           !< performed

    REAL(wp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:,:) !< in 2D cartesian coordinates
                                    !< dim: (nproma,nlev,nblks_e,4,2)

    REAL(wp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:) !< dim: (nproma,8,nlev,nblks_e)

    REAL(wp), INTENT(OUT) :: &      !< reciprocal area of departure region
      &  p_rdreg_area(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

   ! local variables
    REAL(wp) ::                &    !< coordinates of gaussian quadrature points
      &  z_gauss_pts(4,2)           !< in physical space

    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac(4)            !< each gaussian quadrature point.

    REAL(wp) ::                &    !< quadrature vector for single integration point
      &  z_quad_vector(4,10)

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: nlev                !< number of full levels
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: slev, elev          !< vertical start and end level

  !-----------------------------------------------------------------------

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

    ! number of vertical levels
    nlev = p_patch%nlev

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,jg,i_startidx,i_endidx,z_gauss_pts,wgt_t_detjac,z_quad_vector)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = slev, elev

        DO je = i_startidx, i_endidx

          ! get coordinates of the quadrature points in physical space (mapping)
          z_gauss_pts(1,1) = DOT_PRODUCT(shape_func(1:4,1),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(1,2) = DOT_PRODUCT(shape_func(1:4,1),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(2,1) = DOT_PRODUCT(shape_func(1:4,2),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(2,2) = DOT_PRODUCT(shape_func(1:4,2),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(3,1) = DOT_PRODUCT(shape_func(1:4,3),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(3,2) = DOT_PRODUCT(shape_func(1:4,3),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(4,1) = DOT_PRODUCT(shape_func(1:4,4),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(4,2) = DOT_PRODUCT(shape_func(1:4,4),p_coords_dreg_v(je,jk,jb,1:4,2))


          ! get Jacobian determinant for each quadrature point and multiply with
          ! corresponding weights
          ! Note: dbl_eps is added, in order to have a meaningful 'edge value' 
          ! (better: area-average) even when the integration-area tends to zero.
          wgt_t_detjac(1) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(1),eta(1)) &
            &                      * wgt_zeta(1) *  wgt_eta(1) ) + dbl_eps
          wgt_t_detjac(2) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(2),eta(2)) &
            &                      * wgt_zeta(1) *  wgt_eta(2) ) + dbl_eps
          wgt_t_detjac(3) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(3),eta(3)) &
            &                      * wgt_zeta(2) *  wgt_eta(1) ) + dbl_eps
          wgt_t_detjac(4) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(4),eta(4)) &
            &                      * wgt_zeta(2) *  wgt_eta(2) ) + dbl_eps


          ! Get quadrature vector for each integration point and multiply by
          ! corresponding wgt_t_detjac
          DO jg=1, 4
            z_quad_vector(jg,1) = wgt_t_detjac(jg)
            z_quad_vector(jg,2) = wgt_t_detjac(jg) * z_gauss_pts(jg,1)
            z_quad_vector(jg,3) = wgt_t_detjac(jg) * z_gauss_pts(jg,2)
            z_quad_vector(jg,4) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1)**2)
            z_quad_vector(jg,5) = wgt_t_detjac(jg) * (z_gauss_pts(jg,2)**2)
            z_quad_vector(jg,6) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1) * z_gauss_pts(jg,2))
            z_quad_vector(jg,7) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1)**3)
            z_quad_vector(jg,8) = wgt_t_detjac(jg) * (z_gauss_pts(jg,2)**3)
            z_quad_vector(jg,9) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1)**2 * z_gauss_pts(jg,2))
            z_quad_vector(jg,10)= wgt_t_detjac(jg) * (z_gauss_pts(jg,1) * z_gauss_pts(jg,2)**2)
          ENDDO


          ! Sum quadrature vectors over all integration points
          p_quad_vector_sum(je,1,jk,jb) = SUM(z_quad_vector(:,1))
          p_quad_vector_sum(je,2,jk,jb) = SUM(z_quad_vector(:,2))
          p_quad_vector_sum(je,3,jk,jb) = SUM(z_quad_vector(:,3))
          p_quad_vector_sum(je,4,jk,jb) = SUM(z_quad_vector(:,4))
          p_quad_vector_sum(je,5,jk,jb) = SUM(z_quad_vector(:,5))
          p_quad_vector_sum(je,6,jk,jb) = SUM(z_quad_vector(:,6))
          p_quad_vector_sum(je,7,jk,jb) = SUM(z_quad_vector(:,7))
          p_quad_vector_sum(je,8,jk,jb) = SUM(z_quad_vector(:,8))
          p_quad_vector_sum(je,9,jk,jb) = SUM(z_quad_vector(:,9))
          p_quad_vector_sum(je,10,jk,jb)= SUM(z_quad_vector(:,10))


          ! reciprocal area of departure region
          p_rdreg_area(je,jk,jb) = 1._wp/SUM(wgt_t_detjac(1:4))

        ENDDO ! loop over edges

      ENDDO  ! loop over levels

    ENDDO  ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE prep_gauss_quadrature_c



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Computes Jacobian determinant for gaussian quadrature
  !!
  !! Computes Jacobian determinant for gaussian quadrature
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-05-14)
  !!
  !!
  FUNCTION jac(x, y, zeta, eta)  RESULT(det_jac)

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: x(1:4), y(1:4)  !< coordinates of vertices in x-y-system
    REAL(wp), INTENT(IN) :: zeta, eta       !< integration point in \zeta,\eta-system

    ! RETURN VALUE:
    REAL(wp) :: det_jac

    REAL(wp), DIMENSION(2,2) :: jacob

  !-----------------------------------------------------------------------

    jacob(1,1) = -(1._wp - eta) * x(1) + (1._wp - eta) * x(2)  &
      &        +  (1._wp + eta) * x(3) - (1._wp + eta) * x(4)
    jacob(1,2) = -(1._wp - eta) * y(1) + (1._wp - eta) * y(2)  &
      &        +  (1._wp + eta) * y(3) - (1._wp + eta) * y(4)
    jacob(2,1) = -(1._wp - zeta)* x(1) - (1._wp + zeta)* x(2)  &
      &        +  (1._wp + zeta)* x(3) + (1._wp - zeta)* x(4)
    jacob(2,2) = -(1._wp - zeta)* y(1) - (1._wp + zeta)* y(2)  &
      &        +  (1._wp + zeta)* y(3) + (1._wp - zeta)* y(4)

    det_jac = 0.0625_wp * (jacob(1,1)*jacob(2,2) - jacob(1,2)*jacob(2,1))

  END FUNCTION jac


END MODULE mo_advection_utils

