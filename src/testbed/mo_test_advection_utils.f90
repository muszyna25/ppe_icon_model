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
MODULE mo_test_advection_utils

  USE mo_advection_config,    ONLY: shape_func, zeta, eta, wgt_zeta, wgt_eta
  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_interpolation,       ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants,      ONLY: min_rlcell_int, min_rledge_int
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  ! In order to avoid circular dependencies these two pointers
  ! have been moved from mo_advection_stepping to this module.
  REAL(wp), POINTER ::  &
    &  ptr_delp_mc_now(:,:,:) => NULL() !< pointer to old layer thickness
                                        !< at cell center
  REAL(wp), POINTER ::  &
    &  ptr_delp_mc_new(:,:,:) => NULL() !< pointer to new layer thickness
                                        !< at cell center

  PUBLIC :: test_back_traj_o1

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
  SUBROUTINE test_back_traj_o1( ptr_p, ptr_int, p_vn, p_vt, p_dthalf, p_cell_indices, &
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
    
    INTEGER :: downwind_indices(SIZE(p_vn,1),SIZE(p_vn,2)), downwind_index

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

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,z_ntdistv_bary, &
!$OMP downwind_indices, downwind_index, upwind_index)

    DO jb = i_startblk, i_endblk

     ! first find the downwind index
     downwind_indices(:,:) = MERGE(1,2, p_vn(:,:,jb) >= 0._wp)
     
     CALL get_indices_e(ptr_p, jb, i_startblk, i_endblk,        &
                        i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = slev, elev
        DO je = i_startidx, i_endidx

           downwind_index = downwind_indices(je,jk)
!            upwind_index = 3 - downwind_index
          !
          ! Calculate backward trajectories
          !

          ! position of barycenter in normal direction
!           pos_barycenter(1) = - p_vn(je,jk,jb) * p_dthalf

          ! position of barycenter in tangential direction
!           pos_barycenter(2) = - p_vt(je,jk,jb) * p_dthalf

          ! logical auxiliary for MERGE operations: .TRUE. for vn >= 0
!           lvn_pos = p_vn(je,jk,jb) >= 0._wp

          ! If vn > 0 (vn < 0), the upwind cell is cell 1 (cell 2)

          ! line and block indices of neighbor cell with barycenter          
          p_cell_indices(je,jk,jb,1) = ptr_p%edges%cell_idx(je,jb,downwind_index)
          p_cell_indices(je,jk,jb,2) = ptr_p%edges%cell_blk(je,jb,downwind_index)

          ! Calculate the distance cell center --> barycenter for the cell,
          ! in which the barycenter is located. The distance vector points
          ! from the cell center to the barycenter.
!           z_ntdistv_bary(1:2) = pos_barycenter(1:2)               &
!             & - MERGE(ptr_int%pos_on_tplane_e(je,jb,1,1:2),       &
!             &         ptr_int%pos_on_tplane_e(je,jb,2,1:2),lvn_pos)
          ! position of barycenter in normal direction
           z_ntdistv_bary(1) =  - (p_vn(je,jk,jb) * p_dthalf  &
            & + ptr_int%pos_on_tplane_e(je,jb,downwind_index,1))
          ! position of barycenter in tangential direction
          z_ntdistv_bary(2) =  - p_vt(je,jk,jb) * p_dthalf  &
            & - ptr_int%pos_on_tplane_e(je,jb,downwind_index,2)


          ! In a last step, transform this distance vector into a rotated
          ! geographical coordinate system with its origin at the circumcenter
          ! of the upstream cell. Coordinate axes point to local East and local
          ! North.

          ! component in longitudinal direction
          p_distv_bary(je,jk,jb,1) =                                                        &
            &   z_ntdistv_bary(1) * ptr_p%edges%primal_normal_cell(je,jb,downwind_index)%v1 &
            & + z_ntdistv_bary(2) * ptr_p%edges%dual_normal_cell(je,jb,downwind_index)%v1

          ! component in latitudinal direction
          p_distv_bary(je,jk,jb,2) =                                                        &
            &   z_ntdistv_bary(1) * ptr_p%edges%primal_normal_cell(je,jb,downwind_index)%v2 &
            & + z_ntdistv_bary(2) * ptr_p%edges%dual_normal_cell(je,jb,downwind_index)%v2

        ENDDO ! loop over edges
      ENDDO   ! loop over vertical levels
    END DO    ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE test_back_traj_o1




END MODULE mo_test_advection_utils

