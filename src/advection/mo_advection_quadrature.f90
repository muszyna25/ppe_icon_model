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
!! Modification by Daniel Reinert, DWD (2011-05-03)
!! - added quadrature routine for integrating a linear polynomial
!! Modification by Will Sawyer, CSCS (2016-07-14)
!! - added OpenACC support
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
MODULE mo_advection_quadrature

  USE mo_kind,                ONLY: wp, vp
  USE mo_advection_config,    ONLY: shape_func, zeta, eta, wgt_zeta, wgt_eta, &
    &                               shape_func_l, zeta_l, eta_l, wgt_zeta_l,  &
    &                               wgt_eta_l
  USE mo_advection_utils,     ONLY: t_list2D
  USE mo_model_domain,        ONLY: t_patch
  USE mo_loopindices,         ONLY: get_indices_e
  USE mo_impl_constants,      ONLY: min_rledge_int
  USE mo_math_constants,      ONLY: dbl_eps, eps
  USE mo_parallel_config,     ONLY: nproma
#ifdef _OPENACC
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, check_patch_array
  USE mo_mpi,                 ONLY: i_am_accel_node
#endif

  IMPLICIT NONE

  PRIVATE



  PUBLIC :: prep_gauss_quadrature_l
  PUBLIC :: prep_gauss_quadrature_l_list
  PUBLIC :: prep_gauss_quadrature_q
  PUBLIC :: prep_gauss_quadrature_q_list
  PUBLIC :: prep_gauss_quadrature_c
  PUBLIC :: prep_gauss_quadrature_c_list
  PUBLIC :: prep_gauss_quadrature_cpoor

#if defined( _OPENACC )
#if defined(__ADVECTION_QUADRATURE_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
  LOGICAL, PARAMETER ::  acc_validate = .FALSE.   ! ONLY SET TO .TRUE. FOR VALIDATION PHASE
#endif

  
CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of linear tracer subgrid distribution
  !!
  !! Provides tracer-independent parts for a gauss quadrature of order 1.
  !! I.e. a single quadrature point in physical space and the product of weights
  !! and the determinant of the Jacobian for the quadrature point.
  !! This subroutine is specific to a linear polynomial. It needs to be called 
  !! only once per time step, independent of the number of advected fields.
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2012-05-03)
  !!
  !!
  SUBROUTINE prep_gauss_quadrature_l( p_patch, p_coords_dreg_v,         &
    &                                 p_quad_vector_sum, p_dreg_area,   &
    &                                 opt_rlstart, opt_rlend, opt_slev, &
    &                                 opt_elev )

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                             !< performed

    REAL(vp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:,:) !< in 2D cartesian coordinates
                                    !< dim: (nproma,4,2,nlev,nblks_e)

    REAL(vp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:) !< dim: (nproma,3,nlev,nblks_e)

    REAL(vp), INTENT(OUT) :: &      !< area of departure region  [m**2]
      &  p_dreg_area(:,:,:)         !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

   ! local variables
    REAL(wp) ::                &       !< coordinates of gaussian quadrature points
      &  z_gauss_pts_1, z_gauss_pts_2  !< in physical space

    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac               !< gaussian quadrature point.

    REAL(wp) :: z_x(nproma,4), z_y(nproma,4) !< storage for local coordinates

    INTEGER  :: jb, je, jk          !< loop index for blocks and edges, levels
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: slev, elev          !< vertical start and end level

  !-----------------------------------------------------------------------

!$ACC DATA PCOPYIN( p_coords_dreg_v), PCOPYOUT( p_quad_vector_sum, p_dreg_area ), CREATE( z_x, z_y ), &
!$ACC      IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_coords_dreg_v, p_dreg_area ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

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
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx,z_gauss_pts_1,z_gauss_pts_2,wgt_t_detjac,z_x,z_y &
!$OMP ) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG PRIVATE( z_x, z_y )
      DO jk = slev, elev

        !$ACC LOOP VECTOR
        DO je = i_startidx, i_endidx

          z_x(je,1:4) = p_coords_dreg_v(je,1:4,1,jk,jb)
          z_y(je,1:4) = p_coords_dreg_v(je,1:4,2,jk,jb)

          ! get coordinates of the quadrature points in physical space (mapping)
!WS: TODO:  make sure that DOT_PRODUCT is supported in this OpenACC context
          z_gauss_pts_1 = DOT_PRODUCT(shape_func_l(1:4),z_x(je,1:4))
          z_gauss_pts_2 = DOT_PRODUCT(shape_func_l(1:4),z_y(je,1:4))


          ! get Jacobian determinant for each quadrature point and multiply with
          ! corresponding weights
          ! Note: dbl_eps is added, in order to have a meaningful 'edge value' 
          ! (better: area-average) even when the integration-area tends to zero.
          wgt_t_detjac = ( jac(z_x(je,1:4),z_y(je,1:4),zeta_l,eta_l) &
            &                    * wgt_zeta_l * wgt_eta_l ) + dbl_eps


          ! Get quadrature vector for each integration point and multiply by
          ! corresponding wgt_t_detjac. No summation necessary, since a 
          ! single integration point is used.
          p_quad_vector_sum(je,1,jk,jb) = wgt_t_detjac
          p_quad_vector_sum(je,2,jk,jb) = wgt_t_detjac * z_gauss_pts_1
          p_quad_vector_sum(je,3,jk,jb) = wgt_t_detjac * z_gauss_pts_2

          ! area of departure region
          p_dreg_area(je,jk,jb) = wgt_t_detjac

        ENDDO ! loop over edges

      ENDDO  ! loop over levels
!$ACC END PARALLEL

    ENDDO  ! loop over blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( p_dreg_area, p_quad_vector_sum ), IF ( acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA

  END SUBROUTINE prep_gauss_quadrature_l


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of linear tracer subgrid distribution
  !!
  !! Provides tracer-independent parts for a gauss quadrature of order 1.
  !! I.e. a single quadrature point in physical space and the product of weights
  !! and the determinant of the Jacobian for the quadrature point.
  !! This subroutine is specific to a linear polynomial. It needs to be called 
  !! only once per time step, independent of the number of advected fields.
  !!
  !! Index-list based version. Otherwise identical to prep_gauss_quadrature_l
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2013-10-28)
  !!
  !!
  SUBROUTINE prep_gauss_quadrature_l_list( p_patch, p_coords_dreg_v, falist, &
    &                                 p_quad_vector_sum, p_dreg_area,        &
    &                                 opt_rlstart, opt_rlend, opt_slev,      &
    &                                 opt_elev )

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                             !< performed

    REAL(vp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:)   !< in 2D cartesian coordinates
                                    !< dim: (npoints,4,2,nblks_e)

    TYPE(t_list2D), INTENT(IN) :: & !< index list with points for which the standard 
      &  falist                     !< Miura-type treatment of flux areas is 
                                    !< insufficient

    REAL(vp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:)   !< dim: (npoints,3,nblks_e)

    REAL(vp), INTENT(INOUT) :: &    !< total area of departure region  [m**2]
      &  p_dreg_area(:,:,:)         !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

   ! local variables
    REAL(wp) ::                &       !< coordinates of gaussian quadrature points
      &  z_gauss_pts_1, z_gauss_pts_2  !< in physical space

    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac               !< gaussian quadrature point.

    REAL(wp) :: z_x(nproma,4), z_y(nproma,4) !< storage for local coordinates

    INTEGER  :: jb, je, jk          !< loop index for blocks and edges, levels
    INTEGER  :: ie                  !< index list loop counter
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: slev, elev          !< vertical start and end level

  !-----------------------------------------------------------------------

!$ACC DATA PCOPYIN( p_coords_dreg_v, falist), PCOPY( p_dreg_area ),   &
!$ACC      PCOPYOUT( p_quad_vector_sum ), CREATE( z_x, z_y ),         &
!$ACC      IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_coords_dreg_v, p_dreg_area ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

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
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG PRIVATE( z_x, z_y ) 

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,ie,z_gauss_pts_1,z_gauss_pts_2,wgt_t_detjac,z_x,z_y) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

!$ACC LOOP VECTOR
!CDIR NODEP,VOVERTAKE,VOB
      DO ie = 1, falist%len(jb)

        je = falist%eidx(ie,jb)
        jk = falist%elev(ie,jb)

        z_x(ie,1:4) = p_coords_dreg_v(ie,1:4,1,jb)
        z_y(ie,1:4) = p_coords_dreg_v(ie,1:4,2,jb)

        ! get coordinates of the quadrature points in physical space (mapping)
        z_gauss_pts_1 = DOT_PRODUCT(shape_func_l(1:4),z_x(ie,1:4))
        z_gauss_pts_2 = DOT_PRODUCT(shape_func_l(1:4),z_y(ie,1:4))


        ! get Jacobian determinant for each quadrature point and multiply with
        ! corresponding weights
        ! Note: dbl_eps is added, in order to have a meaningful 'edge value' 
        ! (better: area-average) even when the integration-area tends to zero.
        wgt_t_detjac = ( jac(z_x(ie,1:4),z_y(ie,1:4),zeta_l,eta_l) &
          &                     * wgt_zeta_l * wgt_eta_l ) + dbl_eps


        ! Get quadrature vector for each integration point and multiply by
        ! corresponding wgt_t_detjac. No summation necessary, since a 
        ! single integration point is used.
        p_quad_vector_sum(ie,1,jb) = wgt_t_detjac
        p_quad_vector_sum(ie,2,jb) = wgt_t_detjac * z_gauss_pts_1
        p_quad_vector_sum(ie,3,jb) = wgt_t_detjac * z_gauss_pts_2

        ! Add contribution to total area of departure region
        p_dreg_area(je,jk,jb) = p_dreg_area(je,jk,jb) + wgt_t_detjac

      ENDDO ! ie: loop over index list

    ENDDO  ! loop over blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC END PARALLEL

!$ACC UPDATE HOST( p_dreg_area, p_quad_vector_sum ), IF ( acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA

  END SUBROUTINE prep_gauss_quadrature_l_list

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
    &                                 p_quad_vector_sum, p_dreg_area,   &
    &                                 opt_rlstart, opt_rlend, opt_slev, &
    &                                 opt_elev )

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                           !< performed

    REAL(vp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:,:) !< in 2D cartesian coordinates
                                    !< dim: (nproma,4,2,nlev,nblks_e)

    REAL(vp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:) !< dim: (nproma,6,nlev,nblks_e)

    REAL(vp), INTENT(OUT) :: &      !< area of departure region  [m**2]
      &  p_dreg_area(:,:,:)         !< dim: (nproma,nlev,nblks_e)

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
      &  z_gauss_pts(nproma,4,2)    !< in physical space

    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac(nproma,4)     !< each gaussian quadrature point.

    REAL(wp) ::                &    !< quadrature vector for single integration point
      &  z_quad_vector(nproma,4,6)

    REAL(wp) :: z_x(nproma,4), z_y(nproma,4) !< storage for local coordinates
    REAL(wp) :: z_wgt(4), z_eta(4,4)         !< for precomputation of coefficients
    REAL(wp) :: z_area                       !< auxiliary for dreg area

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: slev, elev          !< vertical start and end level

  !-----------------------------------------------------------------------

!$ACC DATA PCOPYIN( p_coords_dreg_v ), PCOPYOUT( p_quad_vector_sum, p_dreg_area ), &
!$ACC      CREATE( z_x, z_y, z_wgt, z_eta ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_coords_dreg_v, p_dreg_area ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

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
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

    z_wgt(1) = 0.0625_wp * wgt_zeta(1) *  wgt_eta(1)
    z_wgt(2) = 0.0625_wp * wgt_zeta(2) *  wgt_eta(2)
    z_wgt(3) = 0.0625_wp * wgt_zeta(3) *  wgt_eta(3)
    z_wgt(4) = 0.0625_wp * wgt_zeta(4) *  wgt_eta(4)

    z_eta(1,1:4) = 1._wp - eta(1:4)
    z_eta(2,1:4) = 1._wp + eta(1:4)
    z_eta(3,1:4) = 1._wp - zeta(1:4)
    z_eta(4,1:4) = 1._wp + zeta(1:4)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,jg,i_startidx,i_endidx,z_gauss_pts,wgt_t_detjac, &
!$OMP z_quad_vector,z_x,z_y,z_area) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG PRIVATE( z_x, z_y, wgt_t_detjac, z_gauss_pts, z_quad_vector )
      DO jk = slev, elev

        !$ACC LOOP VECTOR
        DO je = i_startidx, i_endidx

          z_x(je,1:4) = p_coords_dreg_v(je,1:4,1,jk,jb)
          z_y(je,1:4) = p_coords_dreg_v(je,1:4,2,jk,jb)

          ! get Jacobian determinant for each quadrature point and multiply with
          ! corresponding weights
          ! Note: dbl_eps is added, in order to have a meaningful 'edge value' 
          ! (better: area-average) even when the integration-area tends to zero.
          wgt_t_detjac(je,1:4) = dbl_eps + z_wgt(1:4) * (                                 &
            &   (z_eta(1,1:4)*(z_x(je,2)-z_x(je,1)) + z_eta(2,1:4)*(z_x(je,3)-z_x(je,4))) &
            & * (z_eta(3,1:4)*(z_y(je,4)-z_y(je,1)) - z_eta(4,1:4)*(z_y(je,2)-z_y(je,3))) &
            & - (z_eta(1,1:4)*(z_y(je,2)-z_y(je,1)) + z_eta(2,1:4)*(z_y(je,3)-z_y(je,4))) &
            & * (z_eta(3,1:4)*(z_x(je,4)-z_x(je,1)) - z_eta(4,1:4)*(z_x(je,2)-z_x(je,3))) )


          ! get coordinates of the quadrature points in physical space (mapping)
          z_gauss_pts(je,1,1) = DOT_PRODUCT(shape_func(1:4,1),z_x(je,1:4))
          z_gauss_pts(je,1,2) = DOT_PRODUCT(shape_func(1:4,1),z_y(je,1:4))
          z_gauss_pts(je,2,1) = DOT_PRODUCT(shape_func(1:4,2),z_x(je,1:4))
          z_gauss_pts(je,2,2) = DOT_PRODUCT(shape_func(1:4,2),z_y(je,1:4))
          z_gauss_pts(je,3,1) = DOT_PRODUCT(shape_func(1:4,3),z_x(je,1:4))
          z_gauss_pts(je,3,2) = DOT_PRODUCT(shape_func(1:4,3),z_y(je,1:4))
          z_gauss_pts(je,4,1) = DOT_PRODUCT(shape_func(1:4,4),z_x(je,1:4))
          z_gauss_pts(je,4,2) = DOT_PRODUCT(shape_func(1:4,4),z_y(je,1:4))

          ! Get quadrature vector for each integration point and multiply by
          ! corresponding wgt_t_detjac
          DO jg=1, 4
            z_quad_vector(je,jg,1) = wgt_t_detjac(je,jg)
            z_quad_vector(je,jg,2) = wgt_t_detjac(je,jg) * z_gauss_pts(je,jg,1)
            z_quad_vector(je,jg,3) = wgt_t_detjac(je,jg) * z_gauss_pts(je,jg,2)
            z_quad_vector(je,jg,4) = wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,1)**2)
            z_quad_vector(je,jg,5) = wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,2)**2)
            z_quad_vector(je,jg,6) = wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,1) * z_gauss_pts(je,jg,2))
          ENDDO


          ! Sum quadrature vectors over all integration points
          p_quad_vector_sum(je,1,jk,jb) = SUM(z_quad_vector(je,:,1))
          p_quad_vector_sum(je,2,jk,jb) = SUM(z_quad_vector(je,:,2))
          p_quad_vector_sum(je,3,jk,jb) = SUM(z_quad_vector(je,:,3))
          p_quad_vector_sum(je,4,jk,jb) = SUM(z_quad_vector(je,:,4))
          p_quad_vector_sum(je,5,jk,jb) = SUM(z_quad_vector(je,:,5))
          p_quad_vector_sum(je,6,jk,jb) = SUM(z_quad_vector(je,:,6))


          ! area of departure region
          z_area = SUM(wgt_t_detjac(je,1:4))
          p_dreg_area(je,jk,jb) = SIGN(MAX(eps,ABS(z_area)),z_area)

        ENDDO ! loop over edges

      ENDDO  ! loop over levels
!$ACC END PARALLEL

    ENDDO  ! loop over blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( p_dreg_area, p_quad_vector_sum ), IF ( acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA

  END SUBROUTINE prep_gauss_quadrature_q


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
  !! Index-list based version. Otherwise identical to prep_gauss_quadrature_q
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2013-10-28)
  !!
  !!
  SUBROUTINE prep_gauss_quadrature_q_list( p_patch, p_coords_dreg_v, falist, &
    &                                 p_quad_vector_sum, p_dreg_area,   &
    &                                 opt_rlstart, opt_rlend, opt_slev, &
    &                                 opt_elev )

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                           !< performed

    REAL(vp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:)   !< in 2D cartesian coordinates
                                    !< dim: (npoints,4,2,nblks_e)

    TYPE(t_list2D), INTENT(IN) :: & !< index list with points for which the standard 
      &  falist                     !< Miura-type treatment of flux areas is
                                    !< insufficient

    REAL(vp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:)   !< dim: (npoints,6,nblks_e)

    REAL(vp), INTENT(INOUT) :: &    !< total area of departure region  [m**2]
      &  p_dreg_area(:,:,:)         !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

   ! local variables
    REAL(wp) ::                        &    !< coordinates of gaussian quadrature points
      &  z_gauss_pts(falist%npoints,4,2)    !< in physical space

    REAL(wp) ::                       &     !< weights times determinant of Jacobian for
      &  wgt_t_detjac(falist%npoints,4)     !< each gaussian quadrature point.

    REAL(wp) ::                          &  !< quadrature vector for single integration point
      &  z_quad_vector(falist%npoints,4,6)

    REAL(wp) :: z_x(falist%npoints,4), z_y(falist%npoints,4) !< storage for local coordinates
    REAL(wp) :: z_wgt(4), z_eta(4,4)         !< for precomputation of coefficients

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: ie                  !< index list loop counter
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: slev, elev          !< vertical start and end level

  !-----------------------------------------------------------------------

!$ACC DATA PCOPYIN( p_coords_dreg_v, falist ), PCOPY( p_dreg_area), PCOPYOUT( p_quad_vector_sum ), &
!$ACC      CREATE( z_x, z_y, z_wgt, z_eta ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_coords_dreg_v, p_dreg_area ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

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
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

    z_wgt(1) = 0.0625_wp * wgt_zeta(1) *  wgt_eta(1)
    z_wgt(2) = 0.0625_wp * wgt_zeta(2) *  wgt_eta(2)
    z_wgt(3) = 0.0625_wp * wgt_zeta(3) *  wgt_eta(3)
    z_wgt(4) = 0.0625_wp * wgt_zeta(4) *  wgt_eta(4)

    z_eta(1,1:4) = 1._wp - eta(1:4)
    z_eta(2,1:4) = 1._wp + eta(1:4)
    z_eta(3,1:4) = 1._wp - zeta(1:4)
    z_eta(4,1:4) = 1._wp + zeta(1:4)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,ie,jg,z_gauss_pts,wgt_t_detjac, &
!$OMP z_quad_vector,z_x,z_y) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG VECTOR
      DO ie = 1, falist%len(jb)

        z_x(ie,1:4) = p_coords_dreg_v(ie,1:4,1,jb)
        z_y(ie,1:4) = p_coords_dreg_v(ie,1:4,2,jb)

        ! get Jacobian determinant for each quadrature point and multiply with
        ! corresponding weights
        ! Note: dbl_eps is added, in order to have a meaningful 'edge value' 
        ! (better: area-average) even when the integration-area tends to zero.
        wgt_t_detjac(ie,1:4) = dbl_eps + z_wgt(1:4) * (                                 &
          &   (z_eta(1,1:4)*(z_x(ie,2)-z_x(ie,1)) + z_eta(2,1:4)*(z_x(ie,3)-z_x(ie,4))) &
          & * (z_eta(3,1:4)*(z_y(ie,4)-z_y(ie,1)) - z_eta(4,1:4)*(z_y(ie,2)-z_y(ie,3))) &
          & - (z_eta(1,1:4)*(z_y(ie,2)-z_y(ie,1)) + z_eta(2,1:4)*(z_y(ie,3)-z_y(ie,4))) &
          & * (z_eta(3,1:4)*(z_x(ie,4)-z_x(ie,1)) - z_eta(4,1:4)*(z_x(ie,2)-z_x(ie,3))) )


        ! get coordinates of the quadrature points in physical space (mapping)
!WS: TODO:  check whether intrinsic DOT_PRODUCT is properly supported on GPU
        z_gauss_pts(ie,1,1) = DOT_PRODUCT(shape_func(1:4,1),z_x(ie,1:4))
        z_gauss_pts(ie,1,2) = DOT_PRODUCT(shape_func(1:4,1),z_y(ie,1:4))
        z_gauss_pts(ie,2,1) = DOT_PRODUCT(shape_func(1:4,2),z_x(ie,1:4))
        z_gauss_pts(ie,2,2) = DOT_PRODUCT(shape_func(1:4,2),z_y(ie,1:4))
        z_gauss_pts(ie,3,1) = DOT_PRODUCT(shape_func(1:4,3),z_x(ie,1:4))
        z_gauss_pts(ie,3,2) = DOT_PRODUCT(shape_func(1:4,3),z_y(ie,1:4))
        z_gauss_pts(ie,4,1) = DOT_PRODUCT(shape_func(1:4,4),z_x(ie,1:4))
        z_gauss_pts(ie,4,2) = DOT_PRODUCT(shape_func(1:4,4),z_y(ie,1:4))


        ! Get quadrature vector for each integration point and multiply by
        ! corresponding wgt_t_detjac
        DO jg=1, 4
          z_quad_vector(ie,jg,1) = wgt_t_detjac(ie,jg)
          z_quad_vector(ie,jg,2) = wgt_t_detjac(ie,jg) * z_gauss_pts(ie,jg,1)
          z_quad_vector(ie,jg,3) = wgt_t_detjac(ie,jg) * z_gauss_pts(ie,jg,2)
          z_quad_vector(ie,jg,4) = wgt_t_detjac(ie,jg) * (z_gauss_pts(ie,jg,1)**2)
          z_quad_vector(ie,jg,5) = wgt_t_detjac(ie,jg) * (z_gauss_pts(ie,jg,2)**2)
          z_quad_vector(ie,jg,6) = wgt_t_detjac(ie,jg) * (z_gauss_pts(ie,jg,1) * z_gauss_pts(ie,jg,2))
        ENDDO


        ! Sum quadrature vectors over all integration points
        p_quad_vector_sum(ie,1,jb) = SUM(z_quad_vector(ie,:,1))
        p_quad_vector_sum(ie,2,jb) = SUM(z_quad_vector(ie,:,2))
        p_quad_vector_sum(ie,3,jb) = SUM(z_quad_vector(ie,:,3))
        p_quad_vector_sum(ie,4,jb) = SUM(z_quad_vector(ie,:,4))
        p_quad_vector_sum(ie,5,jb) = SUM(z_quad_vector(ie,:,5))
        p_quad_vector_sum(ie,6,jb) = SUM(z_quad_vector(ie,:,6))

      ENDDO ! ie: loop over index list
!$ACC END PARALLEL


!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG VECTOR
!CDIR NODEP,VOVERTAKE,VOB
      DO ie = 1, falist%len(jb)

        je = falist%eidx(ie,jb)
        jk = falist%elev(ie,jb)

        ! Add contribution to total area of departure region
        p_dreg_area(je,jk,jb) = p_dreg_area(je,jk,jb) + SUM(wgt_t_detjac(ie,1:4))

      ENDDO ! ie: loop over index list
!$ACC END PARALLEL

    ENDDO  ! loop over blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( p_dreg_area, p_quad_vector_sum ), IF ( acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA

  END SUBROUTINE prep_gauss_quadrature_q_list


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
    &                                     p_quad_vector_sum, p_dreg_area,   &
    &                                     opt_rlstart, opt_rlend, opt_slev, &
    &                                     opt_elev                          )

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(IN) :: & !< patch on which computation is
      &  p_patch                           !< performed

    REAL(vp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:,:) !< in 2D cartesian coordinates
                                    !< dim: (nproma,4,2,nlev,nblks_e)

    REAL(vp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:) !< dim: (nproma,8,nlev,nblks_e)

    REAL(vp), INTENT(OUT) :: &      !< area of departure region  [m**2]
      &  p_dreg_area(:,:,:)         !< dim: (nproma,nlev,nblks_e)

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
      &  z_gauss_pts(nproma,4,2)    !< in physical space

    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac(nproma,4)     !< each gaussian quadrature point.

    REAL(wp) ::                &    !< quadrature vector for single integration point
      &  z_quad_vector(nproma,4,8)

    REAL(wp) :: z_x(nproma,4), z_y(nproma,4) !< storage for local coordinates
    REAL(wp) :: z_wgt(4), z_eta(4,4)         !< for precomputation of coefficients
    REAL(wp) :: z_area                       !< auxiliary for dreg area

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: slev, elev          !< vertical start and end level

  !-----------------------------------------------------------------------

!$ACC DATA PCOPYIN( p_coords_dreg_v ), PCOPYOUT( p_quad_vector_sum, p_dreg_area ), &
!$ACC      CREATE( z_gauss_pts, wgt_t_detjac, z_quad_vector, z_x, z_y, z_wgt, z_eta ), &
!$ACC      IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_coords_dreg_v, p_dreg_area ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

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
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

    z_wgt(1) = 0.0625_wp * wgt_zeta(1) *  wgt_eta(1)
    z_wgt(2) = 0.0625_wp * wgt_zeta(1) *  wgt_eta(2)
    z_wgt(3) = 0.0625_wp * wgt_zeta(2) *  wgt_eta(1)
    z_wgt(4) = 0.0625_wp * wgt_zeta(2) *  wgt_eta(2)

    z_eta(1,1:4) = 1._wp - eta(1:4)
    z_eta(2,1:4) = 1._wp + eta(1:4)
    z_eta(3,1:4) = 1._wp - zeta(1:4)
    z_eta(4,1:4) = 1._wp + zeta(1:4)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,jg,i_startidx,i_endidx,z_gauss_pts,wgt_t_detjac,&
!$OMP z_quad_vector,z_x,z_y,z_area) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG PRIVATE( z_x, z_y, wgt_t_detjac, z_gauss_pts, z_quad_vector )
      DO jk = slev, elev

        !$ACC LOOP VECTOR
        DO je = i_startidx, i_endidx

          z_x(je,1:4) = p_coords_dreg_v(je,1:4,1,jk,jb)
          z_y(je,1:4) = p_coords_dreg_v(je,1:4,2,jk,jb)

          ! get Jacobian determinant for each quadrature point and multiply with
          ! corresponding weights
          ! Note: dbl_eps is added, in order to have a meaningful 'edge value' 
          ! (better: area-average) even when the integration-area tends to zero.
          wgt_t_detjac(je,1:4) = dbl_eps + z_wgt(1:4) * (                                 &
            &   (z_eta(1,1:4)*(z_x(je,2)-z_x(je,1)) + z_eta(2,1:4)*(z_x(je,3)-z_x(je,4))) &
            & * (z_eta(3,1:4)*(z_y(je,4)-z_y(je,1)) - z_eta(4,1:4)*(z_y(je,2)-z_y(je,3))) &
            & - (z_eta(1,1:4)*(z_y(je,2)-z_y(je,1)) + z_eta(2,1:4)*(z_y(je,3)-z_y(je,4))) &
            & * (z_eta(3,1:4)*(z_x(je,4)-z_x(je,1)) - z_eta(4,1:4)*(z_x(je,2)-z_x(je,3))) )


          ! get coordinates of the quadrature points in physical space (mapping)
          z_gauss_pts(je,1,1) = DOT_PRODUCT(shape_func(1:4,1),z_x(je,1:4))
          z_gauss_pts(je,1,2) = DOT_PRODUCT(shape_func(1:4,1),z_y(je,1:4))
          z_gauss_pts(je,2,1) = DOT_PRODUCT(shape_func(1:4,2),z_x(je,1:4))
          z_gauss_pts(je,2,2) = DOT_PRODUCT(shape_func(1:4,2),z_y(je,1:4))
          z_gauss_pts(je,3,1) = DOT_PRODUCT(shape_func(1:4,3),z_x(je,1:4))
          z_gauss_pts(je,3,2) = DOT_PRODUCT(shape_func(1:4,3),z_y(je,1:4))
          z_gauss_pts(je,4,1) = DOT_PRODUCT(shape_func(1:4,4),z_x(je,1:4))
          z_gauss_pts(je,4,2) = DOT_PRODUCT(shape_func(1:4,4),z_y(je,1:4))


          ! Get quadrature vector for each integration point and multiply by
          ! corresponding wgt_t_detjac
          DO jg=1, 4
            z_quad_vector(je,jg,1) = wgt_t_detjac(je,jg)
            z_quad_vector(je,jg,2) = wgt_t_detjac(je,jg) * z_gauss_pts(je,jg,1)
            z_quad_vector(je,jg,3) = wgt_t_detjac(je,jg) * z_gauss_pts(je,jg,2)
            z_quad_vector(je,jg,4) = wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,1)**2)
            z_quad_vector(je,jg,5) = wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,2)**2)
            z_quad_vector(je,jg,6) = wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,1) * z_gauss_pts(je,jg,2))
            z_quad_vector(je,jg,7) = wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,1)**3)
            z_quad_vector(je,jg,8) = wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,2)**3)
          ENDDO


          ! Sum quadrature vectors over all integration points
          p_quad_vector_sum(je,1,jk,jb) = SUM(z_quad_vector(je,:,1))
          p_quad_vector_sum(je,2,jk,jb) = SUM(z_quad_vector(je,:,2))
          p_quad_vector_sum(je,3,jk,jb) = SUM(z_quad_vector(je,:,3))
          p_quad_vector_sum(je,4,jk,jb) = SUM(z_quad_vector(je,:,4))
          p_quad_vector_sum(je,5,jk,jb) = SUM(z_quad_vector(je,:,5))
          p_quad_vector_sum(je,6,jk,jb) = SUM(z_quad_vector(je,:,6))
          p_quad_vector_sum(je,7,jk,jb) = SUM(z_quad_vector(je,:,7))
          p_quad_vector_sum(je,8,jk,jb) = SUM(z_quad_vector(je,:,8))


          ! area of departure region
          z_area = SUM(wgt_t_detjac(je,1:4))
          p_dreg_area(je,jk,jb) = SIGN(MAX(eps,ABS(z_area)),z_area)

        ENDDO ! loop over edges

      ENDDO  ! loop over levels
!$ACC END PARALLEL

    ENDDO  ! loop over blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( p_dreg_area, p_quad_vector_sum ), IF (acc_validate .AND.  i_am_accel_node .AND. acc_on )
!$ACC END DATA

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
    &                                 p_quad_vector_sum, p_dreg_area,   &
    &                                 opt_rlstart, opt_rlend, opt_slev, &
    &                                 opt_elev                          )

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                           !< performed

    REAL(vp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:,:) !< in 2D cartesian coordinates
                                    !< dim: (nproma,4,2,nlev,nblks_e)

    REAL(vp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:) !< dim: (nproma,10,nlev,nblks_e)

    REAL(vp), INTENT(OUT) :: &      !< area of departure region  [m**2]
      &  p_dreg_area(:,:,:)         !< dim: (nproma,nlev,nblks_e)

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
      &  z_gauss_pts(nproma,4,2)    !< in physical space

    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac(nproma,4)      !< each gaussian quadrature point.

    REAL(wp) ::                &    !< quadrature vector for single integration point
      &  z_quad_vector(nproma,4,10)

    REAL(wp) :: z_x(nproma,4), z_y(nproma,4) !< storage for local coordinates
    REAL(wp) :: z_wgt(4), z_eta(4,4)         !< for precomputation of coefficients
    REAL(wp) :: z_area                       !< auxiliary for dreg area

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: slev, elev          !< vertical start and end level

  !-----------------------------------------------------------------------

!$ACC DATA PCOPYIN( p_coords_dreg_v ), PCOPYOUT( p_quad_vector_sum, p_dreg_area ), &
!$ACC      CREATE( z_gauss_pts, wgt_t_detjac, z_quad_vector, z_x, z_y, z_wgt, z_eta ), &
!$ACC      IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_coords_dreg_v ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

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
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

    z_wgt(1) = 0.0625_wp * wgt_zeta(1) *  wgt_eta(1)
    z_wgt(2) = 0.0625_wp * wgt_zeta(1) *  wgt_eta(2)
    z_wgt(3) = 0.0625_wp * wgt_zeta(2) *  wgt_eta(1)
    z_wgt(4) = 0.0625_wp * wgt_zeta(2) *  wgt_eta(2)

    z_eta(1,1:4) = 1._wp - eta(1:4)
    z_eta(2,1:4) = 1._wp + eta(1:4)
    z_eta(3,1:4) = 1._wp - zeta(1:4)
    z_eta(4,1:4) = 1._wp + zeta(1:4)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,jg,i_startidx,i_endidx,z_gauss_pts,wgt_t_detjac,&
!$OMP z_quad_vector,z_x,z_y,z_area) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG PRIVATE( z_x, z_y, wgt_t_detjac, z_gauss_pts, z_quad_vector )
      DO jk = slev, elev

        !$ACC LOOP VECTOR
        DO je = i_startidx, i_endidx

          z_x(je,1:4) = p_coords_dreg_v(je,1:4,1,jk,jb)
          z_y(je,1:4) = p_coords_dreg_v(je,1:4,2,jk,jb)

          ! get Jacobian determinant for each quadrature point and multiply with
          ! corresponding weights
          ! Note: dbl_eps is added, in order to have a meaningful 'edge value' 
          ! (better: area-average) even when the integration-area tends to zero.
          wgt_t_detjac(je,1:4) = dbl_eps + z_wgt(1:4) * (                                 &
            &   (z_eta(1,1:4)*(z_x(je,2)-z_x(je,1)) + z_eta(2,1:4)*(z_x(je,3)-z_x(je,4))) &
            & * (z_eta(3,1:4)*(z_y(je,4)-z_y(je,1)) - z_eta(4,1:4)*(z_y(je,2)-z_y(je,3))) &
            & - (z_eta(1,1:4)*(z_y(je,2)-z_y(je,1)) + z_eta(2,1:4)*(z_y(je,3)-z_y(je,4))) &
            & * (z_eta(3,1:4)*(z_x(je,4)-z_x(je,1)) - z_eta(4,1:4)*(z_x(je,2)-z_x(je,3))) )


          ! get coordinates of the quadrature points in physical space (mapping)
          z_gauss_pts(je,1,1) = DOT_PRODUCT(shape_func(1:4,1),z_x(je,1:4))
          z_gauss_pts(je,1,2) = DOT_PRODUCT(shape_func(1:4,1),z_y(je,1:4))
          z_gauss_pts(je,2,1) = DOT_PRODUCT(shape_func(1:4,2),z_x(je,1:4))
          z_gauss_pts(je,2,2) = DOT_PRODUCT(shape_func(1:4,2),z_y(je,1:4))
          z_gauss_pts(je,3,1) = DOT_PRODUCT(shape_func(1:4,3),z_x(je,1:4))
          z_gauss_pts(je,3,2) = DOT_PRODUCT(shape_func(1:4,3),z_y(je,1:4))
          z_gauss_pts(je,4,1) = DOT_PRODUCT(shape_func(1:4,4),z_x(je,1:4))
          z_gauss_pts(je,4,2) = DOT_PRODUCT(shape_func(1:4,4),z_y(je,1:4))


          ! Get quadrature vector for each integration point and multiply by
          ! corresponding wgt_t_detjac
          DO jg=1, 4
            z_quad_vector(je,jg,1) = wgt_t_detjac(je,jg)
            z_quad_vector(je,jg,2) = wgt_t_detjac(je,jg) * z_gauss_pts(je,jg,1)
            z_quad_vector(je,jg,3) = wgt_t_detjac(je,jg) * z_gauss_pts(je,jg,2)
            z_quad_vector(je,jg,4) = wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,1)**2)
            z_quad_vector(je,jg,5) = wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,2)**2)
            z_quad_vector(je,jg,6) = wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,1) * z_gauss_pts(je,jg,2))
            z_quad_vector(je,jg,7) = wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,1)**3)
            z_quad_vector(je,jg,8) = wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,2)**3)
            z_quad_vector(je,jg,9) = wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,1)**2 * z_gauss_pts(je,jg,2))
            z_quad_vector(je,jg,10)= wgt_t_detjac(je,jg) * (z_gauss_pts(je,jg,1) * z_gauss_pts(je,jg,2)**2)
          ENDDO


          ! Sum quadrature vectors over all integration points
          p_quad_vector_sum(je, 1,jk,jb) = SUM(z_quad_vector(je,:,1))
          p_quad_vector_sum(je, 2,jk,jb) = SUM(z_quad_vector(je,:,2))
          p_quad_vector_sum(je, 3,jk,jb) = SUM(z_quad_vector(je,:,3))
          p_quad_vector_sum(je, 4,jk,jb) = SUM(z_quad_vector(je,:,4))
          p_quad_vector_sum(je, 5,jk,jb) = SUM(z_quad_vector(je,:,5))
          p_quad_vector_sum(je, 6,jk,jb) = SUM(z_quad_vector(je,:,6))
          p_quad_vector_sum(je, 7,jk,jb) = SUM(z_quad_vector(je,:,7))
          p_quad_vector_sum(je, 8,jk,jb) = SUM(z_quad_vector(je,:,8))
          p_quad_vector_sum(je, 9,jk,jb) = SUM(z_quad_vector(je,:,9))
          p_quad_vector_sum(je,10,jk,jb) = SUM(z_quad_vector(je,:,10))


          ! area of departure region
          z_area = SUM(wgt_t_detjac(je,1:4))
          p_dreg_area(je,jk,jb) = SIGN(MAX(eps,ABS(z_area)),z_area)

!!$IF (p_dreg_area(je,jk,jb) < 0._wp) THEN
!!$  WRITE(0,*) "ATTENTION: negative areas at je,jk,jb= ", je, jk, jb, p_dreg_area(je,jk,jb)
!!$  WRITE(0,*) "system orientation: ", p_patch%edges%tangent_orientation(je,jb)
!!$  ELSE IF ((p_dreg_area(je,jk,jb) >= 0._wp)) THEN
!!$  WRITE(0,*) "OK for system orientation= ", je, jk, jb, p_patch%edges%tangent_orientation(je,jb)
!!$ENDIF
        ENDDO ! loop over edges

      ENDDO  ! loop over levels
!$ACC END PARALLEL

    ENDDO  ! loop over blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( p_dreg_area, p_quad_vector_sum ), IF ( acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA

  END SUBROUTINE prep_gauss_quadrature_c



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
  !! Index-list based version. Otherwise identical to prep_gauss_quadrature_c
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2013-10-28)
  !!
  !!
  SUBROUTINE prep_gauss_quadrature_c_list( p_patch, p_coords_dreg_v, falist, &
    &                                 p_quad_vector_sum, p_dreg_area,   &
    &                                 opt_rlstart, opt_rlend, opt_slev, &
    &                                 opt_elev                          )

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                           !< performed

    REAL(vp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:)   !< in 2D cartesian coordinates
                                    !< dim: (npoints,4,2,nblks_e)

    TYPE(t_list2D), INTENT(IN) :: & !< index list with points for which the standard 
      &  falist                     !< Miura-type treatment of flux areas is 
                                    !< insufficient

    REAL(vp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:)   !< dim: (npoints,10,nblks_e)

    REAL(vp), INTENT(INOUT) :: &    !< total area of departure region  [m**2]
      &  p_dreg_area(:,:,:)         !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

   ! local variables
    REAL(wp) ::                        &    !< coordinates of gaussian quadrature points
      &  z_gauss_pts(falist%npoints,4,2)    !< in physical space

    REAL(wp) ::                       &     !< weights times determinant of Jacobian for
      &  wgt_t_detjac(falist%npoints,4)     !< each gaussian quadrature point.

    REAL(wp) ::                           & !< quadrature vector for single integration point
      &  z_quad_vector(falist%npoints,4,10)

    REAL(wp) :: z_x(falist%npoints,4), z_y(falist%npoints,4) !< storage for local coordinates
    REAL(wp) :: z_wgt(4), z_eta(4,4)         !< for precomputation of coefficients

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: ie                  !< index list loop counter
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: slev, elev          !< vertical start and end level

  !-----------------------------------------------------------------------

!$ACC DATA PCOPYIN( p_coords_dreg_v, falist ), PCOPY( p_dreg_area ), PCOPYOUT( p_quad_vector_sum ), &
!$ACC      CREATE( z_gauss_pts, wgt_t_detjac, z_quad_vector, z_x, z_y, z_wgt, z_eta ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_coords_dreg_v, p_dreg_area ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )
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
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

    z_wgt(1) = 0.0625_wp * wgt_zeta(1) *  wgt_eta(1)
    z_wgt(2) = 0.0625_wp * wgt_zeta(1) *  wgt_eta(2)
    z_wgt(3) = 0.0625_wp * wgt_zeta(2) *  wgt_eta(1)
    z_wgt(4) = 0.0625_wp * wgt_zeta(2) *  wgt_eta(2)

    z_eta(1,1:4) = 1._wp - eta(1:4)
    z_eta(2,1:4) = 1._wp + eta(1:4)
    z_eta(3,1:4) = 1._wp - zeta(1:4)
    z_eta(4,1:4) = 1._wp + zeta(1:4)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,ie,jg,z_gauss_pts,wgt_t_detjac,&
!$OMP z_quad_vector,z_x,z_y) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG VECTOR
      DO ie = 1, falist%len(jb)

        z_x(ie,1:4) = p_coords_dreg_v(ie,1:4,1,jb)
        z_y(ie,1:4) = p_coords_dreg_v(ie,1:4,2,jb)

        ! get Jacobian determinant for each quadrature point and multiply with
        ! corresponding weights
        ! Note: dbl_eps is added, in order to have a meaningful 'edge value' 
        ! (better: area-average) even when the integration-area tends to zero.
        wgt_t_detjac(ie,1:4) = dbl_eps + z_wgt(1:4) * (                                 &
          &   (z_eta(1,1:4)*(z_x(ie,2)-z_x(ie,1)) + z_eta(2,1:4)*(z_x(ie,3)-z_x(ie,4))) &
          & * (z_eta(3,1:4)*(z_y(ie,4)-z_y(ie,1)) - z_eta(4,1:4)*(z_y(ie,2)-z_y(ie,3))) &
          & - (z_eta(1,1:4)*(z_y(ie,2)-z_y(ie,1)) + z_eta(2,1:4)*(z_y(ie,3)-z_y(ie,4))) &
          & * (z_eta(3,1:4)*(z_x(ie,4)-z_x(ie,1)) - z_eta(4,1:4)*(z_x(ie,2)-z_x(ie,3))) )


        ! get coordinates of the quadrature points in physical space (mapping)
!WS: TODO: make sure DOT_PRODUCT works in this OpenACC context
        z_gauss_pts(ie,1,1) = DOT_PRODUCT(shape_func(1:4,1),z_x(ie,1:4))
        z_gauss_pts(ie,1,2) = DOT_PRODUCT(shape_func(1:4,1),z_y(ie,1:4))
        z_gauss_pts(ie,2,1) = DOT_PRODUCT(shape_func(1:4,2),z_x(ie,1:4))
        z_gauss_pts(ie,2,2) = DOT_PRODUCT(shape_func(1:4,2),z_y(ie,1:4))
        z_gauss_pts(ie,3,1) = DOT_PRODUCT(shape_func(1:4,3),z_x(ie,1:4))
        z_gauss_pts(ie,3,2) = DOT_PRODUCT(shape_func(1:4,3),z_y(ie,1:4))
        z_gauss_pts(ie,4,1) = DOT_PRODUCT(shape_func(1:4,4),z_x(ie,1:4))
        z_gauss_pts(ie,4,2) = DOT_PRODUCT(shape_func(1:4,4),z_y(ie,1:4))


        ! Get quadrature vector for each integration point and multiply by
        ! corresponding wgt_t_detjac
        DO jg=1, 4
          z_quad_vector(ie,jg,1) = wgt_t_detjac(ie,jg)
          z_quad_vector(ie,jg,2) = wgt_t_detjac(ie,jg) * z_gauss_pts(ie,jg,1)
          z_quad_vector(ie,jg,3) = wgt_t_detjac(ie,jg) * z_gauss_pts(ie,jg,2)
          z_quad_vector(ie,jg,4) = wgt_t_detjac(ie,jg) * (z_gauss_pts(ie,jg,1)**2)
          z_quad_vector(ie,jg,5) = wgt_t_detjac(ie,jg) * (z_gauss_pts(ie,jg,2)**2)
          z_quad_vector(ie,jg,6) = wgt_t_detjac(ie,jg) * (z_gauss_pts(ie,jg,1) * z_gauss_pts(ie,jg,2))
          z_quad_vector(ie,jg,7) = wgt_t_detjac(ie,jg) * (z_gauss_pts(ie,jg,1)**3)
          z_quad_vector(ie,jg,8) = wgt_t_detjac(ie,jg) * (z_gauss_pts(ie,jg,2)**3)
          z_quad_vector(ie,jg,9) = wgt_t_detjac(ie,jg) * (z_gauss_pts(ie,jg,1)**2 * z_gauss_pts(ie,jg,2))
          z_quad_vector(ie,jg,10)= wgt_t_detjac(ie,jg) * (z_gauss_pts(ie,jg,1) * z_gauss_pts(ie,jg,2)**2)
        ENDDO


        ! Sum quadrature vectors over all integration points
        p_quad_vector_sum(ie, 1,jb) = SUM(z_quad_vector(ie,:,1))
        p_quad_vector_sum(ie, 2,jb) = SUM(z_quad_vector(ie,:,2))
        p_quad_vector_sum(ie, 3,jb) = SUM(z_quad_vector(ie,:,3))
        p_quad_vector_sum(ie, 4,jb) = SUM(z_quad_vector(ie,:,4))
        p_quad_vector_sum(ie, 5,jb) = SUM(z_quad_vector(ie,:,5))
        p_quad_vector_sum(ie, 6,jb) = SUM(z_quad_vector(ie,:,6))
        p_quad_vector_sum(ie, 7,jb) = SUM(z_quad_vector(ie,:,7))
        p_quad_vector_sum(ie, 8,jb) = SUM(z_quad_vector(ie,:,8))
        p_quad_vector_sum(ie, 9,jb) = SUM(z_quad_vector(ie,:,9))
        p_quad_vector_sum(ie,10,jb) = SUM(z_quad_vector(ie,:,10))

      ENDDO ! ie: loop over index list
!$ACC END PARALLEL

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG VECTOR
!CDIR NODEP,VOVERTAKE,VOB
      DO ie = 1, falist%len(jb)

        je = falist%eidx(ie,jb)
        jk = falist%elev(ie,jb)

        ! Add contribution to total area of departure region
        p_dreg_area(je,jk,jb) = p_dreg_area(je,jk,jb) + SUM(wgt_t_detjac(ie,1:4))

      ENDDO ! ie: loop over index list
!$ACC END PARALLEL

    ENDDO  ! loop over blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( p_dreg_area, p_quad_vector_sum ), IF ( acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA

  END SUBROUTINE prep_gauss_quadrature_c_list


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
!$ACC ROUTINE SEQ
    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: x(1:4), y(1:4)  !< coordinates of vertices in x-y-system
    REAL(wp), INTENT(IN) :: zeta, eta       !< integration point in \zeta,\eta-system

    ! RETURN VALUE:
    REAL(wp) :: det_jac

    REAL(wp), DIMENSION(2,2) :: jacob

  !-----------------------------------------------------------------------

    jacob(1,1) = (1._wp - eta) * ( x(2) - x(1))   &
      &        + (1._wp + eta) * ( x(3) - x(4))
    jacob(1,2) = (1._wp - eta) * ( y(2) - y(1))   &
      &        + (1._wp + eta) * ( y(3) - y(4))
    jacob(2,1) = (1._wp - zeta)* ( x(4) - x(1))   &
      &        - (1._wp + zeta)* ( x(2) - x(3))
    jacob(2,2) = (1._wp - zeta)* ( y(4) - y(1))   &
      &        - (1._wp + zeta)* ( y(2) - y(3))

    det_jac = 0.0625_wp * (jacob(1,1)*jacob(2,2) - jacob(1,2)*jacob(2,1))

  END FUNCTION jac


END MODULE mo_advection_quadrature

