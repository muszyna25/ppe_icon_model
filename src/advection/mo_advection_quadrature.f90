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
MODULE mo_advection_quadrature

  USE mo_advection_config,    ONLY: shape_func, zeta, eta, wgt_zeta, wgt_eta
  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_loopindices,         ONLY: get_indices_e
  USE mo_impl_constants,      ONLY: min_rledge_int
  USE mo_math_constants,      ONLY: dbl_eps


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prep_gauss_quadrature_q
  PUBLIC :: prep_gauss_quadrature_cpoor
  PUBLIC :: prep_gauss_quadrature_c

  
  CHARACTER(len=*), PARAMETER :: version = '$Id$'



CONTAINS


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


END MODULE mo_advection_quadrature

