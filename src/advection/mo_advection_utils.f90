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
MODULE mo_advection_utils

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_impl_constants,      ONLY: min_rlcell_int
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: laxfr_upflux
  PUBLIC :: laxfr_upflux_v
  PUBLIC :: ptr_delp_mc_now
  PUBLIC :: ptr_delp_mc_new

  
  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  ! In order to avoid circular dependencies these two pointers
  ! have been moved from mo_advection_stepping to this module.
  REAL(wp), POINTER ::  &
    &  ptr_delp_mc_now(:,:,:) => NULL() !< pointer to old layer thickness
                                        !< at cell center
  REAL(wp), POINTER ::  &
    &  ptr_delp_mc_new(:,:,:) => NULL() !< pointer to new layer thickness
                                        !< at cell center

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
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
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
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE tupdate_tracer


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


END MODULE mo_advection_utils

