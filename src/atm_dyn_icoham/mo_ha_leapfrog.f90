!>
!!  Contains the subroutine for the explicit leap frog  time integration scheme.
!!
!!
!! @par Revision History
!!  <i>SUBROUTINE asselin</i> by L.Bonaventura, Polimi (2006).
!!  <i>SUBROUTINE step_leapfrog_expl</i> by Hui Wan, MPI-M (2006-09-10)
!!  Modification by Almut Gassmann, MPI-M (2008-09-19)
!!  - Code restructuring, remove loop over grid levels etc.
!!  Modification by Marco Giorgetta, MPI-M (2009-03-27)
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ha_leapfrog

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE m_dyn,                  ONLY: dyn_theta
  USE mo_ha_dynamics,         ONLY: dyn_temp
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  step_leapfrog_expl, asselin, leapfrog_update_prog

CONTAINS

  !>
  !!
  !! @par Revision History
  !!  Original version  by Hui Wan, MPI-M (2006-09-10).
  !!  Modification by Hui Wan, MPI-M (2007-07-30)
  !!  - pointers to *grid* removed
  !!  - second dimension of pointers (i.e. index of patches on the same
  !!    grid level) removed.
  !!  Modification by Hui Wan, MPI-M (2007-09-16)
  !!  - horizontal diffusion moved to a separate module and called from
  !!    the main program.
  !!  Code restructuring by Almut Gassmann, MPI-M (2008-09-19)
  !!  Modification by Hui Wan, MPI-M (2009-08-06):
  !!  - Packed the updating of prognostic variables into the subroutine
  !!    leapfrog_update_prog.
  !!
  SUBROUTINE step_leapfrog_expl( pdtime, dtime_bdy,      & ! input
                                 ltheta_dyn,             & ! input
                                 curr_patch, p_int_state,& ! input
                                 p_old,                  & ! input
                                 p_ext_data,             & ! input
                                 p_now,                  & ! in and out
                                 p_diag,                 & ! in and out
                                 p_new,                  & ! in and out
                                 p_tend_dyn             )  ! in and out

    REAL(wp),INTENT(IN) :: pdtime      ! time step in seconds
    REAL(wp),INTENT(IN) :: dtime_bdy   ! time step for boundary tendencies
    LOGICAL, INTENT(IN) :: ltheta_dyn

    TYPE(t_patch), TARGET, INTENT(INOUT) :: curr_patch
    TYPE(t_int_state),TARGET,INTENT(IN) :: p_int_state
    TYPE(t_hydro_atm_prog),  INTENT(IN) :: p_old
    TYPE(t_external_data),   INTENT(INOUT) :: p_ext_data !< external data

    TYPE(t_hydro_atm_diag),INTENT(INOUT) :: p_diag
    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: p_now
    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: p_new
    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: p_tend_dyn

    REAL(wp) :: zdt2    ! double time step

   !-----------------------------------------------------------------------

   zdt2 = 2._wp*pdtime

   ! calculate tendency of the prognostic variables

    IF (ltheta_dyn) THEN
       CALL dyn_theta( curr_patch, p_int_state, p_ext_data, p_now, &! in
                       p_diag, p_tend_dyn )                         ! out
    ELSE
       CALL dyn_temp(  curr_patch, p_int_state, p_now, p_ext_data, &! in
                       p_diag, p_tend_dyn )                         ! out
    END IF

    ! update prognostic variables

    CALL leapfrog_update_prog( p_new, p_now, p_old,      &! inout,in,in
                               p_tend_dyn,               &! in
                               zdt2, dtime_bdy, .FALSE., &! in. Do not touch tracers
                               ltheta_dyn, curr_patch    )! in

  END SUBROUTINE step_leapfrog_expl

  !-------------------------------------------------------------------------
  !>
  !!  Update the prognostic variables using the given tendencies and  time steps.
  !!
  !! The new values are store in p_new.
  !!  Nesting boundaries and domain interior are treated separately.
  !!
  !! @par Revision History
  !!  Separated from step_leapfrog_expl by Hui Wan, MPI-M (2009-08-06).
  !!  Switch ltracer added by Hui Wan, MPI-M (2010-07-27)
  !!
  SUBROUTINE leapfrog_update_prog( p_new, p_now, p_old,  & ! inout,in,in
                                   p_tend,               & ! in
                                   pdt2, dtime_bdy,      & ! in
                                   ltracer,              & ! in
                                   ltheta_dyn,           & ! in
                                   curr_patch            ) ! in

   REAL(wp), INTENT(IN) :: pdt2       !< time step in seconds for the interior
   REAL(wp), INTENT(IN) :: dtime_bdy  !< time step for boundary tendencies
   LOGICAL,  INTENT(IN) :: ltracer    !< if .TRUE., update tracer fields
   LOGICAL,  INTENT(IN) :: ltheta_dyn !< if .TRUE., update tracer fields

   TYPE(t_patch),        INTENT(IN) :: curr_patch  !< domain info.
   TYPE(t_hydro_atm_prog),INTENT(IN) :: p_old
   TYPE(t_hydro_atm_prog),INTENT(IN) :: p_now
   TYPE(t_hydro_atm_prog),INTENT(IN) :: p_tend

   TYPE(t_hydro_atm_prog),INTENT(INOUT) :: p_new  !< intent = inout to keep the
                                                  !< tracer component intact when
                                                  !< ltracer = .FALSE.
   INTEGER :: jb, jbs, jbe, is, ie

!$OMP PARALLEL PRIVATE(jbs,jbe)
    !-----------------------------------
    !step forward - model interior
    !-----------------------------------
    jbs = curr_patch%cells%start_blk(grf_bdywidth_c+1,1)
    jbe = curr_patch%nblks_c
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs, jbe
      CALL get_indices_c(curr_patch, jb,jbs,jbe, is,ie, grf_bdywidth_c+1)

      p_new%pres_sfc(is:ie,jb) = p_old%pres_sfc(is:ie,jb) + pdt2*p_tend%pres_sfc(is:ie,jb)

      IF (ltheta_dyn) THEN
        p_new%theta(is:ie,:,jb) = p_old%theta(is:ie,:,jb) + pdt2*p_tend%temp(is:ie,:,jb)
      ELSE
        p_new% temp(is:ie,:,jb) = p_old% temp(is:ie,:,jb) + pdt2*p_tend%temp(is:ie,:,jb)
      ENDIF

      IF (ltracer) THEN
        p_new%tracer(is:ie,:,jb,:) =    p_old%tracer(is:ie,:,jb,:)      &
                                   & + p_tend%tracer(is:ie,:,jb,:)*pdt2
      ENDIF
    ENDDO
!$OMP END DO

    jbs = curr_patch%edges%start_blk(grf_bdywidth_e+1,1)
    jbe = curr_patch%nblks_e
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs,jbe
       CALL get_indices_e(curr_patch, jb,jbs,jbe, is,ie, grf_bdywidth_e+1)
       p_new%vn(is:ie,:,jb) = p_old%vn(is:ie,:,jb) + pdt2*p_tend%vn(is:ie,:,jb)
    ENDDO ! block loop
!$OMP END DO

    !-----------------------------------
    ! Update of nest boundaries
    !-----------------------------------
    ! Note: the boundary tendencies are nnow -> nnew
    jbs = curr_patch%cells%start_blk(1,1)
    jbe = curr_patch%cells%end_blk(grf_bdywidth_c,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs,jbe
       CALL get_indices_c(curr_patch, jb,jbs,jbe, is,ie, 1, grf_bdywidth_c)

       p_new%pres_sfc(is:ie,jb) =    p_now%pres_sfc(is:ie,jb)           &
                                & + p_tend%pres_sfc(is:ie,jb)*dtime_bdy

       IF (ltheta_dyn) THEN
          p_new%theta(is:ie,:,jb) =    p_now%theta(is:ie,:,jb)          &
                                  & + p_tend% temp(is:ie,:,jb)*dtime_bdy
       ELSE
          p_new% temp(is:ie,:,jb) =    p_now% temp(is:ie,:,jb)          &
                                  & + p_tend% temp(is:ie,:,jb)*dtime_bdy
       ENDIF

       IF (ltracer) THEN
          p_new%tracer(is:ie,:,jb,:) =    p_now%tracer(is:ie,:,jb,:)           &
                                     & + p_tend%tracer(is:ie,:,jb,:)*dtime_bdy
       ENDIF
    ENDDO
!$OMP END DO

    jbs = curr_patch%edges%start_blk(1,1)
    jbe = curr_patch%edges%end_blk(grf_bdywidth_e,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs,jbe
       CALL get_indices_e(curr_patch, jb,jbs,jbe, is,ie, 1, grf_bdywidth_e)
       p_new%vn(is:ie,:,jb) = p_now%vn(is:ie,:,jb) + dtime_bdy*p_tend%vn(is:ie,:,jb)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE leapfrog_update_prog
  !-------------
  !>
  !!
  !! @par Revision History
  !! Developed  by L.Bonaventura, Polimi (2006).
  !! Code restructuring by Almut Gassmann, MPI-M, (2008-09-19)
  !!
  SUBROUTINE asselin( asselin_coeff, ltheta_dyn, &! in
                      p_prog_old, p_prog_new,    &! in
                      p_prog_now                 )! inout

    REAL(wp),INTENT(IN) :: asselin_coeff
    LOGICAL, INTENT(IN) :: ltheta_dyn

    TYPE(t_hydro_atm_prog), INTENT(INOUT) :: p_prog_now
    TYPE(t_hydro_atm_prog), INTENT(IN)    :: p_prog_new
    TYPE(t_hydro_atm_prog), INTENT(IN)    :: p_prog_old


!$OMP PARALLEL
    IF (ltheta_dyn) THEN
!$OMP WORKSHARE
      p_prog_now%theta = p_prog_now%theta + asselin_coeff* &
                        (p_prog_new%theta - 2._wp*p_prog_now%theta + p_prog_old%theta)
!$OMP END WORKSHARE
    ENDIF

!$OMP WORKSHARE
    p_prog_now%temp = p_prog_now%temp + asselin_coeff* &
                     (p_prog_new%temp -2._wp*p_prog_now%temp +p_prog_old%temp)

    p_prog_now%vn = p_prog_now%vn + asselin_coeff* &
                   (p_prog_new%vn -2._wp*p_prog_now%vn +p_prog_old%vn)

    p_prog_now%pres_sfc = p_prog_now%pres_sfc + asselin_coeff*          &
                         (p_prog_new%pres_sfc-2._wp*p_prog_now%pres_sfc &
                         +p_prog_old%pres_sfc)
!$OMP END WORKSHARE
!$OMP END PARALLEL

  END SUBROUTINE asselin

END MODULE mo_ha_leapfrog


