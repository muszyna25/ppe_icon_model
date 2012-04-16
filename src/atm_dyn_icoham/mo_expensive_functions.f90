!>
!! Conversion between potential temperature and temperature.
!!
!! @par Revision History
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
MODULE mo_expensive_functions
!RJ: For some strange reasons the Intel compiler produces wrong code when
!optimizing convert_theta2t_lin, therefore we switch off optimization here:
!DEC$ NOOPTIMIZE
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_run_config,         ONLY: nlevp1, nlev
  USE mo_physical_constants, ONLY: rd, cpd, p0ref
  USE mo_vertical_coord_table,ONLY: nplvp1, delpr, nplev
  USE mo_model_domain,       ONLY: t_patch
  USE mo_icoham_dyn_types,   ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_eta_coord_diag,     ONLY: half_level_pressure,  &
                                   full_level_pressure
  USE mo_timer,              ONLY: timer_start, timer_stop,ltimer,  &
    & timer_con_l_theta2t, timer_con_l_t2theta, timer_con_theta2t, timer_con_t2theta


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC ::  convert_t2theta,     convert_theta2t,    &
             convert_t2theta_lin, convert_theta2t_lin

CONTAINS
!-------------------------------------------------------------------------
!
!
  !>
  !!               Converts temperature into potential temperature times delta p.
  !!
  !!
  !! @par Revision History
  !!   Created by Guenther Zaengl 2009-01-13
  !!
  SUBROUTINE convert_t2theta( pt_patch, pt_prog, pt_diag)

  IMPLICIT NONE

  TYPE(t_patch),           INTENT(IN)    :: pt_patch !< grid/patch info.
  TYPE(t_hydro_atm_prog),INTENT(INOUT) :: pt_prog  !< prognostic variables
  TYPE(t_hydro_atm_diag),INTENT(INOUT) :: pt_diag  !< diagnostic variables

  ! Local array bounds:

  INTEGER :: nlen          !< arb. field length
  INTEGER :: nblks_c       !< number of blocks for cells
  INTEGER :: npromz_c      !< length of last block line

  ! Local scalars:

  INTEGER  :: jb, jc, jk   !< loop indices
  REAL(wp) :: rovcp        !< R/cp

  IF (ltimer) CALL timer_start(timer_con_t2theta)
  
  nblks_c   = pt_patch%nblks_int_c
  npromz_c  = pt_patch%npromz_int_c

  rovcp = rd/cpd

  ! Diagnose p and delta p (copied from m_dyn), and convert temp into theta

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,nlen) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1,nblks_c

    IF (jb /= nblks_c) THEN
       nlen = nproma
    ELSE
       nlen = npromz_c
    ENDIF

    CALL half_level_pressure( pt_prog%pres_sfc(:,jb), nproma, nlen,  &! in
                              pt_diag%pres_ic(:,:,jb)               ) ! out

    DO jk = 1, nplev
      DO jc = 1, nlen
        pt_diag%delp_c(jc,jk,jb) = delpr(jk)
      END DO
    END DO


    DO jk = nplvp1, nlev
      DO jc = 1, nlen
       pt_diag%delp_c(jc,jk,jb) = pt_diag%pres_ic(jc,jk+1,jb) - &
                                      pt_diag%pres_ic(jc,jk,jb)
      END DO
    END DO


    CALL full_level_pressure( pt_diag%pres_ic(:,:,jb), nproma, nlen, &! in
                              pt_diag%pres_mc(:,:,jb)              )  ! out

    pt_prog%theta(1:nlen,1:nlev,jb) = pt_prog%temp(1:nlen,1:nlev,jb)*   &
      EXP(rovcp*LOG(1.e5_wp/pt_diag%pres_mc(1:nlen,1:nlev,jb)))*        &
      pt_diag%delp_c(1:nlen,1:nlev,jb)

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  
  IF (ltimer) CALL timer_stop(timer_con_t2theta)

END SUBROUTINE convert_t2theta



  !-------------------------------------------------------------------------
  !>
  !!               Converts potential temperature times delta p into temperature.
  !!
  !!
  !! @par Revision History
  !!   Created by Guenther Zaengl 2009-01-13
  !!
  SUBROUTINE convert_theta2t( pt_patch, pt_prog, pt_diag)
!
  IMPLICIT NONE

  TYPE(t_patch),           INTENT(IN)    :: pt_patch  !< grid/patch info.
  TYPE(t_hydro_atm_prog),INTENT(INOUT) :: pt_prog   !< prognostic variables
  TYPE(t_hydro_atm_diag),INTENT(INOUT) :: pt_diag   !< diagnostic variables

  ! Local array bounds:

  INTEGER :: nlen          ! arb. field length
  INTEGER :: nblks_c       ! number of blocks for cells / edges
  INTEGER :: npromz_c      ! length of last block line

  ! Local scalars:

  INTEGER  :: jb, jc, jk           ! loop indices
  REAL(wp) :: rovcp                ! R/cp

  IF (ltimer) CALL timer_start(timer_con_theta2t)
  
  nblks_c   = pt_patch%nblks_int_c
  npromz_c  = pt_patch%npromz_int_c

  rovcp = rd/cpd

  ! Diagnose p and delta p (copied from m_dyn), and convert theta into temp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,nlen) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1,nblks_c

    IF (jb /= nblks_c) THEN
       nlen = nproma
    ELSE
       nlen = npromz_c
    ENDIF

    CALL half_level_pressure( pt_prog%pres_sfc(:,jb), nproma, nlen,  &! in
                              pt_diag%pres_ic(:,:,jb)               ) ! out

    DO jk = 1, nplev
      DO jc = 1, nlen
        pt_diag%delp_c(jc,jk,jb) = delpr(jk)
      END DO
    END DO


    DO jk = nplvp1, nlev
      DO jc = 1, nlen
       pt_diag%delp_c(jc,jk,jb) = pt_diag%pres_ic(jc,jk+1,jb) - &
                                  pt_diag%pres_ic(jc,jk,jb)
      END DO
    END DO

    CALL full_level_pressure( pt_diag%pres_ic(:,:,jb), nproma, nlen, &! in
                              pt_diag%pres_mc(:,:,jb)              )  ! out

    pt_prog%temp(1:nlen,1:nlev,jb) = (pt_prog%theta(1:nlen,1:nlev,jb)/ &
      pt_diag%delp_c(1:nlen,1:nlev,jb))*                               &
      EXP(rovcp*LOG(pt_diag%pres_mc(1:nlen,1:nlev,jb)/p0ref))

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  
  IF (ltimer) CALL timer_stop(timer_con_theta2t)

END SUBROUTINE convert_theta2t


  !-------------------------------------------------------------------------
  !>
  !!               Converts temperature into potential temperature times delta p
  !!               Linearized version for efficiency improvement
  !!
  !! @par Revision History
  !!   Created by Guenther Zaengl 2009-02-11
  !!
  SUBROUTINE convert_t2theta_lin( pt_patch, pt_prog, pt_diag, temp_save)
!
  IMPLICIT NONE

  TYPE(t_patch),       INTENT(IN)    :: pt_patch !< grid/patch info.
  TYPE(t_hydro_atm_prog),INTENT(INOUT) :: pt_prog  !< prognostic variables
  TYPE(t_hydro_atm_diag),INTENT(INOUT) :: pt_diag  !< diagnostic variables

  REAL(wp), DIMENSION ( nproma, nlev, pt_patch%nblks_c ),  INTENT(IN) :: temp_save

  ! Local array bounds:

  INTEGER :: nlen          ! arb. field length
  INTEGER :: nblks_c       ! number of blocks for cells / edges
  INTEGER :: npromz_c      ! length of last block line

  ! Local scalars:

  INTEGER  :: jb, jc, jk           ! loop indices
  REAL(wp) :: rovcp                ! R/cp

  REAL(wp), DIMENSION ( nproma, nlevp1, pt_patch%nblks_c ) :: pres_old
  
  IF (ltimer) CALL timer_start(timer_con_l_t2theta)

  nblks_c   = pt_patch%nblks_int_c
  npromz_c  = pt_patch%npromz_int_c

  rovcp = rd/cpd

  ! Diagnose p and delta p (copied from m_dyn), and convert temp into theta

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,nlen) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1,nblks_c

    IF (jb /= nblks_c) THEN
       nlen = nproma
    ELSE
       nlen = npromz_c
    ENDIF

    pres_old(1:nlen,1:nlevp1,jb) = pt_diag%pres_ic_new(1:nlen,1:nlevp1,jb)

    CALL half_level_pressure( pt_prog%pres_sfc(:,jb), nproma, nlen,  &! in
                              pt_diag%pres_ic_new(:,:,jb)           ) ! out

    ! Only delp is updated here because rdelp is required to carry
    ! the value corresponding to pt_prog_old

    DO jk = nplvp1, nlev
      DO jc = 1, nlen
       pt_diag%delp_c_new(jc,jk,jb) = pt_diag%pres_ic_new(jc,jk+1,jb) - &
                                      pt_diag%pres_ic_new(jc,jk,jb)
      END DO
    END DO


    pt_prog%theta(1:nlen,1:nlev,jb) = ( pt_prog%theta(1:nlen,1:nlev,jb)*       &
      pt_diag%rdelp_c_new(1:nlen,1:nlev,jb) + (pt_prog%temp(1:nlen,1:nlev,jb) -    &
      temp_save(1:nlen,1:nlev,jb)*(1._wp + rovcp - rovcp*                      &
      (pres_old(1:nlen,1:nlev,jb) + pres_old(1:nlen,2:nlevp1,jb)) /            &
      (pt_diag%pres_ic_new(1:nlen,1:nlev,jb)+pt_diag%pres_ic_new(1:nlen,2:nlevp1,jb))))/ &
      pt_diag%exner(1:nlen,1:nlev,jb)) * pt_diag%delp_c_new(1:nlen,1:nlev,jb)

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  IF (ltimer) CALL timer_stop(timer_con_l_t2theta)

END SUBROUTINE convert_t2theta_lin



  !-------------------------------------------------------------------------
  !>
  !!               Converts potential temperature times delta p into temperature
  !!               Linearized version for efficiency improvement
  !!
  !! @par Revision History
  !!   Created by Guenther Zaengl 2009-02-11
  !!
  SUBROUTINE convert_theta2t_lin( pt_patch, pt_prog_old, pt_prog, pt_diag)
!
  IMPLICIT NONE

  TYPE(t_patch),           INTENT(IN)    :: pt_patch    !< grid/patch info.
  TYPE(t_hydro_atm_prog),INTENT(INOUT) :: pt_prog_old !< prognostic variables
  TYPE(t_hydro_atm_prog),INTENT(INOUT) :: pt_prog     !< prognostic variables
  TYPE(t_hydro_atm_diag),INTENT(INOUT) :: pt_diag     !< diagnostic variables

  ! Local array bounds:

  INTEGER :: nlen          ! arb. field length
  INTEGER :: nblks_c       ! number of blocks for cells / edges
  INTEGER :: npromz_c      ! length of last block line

  ! Local scalars:

  INTEGER  :: jb, jc, jk           ! loop indices
  REAL(wp) :: rovcp                !R/cp

  IF (ltimer) CALL timer_start(timer_con_l_theta2t)
  
  nblks_c   = pt_patch%nblks_int_c
  npromz_c  = pt_patch%npromz_int_c

  rovcp = rd/cpd

  ! Diagnose p and delta p (copied from m_dyn), and convert theta into temp
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,nlen) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1,nblks_c

    IF (jb /= nblks_c) THEN
       nlen = nproma
    ELSE
       nlen = npromz_c
    ENDIF

    CALL half_level_pressure( pt_prog%pres_sfc(:,jb), nproma, nlen,  &! in
                              pt_diag%pres_ic_new(:,:,jb)           ) ! out

    ! Only delp is updated here because rdelp is required to carry
    ! the value corresponding to pt_prog_old

    DO jk = 1, nlev
      DO jc = 1, nlen
       pt_diag%delp_c_new(jc,jk,jb) = pt_diag%pres_ic_new(jc,jk+1,jb) - &
                                      pt_diag%pres_ic_new(jc,jk,jb)
      END DO
    END DO

    pt_prog%temp(1:nlen,1:nlev,jb) = pt_prog_old%temp(1:nlen,1:nlev,jb) +             &
      (pt_prog%theta(1:nlen,1:nlev,jb)/pt_diag%delp_c_new(1:nlen,1:nlev,jb) -         &
       pt_prog_old%theta(1:nlen,1:nlev,jb)*pt_diag%rdelp_c(1:nlen,1:nlev,jb)*         &
      (1._wp + rovcp - rovcp*(pt_diag%pres_ic_new(1:nlen,1:nlev,jb) +                 &
       pt_diag%pres_ic_new(1:nlen,2:nlevp1,jb))/(pt_diag%pres_ic(1:nlen,1:nlev,jb) +  &
       pt_diag%pres_ic(1:nlen,2:nlevp1,jb))))*pt_diag%exner(1:nlen,1:nlev,jb)

    ! Now, rdelp and exner are updated for use in the back-conversion
    ! after the si correction

    pt_diag%rdelp_c_new(1:nlen,1:nlev,jb) = 1._wp/pt_diag%delp_c_new(1:nlen,1:nlev,jb)
    pt_diag%exner(1:nlen,1:nlev,jb) = pt_prog%temp(1:nlen,1:nlev,jb)/        &
      (pt_prog%theta(1:nlen,1:nlev,jb)*pt_diag%rdelp_c_new(1:nlen,1:nlev,jb))

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  IF (ltimer) CALL timer_stop(timer_con_l_theta2t)

END SUBROUTINE convert_theta2t_lin


END MODULE mo_expensive_functions

