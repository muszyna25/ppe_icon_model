!>
!! This module contains the interface between ICONAM dynamics and Held-Suarez forcing
!!
!!
!! @author Pilar Ripodas, DWD
!!  based in mo_held_suarez_interface from Hui Wan for the Hydrostatic core
!!
!!
!! @par Revision History
!! First version for non-hydrostatic core by Pilar Ripodas, DWD (2010-09-14)
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

MODULE mo_nh_held_suarez_interface

  USE mo_kind,                  ONLY: wp
  USE mo_parallel_config,       ONLY: nproma
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_intp,                  ONLY: cells2edges_scalar
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e
  USE mo_hs_test,               ONLY: held_suarez_forcing_temp, &
                                      & held_suarez_forcing_vn
  USE mo_nh_diagnose_pres_temp, ONLY: diagnose_pres_temp
  USE mo_physical_constants,    ONLY: rd_o_cpd
!!$  USE mo_impl_constants,        ONLY: MAX_CHAR_LENGTH
  USE mo_nh_testcases,          ONLY: lhs_fric_heat
  USE mo_timer,                 ONLY: ltimer, timer_start, timer_stop, timer_held_suarez_intr


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC  :: held_suarez_nh_interface

CONTAINS

  !>
  !! SUBROUTINE held_suarez_nh_interface -- interface between ICONAN dynamics
  !! and Held-Suarez forcing
  !!
  !!
  !! @par Revision History
  !! fisrt release by Pilar Ripodas, DWD (2010-09-14)
  !!
  SUBROUTINE held_suarez_nh_interface (p_nh_prog,p_patch,p_int_state,p_metrics,p_nh_diag)

    TYPE(t_patch),TARGET,INTENT(in):: p_patch    !< single patch
    TYPE(t_int_state),INTENT(in)  :: p_int_state!< single interpolation state
    TYPE(t_nh_metrics),INTENT(in) :: p_metrics  !< single metrics state
    TYPE(t_nh_prog), INTENT(inout)   :: p_nh_prog  !< single nh prognostic state
    TYPE(t_nh_diag), TARGET, INTENT(inout):: p_nh_diag  !< single nh diagnostic state

    ! Local scalar

    INTEGER :: jk, jb, jbs, is, ie
    INTEGER :: nblks_c, nblks_e
    INTEGER :: nlev              !< number of full levels

    ! Local arrays

    REAL(wp) :: zlat(nproma)      ! latitude

    REAL(wp) ::   & !< pressure @ cells
      &  zsigma_mc(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) ::   & !< pressure @ edges
      &  zsigma_me(nproma,p_patch%nlev,p_patch%nblks_e)

    REAL(wp) ::   & !< tendency of temp due to HS forcing
      &  ddt_temp (nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp), DIMENSION(:,:,:), POINTER :: ptr_ddt_vn
    REAL(wp), DIMENSION(:,:,:), POINTER :: ptr_ddt_exner

!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine =  &
!!$                                   '(mo_nh_held_suarez_interface) held_suarez_nh_interface:'

    !-------------------------------------------------------------------------
    ! Dimension parameters related to refinement and MPI parallelisation
    IF (ltimer) CALL timer_start(timer_held_suarez_intr)

    nblks_e = p_patch%nblks_int_e
    nblks_c = p_patch%nblks_int_c

    ! number of vertical levels
    nlev = p_patch%nlev

    !-------------------------------------------------------------------------
    ! First the surface pressure, pressure and temperature must be diagnosed

    CALL diagnose_pres_temp ( p_metrics, p_nh_prog,               &
      &                       p_nh_prog, p_nh_diag,               &
      &                       p_patch,                            &
      &                       opt_calc_temp=.TRUE.,               &
      &                       opt_calc_pres=.TRUE. )

    !-------------------------------------------------------------------------

    ptr_ddt_vn    => p_nh_diag%ddt_vn_phy
    ptr_ddt_exner => p_nh_diag%ddt_exner_phy

    !-------------------------------------------------------------------------
    ! Newtonian cooling (and optionallly frictional heating due to Rayleigh friction)
    !-------------------------------------------------------------------------

    jbs = p_patch%cells%start_blk( grf_bdywidth_c+1,1 )
!$OMP PARALLEL
!$OMP WORKSHARE
    ddt_temp(:,:,:)=0.0_wp
!$OMP END WORKSHARE

!$OMP DO PRIVATE(jb,is,ie,jk,zlat) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs,nblks_c
       CALL get_indices_c( p_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1 )

       DO jk=1,nlev
         zsigma_mc(is:ie,jk,jb) = p_nh_diag%pres(is:ie,jk,jb)/p_nh_diag%pres_sfc(is:ie,jb)
       ENDDO

       zlat(is:ie) = p_patch%cells%center(is:ie,jb)%lat

       ! last 2 inputs in case of additional computation of frictional heating
       CALL held_suarez_forcing_temp( p_nh_diag%temp(:,:,jb),     &! in
                                    & p_nh_diag%pres(:,:,jb),     &! in
                                    & zsigma_mc(:,:,jb), zlat(:), &! in
                                    & nlev, nproma, is, ie,       &! in
                                    & ddt_temp(:,:,jb),           &! out
                                    & p_nh_diag%e_kinh(:,:,jb),   &! optional in
                                    & lhs_fric_heat)               ! optional in

       ! the tendency in temp must be transfromed to a tendency in the exner function
       ! For this it is assumed that the density is constant
       ptr_ddt_exner(is:ie,:,jb)=rd_o_cpd/p_nh_prog%theta_v(is:ie,:,jb)*ddt_temp(is:ie,:,jb)

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !-------------------------------------------------------------------------
    ! Rayleigh friction
    !-------------------------------------------------------------------------
    ! First interpolate sigma values from cells to edges

    CALL cells2edges_scalar( zsigma_mc,                    &! in
                           & p_patch, p_int_state%c_lin_e, &! in
                           & zsigma_me )                    ! out

    ! Now compute the velocity tendency due to friction

    jbs = p_patch%edges%start_blk( grf_bdywidth_e+1,1 )
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs,nblks_e
       CALL get_indices_e( p_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1 )

       CALL held_suarez_forcing_vn( p_nh_prog%vn(:,:,jb),  &! in
                                  & zsigma_me(:,:,jb),     &! in
                                  & nlev, nproma, is, ie,  &! in
                                  & ptr_ddt_vn(:,:,jb)  )   ! inout
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !--------------------------------------------------------------------------
    IF (ltimer) CALL timer_stop(timer_held_suarez_intr)

  END SUBROUTINE held_suarez_nh_interface

END MODULE mo_nh_held_suarez_interface

