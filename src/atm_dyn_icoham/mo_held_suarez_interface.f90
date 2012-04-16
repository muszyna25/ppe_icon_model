!>
!! Interface between ICOHAM dynamics and Held-Suarez forcing
!!
!! @author Hui Wan (MPI-M)
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
!!      violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!      copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!      an according license agreement with DWD and MPI-M.
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
MODULE mo_held_suarez_interface

  USE mo_kind,               ONLY: wp
  USE mo_icoham_dyn_types,   ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_model_domain,       ONLY: t_patch
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_intp,               ONLY: cells2edges_scalar
  USE mo_parallel_config,    ONLY: nproma
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_hs_test,            ONLY: held_suarez_forcing_temp,  &
                                 & held_suarez_forcing_vn

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: held_suarez_interface

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !! SUBROUTINE held_suarez -- the Interface between ICON dynamics and
  !! Held-Suarez forcing
  !!
  !! This subroutine is called in the time loop of the ICOHAM model.
  !! It takes the following as input:
  !! <ol>
  !! <li> information about the dynamics grid;
  !! <li> interplation coefficients;
  !! <li> prognostic and diagnostic variables of the dynamical core;
  !! </ol>
  !!
  !! The output includes tendencies of the prognostic variables as specified
  !! by Held and Suarez (1994, BAMS)
  !!
  !! Note that each call of this subroutine deals with a single grid level
  !! rather than the entire grid tree.

  SUBROUTINE held_suarez_interface( p_patch, p_int_state, &
                                  & dyn_prog, dyn_diag,   &
                                  & phy_tend              )

    ! Arguments

    TYPE(t_patch),TARGET,INTENT(IN) :: p_patch
    TYPE(t_int_state),INTENT(IN) :: p_int_state

    TYPE(t_hydro_atm_prog),INTENT(IN)    :: dyn_prog
    TYPE(t_hydro_atm_diag),INTENT(IN)    :: dyn_diag
    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: phy_tend

    ! Local scalar

    INTEGER :: jk, jb, jbs, is, ie
    INTEGER :: nblks_c, nblks_e, nlev

    ! Local arrays

    REAL(wp) :: zsigma (nproma,p_patch%nlev) ! sigma = pres/pres_sfc
    REAL(wp) :: zlat   (nproma)      ! latitude

    REAL(wp) :: zpres_me    (nproma,p_patch%nlev,p_patch%nblks_e) !< pressure @ edges
    REAL(wp) :: zpres_sfc_c (nproma,           1,p_patch%nblks_c) !< sfc. pressure @ cells
    REAL(wp) :: zpres_sfc_e (nproma,           1,p_patch%nblks_e) !< sfc. pressure @ edges

    !-------------------------------------------------------------------------
    ! Dimension parameters related to refinement and MPI parallelisation

    nblks_e = p_patch%nblks_int_e
    nblks_c = p_patch%nblks_int_c
    nlev    = p_patch%nlev

    !-------------------------------------------------------------------------
    ! Newtonian coolding
    !-------------------------------------------------------------------------

    jbs = p_patch%cells%start_blk( grf_bdywidth_c+1,1 )
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk,zsigma,zlat) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs,nblks_c
       CALL get_indices_c( p_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1 )
       DO jk=1,nlev
          zsigma(is:ie,jk) = dyn_diag%pres_mc (is:ie,jk,jb) &
                            /dyn_prog%pres_sfc(is:ie,   jb)
       ENDDO

       zlat(is:ie) = p_patch%cells%center(is:ie,jb)%lat

       CALL held_suarez_forcing_temp( dyn_prog%temp(:,:,jb),     &! in
                                    & dyn_diag%pres_mc(:,:,jb),  &! in
                                    & zsigma(:,:), zlat(:),      &! in
                                    & nlev, nproma, is, ie,      &! in
                                    & phy_tend%temp(:,:,jb)   )   ! inout
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !-------------------------------------------------------------------------
    ! Rayleigh friction
    !-------------------------------------------------------------------------
    ! First interpolate surface pressure from cells to edges in order to
    ! compute sigma at edge centers

    CALL cells2edges_scalar( dyn_diag%pres_mc,             &! in
                           & p_patch, p_int_state%c_lin_e, &! in
                           & zpres_me )                     ! out

    ! Interpolate surface pressure

!$OMP PARALLEL WORKSHARE
    zpres_sfc_c(:,1,:) = dyn_prog%pres_sfc(:,:)
!$OMP END PARALLEL WORKSHARE

    CALL cells2edges_scalar( zpres_sfc_c,                  &! in
                           & p_patch, p_int_state%c_lin_e, &! in
                           & zpres_sfc_e,                  &! out
                           & 1,1    ) ! start and end indices of vertical layer

    ! Now compute the velocity tendency due to friction

    jbs = p_patch%edges%start_blk( grf_bdywidth_e+1,1 )
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk,zsigma) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs,nblks_e
       CALL get_indices_e( p_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1 )

       DO jk=1,nlev
          zsigma(is:ie,jk) = zpres_me    (is:ie,jk,jb) &
                            /zpres_sfc_e (is:ie, 1,jb)
       ENDDO

       CALL held_suarez_forcing_vn( dyn_prog%vn(:,:,jb),  &! in
                                  & zsigma(:,:),          &! in
                                  & nlev, nproma, is, ie, &! in
                                  & phy_tend%vn(:,:,jb)  ) ! inout
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  !--------------------------
  END SUBROUTINE held_suarez_interface

END MODULE mo_held_suarez_interface
