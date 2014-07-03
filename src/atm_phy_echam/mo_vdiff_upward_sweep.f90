!>
!! @brief Second half of the driver routine for turbulent mixing.
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI-M (2011-04)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
#ifdef __xlC__
@PROCESS HOT
#endif

MODULE mo_vdiff_upward_sweep

  USE mo_kind,               ONLY: wp
  USE mo_vdiff_solver,       ONLY: nvar_vdiff, nmatrix, ih, iqv,     &
                                 & rhs_bksub, vdiff_tendencies
#ifdef __ICON__
  USE mo_echam_vdiff_params, ONLY: tpfac2, itop
#else
  USE mo_physc2,             ONLY: tpfac2, itop
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: vdiff_up

CONTAINS
  !>
  !!
  !!
  SUBROUTINE vdiff_up( lsfc_heat_flux,                                   &! in
                       kproma, kbdim, klev, klevm1, klevp1, ktrac,       &! in
                       ksfc_type, idx_wtr, idx_ice, idx_lnd,             &! in
                       pdtime, pstep_len, pfrc,                          &! in
                       pcfm_tile,  pcfh_tile,   pqsat_tile,              &! in 
                       aa,                                               &! in
                       ihpbl,      pcptgz,      prhoh,       pqshear,    &! in
                       pum1,       pvm1,        ptm1,        pqm1,       &! in
                       pxlm1,      pxim1,       pxtm1,                   &! in
                       pcd,        pcv,                                  &! in
                       pdelpm1,    pgeom1,      pztkevn,                 &! in
#ifdef __ICON__
                       ptkem1,                                           &! in
#else
                       ptkem1, ptkem0,                                   &! in/inout
#endif
                       ptte_corr,                                        &! in
                       bb,                                               &! inout
                       pzthvvar,   pxvar,       pz0m_tile,   pkedisp,    &! inout
                       pute,       pvte,        ptte,        pqte,       &! inout
                       pxlte,      pxite,       pxtte,                   &! inout
                       pute_vdf,   pvte_vdf,    ptte_vdf,                &! out
                       pqte_vdf,   pxlte_vdf,   pxite_vdf,   pxtte_vdf,  &! out
                       pxvarprod,  pvmixtau,    pz0m,                    &! out
                       pthvvar,    pthvsig,     ptke,                    &! out
                       psh_vdiff,  pqv_vdiff                             )! out

    LOGICAL, INTENT(IN) :: lsfc_heat_flux
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, klevm1, klevp1, ktrac
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr, idx_ice, idx_lnd
    REAL(wp),INTENT(IN) :: pdtime, pstep_len
  
    REAL(wp),INTENT(IN) ::           &
      & pfrc      (kbdim,ksfc_type), &!< area fraction of each surface type
      & pcfm_tile (kbdim,ksfc_type), &!< exchange coeff
      & pcfh_tile (kbdim,ksfc_type), &!< exchange coeff for heat and tracers
      & pqsat_tile(kbdim,ksfc_type)   !< surface specific humidity at saturation

    REAL(wp),INTENT(IN) :: aa    (kbdim,klev,3,nmatrix) !< for all variables


    ! The input variables below are needed only by "vdiff_tendencies"

    INTEGER, INTENT(IN) :: ihpbl   (kbdim)
    REAL(wp),INTENT(IN) :: pcptgz  (kbdim,klev)  !< dry static energy
    REAL(wp),INTENT(IN) :: prhoh   (kbdim,klev)  !< air density at half levels
    REAL(wp),INTENT(IN) :: pqshear (kbdim,klev)
  
    REAL(wp),INTENT(IN) :: pum1    (kbdim,klev)  !< u-wind at step t-dt
    REAL(wp),INTENT(IN) :: pvm1    (kbdim,klev)  !< q-wind at step t-dt
    REAL(wp),INTENT(IN) :: ptm1    (kbdim,klev)  !< temperature at step t-dt
    REAL(wp),INTENT(IN) :: pqm1    (kbdim,klev)  !< specific humidity at step t-dt
    REAL(wp),INTENT(IN) :: pxlm1   (kbdim,klev)  !< cloud water concentration at step t-dt
    REAL(wp),INTENT(IN) :: pxim1   (kbdim,klev)  !< cloud ice   concentration at step t-dt
    REAL(wp),INTENT(IN) :: pxtm1   (kbdim,klev,ktrac) !< specific density of other tracers at step t-dt
  
    REAL(wp),INTENT(IN) :: pcd                   !< specific heat of dry air
    REAL(wp),INTENT(IN) :: pcv                   !< specific heat of water vapor

    REAL(wp),INTENT(IN) :: pdelpm1(kbdim,klev)   !< layer thickness [Pa]
    REAL(wp),INTENT(IN) :: pgeom1 (kbdim,klev)   !< geopotential above ground
    REAL(wp),INTENT(IN) :: pztkevn(kbdim,klev)   !< intermediate value of tke
#ifdef __ICON__
    REAL(wp),INTENT(IN) :: ptkem1(kbdim,klev)    !< TKE at time step t-dt 
#else
    REAL(wp),INTENT(INOUT) :: ptkem1(kbdim,klev)
    REAL(wp),INTENT(INOUT) :: ptkem0(kbdim,klev)
#endif
    REAL(wp),INTENT(IN)    :: ptte_corr(kbdim)   !< tte correction for snow melt over land

    REAL(wp),INTENT(INOUT) :: bb    (kbdim,klev,nvar_vdiff)  !<

    REAL(wp),INTENT(INOUT) :: pzthvvar (kbdim,klev) !< intermediate value of thvvar
    REAL(wp),INTENT(INOUT) :: pxvar    (kbdim,klev) !< distribution width (b-a) 
                                                    !< in: step t-dt, out: modified 
                                                    !< due to vertical diffusion

    ! Roughness length
    ! In: values over each surface type.
    ! Out: z0m_tile over the ocean is updated using the time average ("hat" value)
    ! of u and v and the t-dt value of some other variables.

    REAL(wp),INTENT(INOUT) :: pz0m_tile (kbdim,ksfc_type)

    ! Temporally and vertically integrated dissipation of kinetic energy

    REAL(wp),INTENT(INOUT) :: pkedisp(kbdim) 

    ! Tendencies

    REAL(wp),INTENT(INOUT) :: pute (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pvte (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: ptte (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pqte (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pxlte(kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pxite(kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pxtte(kbdim,klev,ktrac)
  
    REAL(wp),INTENT(INOUT) :: pute_vdf (kbdim,klev)  ! OUT
    REAL(wp),INTENT(INOUT) :: pvte_vdf (kbdim,klev)  ! OUT
    REAL(wp),INTENT(INOUT) :: ptte_vdf (kbdim,klev)  ! OUT
    REAL(wp),INTENT(INOUT) :: pqte_vdf (kbdim,klev)  ! OUT
    REAL(wp),INTENT(INOUT) :: pxlte_vdf(kbdim,klev)  ! OUT
    REAL(wp),INTENT(INOUT) :: pxite_vdf(kbdim,klev)  ! OUT
    REAL(wp),INTENT(INOUT) :: pxtte_vdf(kbdim,klev,ktrac)  ! OUT
  
    ! Some other diagnostics
  
    REAL(wp),INTENT(INOUT) :: pxvarprod    (kbdim,klev) !< shear production
                                                      !< of the variance of total water    out
    REAL(wp),INTENT(INOUT) :: pvmixtau     (kbdim,klev) !< vertical mixing time scale    out
    REAL(wp),INTENT(INOUT) :: pz0m         (kbdim)      !< grid-box mean roughness height    out
    REAL(wp),INTENT(INOUT) :: pthvvar      (kbdim,klev) !< variance of virtual potential temperature
                                                      !< at the new time step t    out
    REAL(wp),INTENT(INOUT) :: pthvsig      (kbdim)      !< sqrt( variance of theta_v )    out
    REAL(wp),INTENT(INOUT) :: ptke       (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: psh_vdiff (kbdim)         ! sens. heat flux
    REAL(wp),INTENT(INOUT) :: pqv_vdiff (kbdim)         ! qv flux

  
    !-----------------------------------------------------------------------
    ! 6. Obtain solution of the tri-diagonal system by back-substitution. 
    !    Then compute tendencies and diagnose moisture flux etc.
    !-----------------------------------------------------------------------
    CALL rhs_bksub( kproma, kbdim, itop, klev, aa, bb ) ! in,...,in, inout
  
    CALL vdiff_tendencies( kproma, kbdim, itop, klev, klevm1, klevp1,   &! in
                         & ktrac, ksfc_type, idx_lnd, idx_wtr, idx_ice, &! in
                         & pdtime, pstep_len,                           &! in
                         & pum1, pvm1, ptm1, pqm1, pxlm1, pxim1,        &! in
                         & pxtm1, pgeom1, pdelpm1, pcptgz,              &! in
                         & pcd, pcv,                                    &! in
#ifdef __ICON__
                         & ptkem1, pztkevn, pzthvvar, prhoh,            &! in
#else
                         & ptkem1, ptkem0, pztkevn, pzthvvar, prhoh,    &! in
#endif
                         & pqshear, ihpbl, pcfh_tile, pqsat_tile,       &! in
                         & pcfm_tile, pfrc, ptte_corr, bb,              &! in
                         & pkedisp(:),                                  &! inout ("pvdis" in echam)
                         & pxvar(:,:), pz0m_tile(:,:),                  &! inout
                         & pute, pvte, ptte, pqte, pxlte, pxite, pxtte, &! inout
                         & pute_vdf, pvte_vdf, ptte_vdf, pqte_vdf,      &! out
                         & pxlte_vdf, pxite_vdf, pxtte_vdf,             &! out
                         & pxvarprod,                                   &! out ("pvdiffp" in echam)
                         & pz0m, ptke, pthvvar, pthvsig, pvmixtau,      &
                         & psh_vdiff, pqv_vdiff                         )! out

    ! Note: computation of additional diagnostics, e.g., surface sensible heat flux,
    !       wind stress, 10m wind, 2m temperature etc., has not been implemented yet.

  END SUBROUTINE vdiff_up
  !-------------

END MODULE mo_vdiff_upward_sweep
