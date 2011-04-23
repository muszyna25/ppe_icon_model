!>
!! @brief Second half of the driver routine for turbulent mixing.
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI-M (2011-04)
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
#ifdef __xlC__
@PROCESS HOT
#endif

MODULE mo_vdiff_upward_sweep

  USE mo_kind,               ONLY: wp
  USE mo_vdiff_solver,       ONLY: nvar_vdiff, nmatrix, ih, iqv,         &
                                 & sfc_solve, rhs_bksub, vdiff_tendencies
#ifdef __ICON__
  USE mo_echam_vdiff_params, ONLY: tpfac2, itop, lsfc_heat_flux
#else
  USE mo_physc2,             ONLY: tpfac2, itop
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: vdiff_up

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
  !!
  SUBROUTINE vdiff_up( kproma, kbdim, klev, klevm1, klevp1, ktrac,       &! in
                       ksfc_type, idx_wtr, idx_ice, idx_lnd, idx_gbm,    &! in
                       pdtime, pstep_len, pfrc, pocu, pocv,              &! in
                       pcfm_sfc,   pcfh_sfc,    pqsat_sfc,               &! in 
                       aa, aa_btm, pprfac_sfc,  pcpt_sfc,                &! in
                       ihpbl,      pcptgz,      prhoh,       pqshear,    &! in
                       pum1,       pvm1,        ptm1,        pqm1,       &! in
                       pxlm1,      pxim1,       pxtm1,                   &! in
                       pdelpm1,    pgeom1,      pztkevn,                 &! in
#ifdef __ICON__
                       ptkem1,                                           &! in
#else
                       ptkem1, ptkem0,                                   &! in/inout
#endif
                       bb, bb_btm,                                       &! inout
                       pzthvvar,   pxvar,       pz0m,        pkedisp,    &! inout
                       pute,       pvte,        ptte,        pqte,       &! inout
                       pxlte,      pxite,       pxtte,                   &! inout
                       pute_vdf,   pvte_vdf,    ptte_vdf,                &! out
                       pqte_vdf,   pxlte_vdf,   pxite_vdf,   pxtte_vdf,  &! out
                       pxvarprod,  pvmixtau,    pqv_mflux_sfc,           &! out
                       pthvvar,    pthvsig,     ptke                     )! out

    INTEGER, INTENT(IN) :: kproma, kbdim, klev, klevm1, klevp1, ktrac
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr, idx_ice, idx_lnd, idx_gbm
    REAL(wp),INTENT(IN) :: pdtime, pstep_len
  
    REAL(wp),INTENT(IN) :: &
      & pfrc      (kbdim,ksfc_type) ,&!< area fraction of each surface type
      & pocu      (kbdim)           ,&!< eastward  velocity of ocean sfc current
      & pocv      (kbdim)           ,&!< northward velocity of ocean sfc current
      & pcfm_sfc  (kbdim,ksfc_type) ,&!< exchange coeff
      & pcfh_sfc  (kbdim,ksfc_type) ,&!< exchange coeff for heat and tracers
      & pqsat_sfc (kbdim,ksfc_type)   !< surface specific humidity at saturation

    REAL(wp),INTENT(IN) :: aa    (kbdim,klev,3,nmatrix) !< for all variables
    REAL(wp),INTENT(IN) :: aa_btm(kbdim,3,ksfc_type)    !< for heat and moisture

    REAL(wp),INTENT(IN) :: pprfac_sfc (kbdim)    !< prefactor for the exchange coefficients
    REAL(wp),INTENT(IN) :: pcpt_sfc(kbdim,ksfc_type) !< dry static energy at surface

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
  
    REAL(wp),INTENT(IN) :: pdelpm1(kbdim,klev)   !< layer thickness [Pa]
    REAL(wp),INTENT(IN) :: pgeom1 (kbdim,klev)   !< geopotential above ground
    REAL(wp),INTENT(IN) :: pztkevn(kbdim,klev)   !< intermediate value of tke
#ifdef __ICON__
    REAL(wp),INTENT(IN) :: ptkem1(kbdim,klev)    !< TKE at time step t-dt 
#else
    REAL(wp),INTENT(INOUT) :: ptkem1(kbdim,klev)
    REAL(wp),INTENT(INOUT) :: ptkem0(kbdim,klev)
#endif

    REAL(wp),INTENT(INOUT) :: bb    (kbdim,klev,nvar_vdiff)  !<
    REAL(wp),INTENT(INOUT) :: bb_btm(kbdim,ksfc_type,ih:iqv) !< for heat and moisture

    REAL(wp),INTENT(INOUT) :: pzthvvar (kbdim,klev) !< intermediate value of thvvar
    REAL(wp),INTENT(INOUT) :: pxvar    (kbdim,klev) !< distribution width (b-a) 
                                                    !< in: step t-dt, out: modified 
                                                    !< due to vertical diffusion

    ! Roughness length
    ! In: values over each surface type computed during the previous time step;
    ! Out: z0m over the ocean is updated using the time average ("hat" value)
    ! of u and v and the t-dt value of some other variables. The grid-box-mean
    ! is then updated.

    REAL(wp),INTENT(INOUT) :: pz0m (kbdim,idx_gbm:ksfc_type)

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
  
    REAL(wp),INTENT(OUT) :: pute_vdf (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pvte_vdf (kbdim,klev)
    REAL(wp),INTENT(OUT) :: ptte_vdf (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pqte_vdf (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxlte_vdf(kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxite_vdf(kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxtte_vdf(kbdim,klev,ktrac)
  
    ! Some other diagnostics
  
    REAL(wp),INTENT(OUT) :: pxvarprod    (kbdim,klev) !< shear production 
                                                      !< of the variance of total water
    REAL(wp),INTENT(OUT) :: pvmixtau     (kbdim,klev) !< vertical mixing time scale
    REAL(wp),INTENT(OUT) :: pqv_mflux_sfc(kbdim)      !< surface mass flux of water vapour
    REAL(wp),INTENT(OUT) :: pthvvar      (kbdim,klev) !< variance of virtual potential temperature 
                                                      !< at the new time step t
    REAL(wp),INTENT(OUT) :: pthvsig      (kbdim)      !< sqrt( variance of theta_v )
  
    REAL(wp),INTENT(OUT) :: ptke  (kbdim,klev)
  
    !-----------------------------------------------------------------------
    ! 5. Handle fractional surfaces (tiles) and obtain solution of 
    !    static energy and water vapor concentration in the lowest model 
    !    layer (the klev-th full level). 
    !    This could be the place where land surface comes in.
    !-----------------------------------------------------------------------
    !    - boundary condition at the air-sea/ice/land interface
  
    ! lnd/oce/ice model should provide zcpt_sfc and pqsat_sfc at new 
    ! time step as well.
  
    ! (lnd model takes the old values as input)
  
    CALL sfc_solve( kproma, kbdim, klev, klevm1,    &! in
                  & ksfc_type, idx_wtr, idx_ice,    &! in
                  & lsfc_heat_flux, tpfac2, pfrc,   &! in
                  & pocu, pocv, pcpt_sfc, pqsat_sfc,&! in
                  & pcfh_sfc, pprfac_sfc,           &! in
                  & aa, aa_btm,                     &! in
                  & bb, bb_btm                      )! inout
  
    !-----------------------------------------------------------------------
    ! 6. Obtain solution of the tri-diagonal system by back-substitution. 
    !    Then compute tendencies and diagnose moisture flux etc.
    !-----------------------------------------------------------------------
    CALL rhs_bksub( kproma, kbdim, itop, klev, aa, bb ) ! in,...,in, inout
  
    CALL vdiff_tendencies( kproma, kbdim, itop, klev, klevm1, klevp1,   &! in
                         & ktrac, ksfc_type, idx_lnd, idx_wtr,          &! in
                         & idx_gbm, pdtime, pstep_len,                  &! in
                         & pum1, pvm1, ptm1, pqm1, pxlm1, pxim1,        &! in
                         & pxtm1, pgeom1, pdelpm1, pcptgz,              &! in
#ifdef __ICON__
                         & ptkem1, pztkevn, pzthvvar, prhoh,            &! in
#else
                         & ptkem1, ptkem0, pztkevn, pzthvvar, prhoh,    &! in
#endif
                         & pqshear, ihpbl, pcfh_sfc, pqsat_sfc,         &! in
                         & pcfm_sfc, pfrc, bb,                          &! in
                         & pkedisp(:),                                  &! inout ("pvdis" in echam)
                         & pxvar(:,:), pz0m(:,:),                       &! inout
                         & pute, pvte, ptte, pqte, pxlte, pxite, pxtte, &! inout
                         & pute_vdf, pvte_vdf, ptte_vdf, pqte_vdf,      &! out
                         & pxlte_vdf, pxite_vdf, pxtte_vdf,             &! out
                         & pxvarprod,                                   &! out ("pvdiffp" in echam)
                         & ptke, pthvvar, pthvsig, pvmixtau,            &! out
                         & pqv_mflux_sfc                                )! out ("pqhfla" in echam)
  
    ! Note: computation of additional diagnostics, e.g., surface sensible heat flux,
    !       wind stress, 10m wind, 2m temperature etc., has not been implemented yet.

  END SUBROUTINE vdiff_up
  !-------------

END MODULE mo_vdiff_upward_sweep

