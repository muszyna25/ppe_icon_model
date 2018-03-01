!>
!! @brief Subroutine echam_phy_main calls all the parameterization schemes
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
!!  Modified for ICONAM by Marco Giorgetta (2014)
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
#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
@PROCESS SPILLSIZE(5000)
#endif
!OCL NOALIAS

MODULE mo_interface_echam_vdiff_surface

  !$ser verbatim USE mo_ser_echam_vdiff_down, ONLY: serialize_vdiff_down_input => serialize_input, &
  !$ser&                                            serialize_vdiff_down_output => serialize_output
  !$ser verbatim USE mo_ser_echam_vdiff_up, ONLY: serialize_vdiff_up_input => serialize_input, &
  !$ser&                                            serialize_vdiff_up_output => serialize_output
  !$ser verbatim USE mo_ser_echam_surface, ONLY: serialize_surface_input => serialize_input, &
  !$ser&                                         serialize_surface_output => serialize_output
  USE mo_kind                ,ONLY: wp

  USE mo_model_domain        ,ONLY: t_patch
  USE mo_loopindices         ,ONLY: get_indices_c

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config          ,ONLY: ntracer, nlev, nlevm1, nlevp1, iqv, iqc, iqi, iqt
  USE mo_mpi_phy_config      ,ONLY: mpi_phy_config

  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, t_echam_phy_tend
  USE mo_echam_sfc_indices   ,ONLY: nsfc_type, iwtr, iice, ilnd

  USE mo_vdiff_downward_sweep,ONLY: vdiff_down
  USE mo_vdiff_upward_sweep  ,ONLY: vdiff_up
  USE mo_vdiff_solver        ,ONLY: nvar_vdiff, nmatrix, imh, imqv, ih_vdiff=>ih, iqv_vdiff=>iqv
  
  USE mo_surface             ,ONLY: update_surface
  USE mo_surface_diag        ,ONLY: nsurf_diag

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop,  &
       &                            timer_vdiff_down, timer_vdiff_up, &
       &                            timer_surface
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_vdiff_surface

CONTAINS

  !-------------------------------------------------------------------
  SUBROUTINE interface_echam_vdiff_surface(is_in_sd_ed_interval,    &
       &                                   is_active,               &
       &                                   patch, rl_start, rl_end, &
       &                                   field, tend,             &
       &                                   pdtime                   )

    LOGICAL                 ,INTENT(in)    :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in)    :: is_active
    TYPE(t_patch)   ,TARGET ,INTENT(in)    :: patch
    INTEGER                 ,INTENT(in)    :: rl_start, rl_end
    TYPE(t_echam_phy_field) ,POINTER       :: field    
    TYPE(t_echam_phy_tend)  ,POINTER       :: tend
    REAL(wp)                ,INTENT(in)    :: pdtime

    INTEGER  :: jg             
    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block
    
    REAL(wp) :: zxt_emis(nproma,ntracer-iqt+1)  !< tracer tendency due to surface emission
                                                !< and dry deposition. "zxtems" in ECHAM5

    jg         = patch%id
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    !-------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(jcs,jce, zxt_emis)
    DO jb = i_startblk,i_endblk
       !
       CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
       !
!!$       ! Emission of aerosols or other tracers (not implemented yet)
!!$       IF (ntrac>0) THEN
!!$          CALL tracer_emission()
!!$       ENDIF
       zxt_emis(jcs:jce,:) = 0._wp
       !
!!$       ! Dry deposition of aerosols or other tracers (not implemented yet)
!!$       CALL dry_deposition()
       !
       CALL echam_vdiff_down_surf_up(is_in_sd_ed_interval,          &
            &                        is_active,                     &
            &                        jg, jb,jcs,jce, nproma,        &
            &                        field,  tend,                  &
            &                        zxt_emis,                      &
            &                        pdtime                         )
       !
    END DO
!$OMP END PARALLEL DO 
    !-------------------------------------------------------------------

  END SUBROUTINE interface_echam_vdiff_surface
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  SUBROUTINE echam_vdiff_down_surf_up(is_in_sd_ed_interval,  &
       &                              is_active,             &
       &                              jg, jb,jcs,jce, nbdim, &
       &                              field,  tend,          &
       &                              zxt_emis,              &
       &                              pdtime                 )

    LOGICAL                 ,INTENT(in)    :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in)    :: is_active
    INTEGER                 ,INTENT(in)    :: jg             
    INTEGER                 ,INTENT(in)    :: jb                  !< block index
    INTEGER                 ,INTENT(in)    :: jcs, jce            !< start/end column index within this block
    INTEGER                 ,INTENT(in)    :: nbdim               !< size of this block
    TYPE(t_echam_phy_field) ,POINTER       :: field
    TYPE(t_echam_phy_tend)  ,POINTER       :: tend
    REAL(wp)                ,INTENT(in)    :: zxt_emis(nbdim,ntracer-iqt+1)  !< tracer tendency due to surface emission
    REAL(wp)                ,INTENT(in)    :: pdtime

    ! local
    REAL(wp) :: zcpt_sfc_tile(nbdim,nsfc_type)  !< dry static energy at surface
   
    ! Coefficient matrices and right-hand-side vectors for the turbulence solver
    ! _btm refers to the lowest model level (i.e., full level "klev", not the surface)
    REAL(wp) :: zri_tile(nbdim,nsfc_type)           !< Richardson number
    REAL(wp) :: zaa    (nbdim,nlev,3,nmatrix)       !< coeff. matrices, all variables
    REAL(wp) :: zaa_btm(nbdim,3,nsfc_type,imh:imqv) !< last row of coeff. matrix of heat and moisture
    REAL(wp) :: zbb    (nbdim,nlev,nvar_vdiff)              !< r.h.s., all variables
    REAL(wp) :: zbb_btm(nbdim,nsfc_type,ih_vdiff:iqv_vdiff) !< last row of r.h.s. of heat and moisture

    ! Temporary arrays used by VDIFF, JSBACH
    REAL(wp) :: zfactor_sfc(nbdim)

    REAL(wp) :: zcptgz   (nbdim,nlev)        !< dry static energy
    REAL(wp) :: zthvvar  (nbdim,nlev)        !< intermediate value of thvvar
    REAL(wp) :: ztkevn   (nbdim,nlev)        !< intermediate value of tke
    REAL(wp) :: zch_tile (nbdim,nsfc_type)   !< for "nsurf_diag"
!!$    REAL(wp) :: zchn_tile(nbdim,nsfc_type)   !< for "nsurf_diag"
!!$    REAL(wp) :: zcdn_tile(nbdim,nsfc_type)   !< for "nsurf_diag"
!!$    REAL(wp) :: zcfnc_tile(nbdim,nsfc_type)  !< for "nsurf_diag"
    REAL(wp) :: zbn_tile (nbdim,nsfc_type)   !< for "nsurf_diag"
    REAL(wp) :: zbhn_tile(nbdim,nsfc_type)   !< for "nsurf_diag"
    REAL(wp) :: zbm_tile (nbdim,nsfc_type)   !< for "nsurf_diag"
    REAL(wp) :: zbh_tile (nbdim,nsfc_type)   !< for "nsurf_diag"

    REAL(wp) :: zq_snocpymlt(nbdim)          !< heating by melting of snow on the canopy [W/m2]
    !                                        !  which warms the lowermost atmospheric layer (JSBACH)

    INTEGER  :: ntrac

    ntrac = ntracer-iqt+1  !# of tracers excluding water vapour and hydrometeors

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          IF (ltimer) CALL timer_start(timer_vdiff_down)
          !
          ! Turbulent mixing, part I:
          ! - computation of exchange coefficients in the atmosphere and at the surface;
          ! - build up the tridiagonal linear algebraic system;
          ! - downward sweep (Gaussian elimination from top till level nlev-1)
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 input fields serialization
          !$ser verbatim call serialize_vdiff_down_input(jb, jce, nbdim, nlev,&
          !$ser verbatim   nlevm1, nlevp1, ntrac, nsfc_type, iwtr, iice, ilnd,&
          !$ser verbatim   pdtime, field)
          !
          CALL vdiff_down(jce, nbdim, nlev, nlevm1, nlevp1,&! in
               &          ntrac, nsfc_type,                &! in
               &          iwtr, iice, ilnd,                &! in, indices of different surface types
               &          pdtime,                          &! in, time step
               &          field%coriol(:,jb),              &! in, Coriolis parameter
               &          field%   zf(:,:,jb),             &! in, geopot. height above sea level, full level
               &          field%   zh(:,:,jb),             &! in, geopot. height above sea level, half level
               &          field%frac_tile(:,jb,:),         &! in, area fraction of each sfc type
               &          field% ts_tile(:,jb,:),          &! in, surface temperature
               &          field% ocu (:,jb),               &! in, ocean sfc velocity, u-component
               &          field% ocv (:,jb),               &! in, ocean sfc velocity, v-component
               &          field% presi_old(:,nlevp1,jb),   &! in, sfc pressure
               &          field%   ua(:,:,jb),             &! in, um1
               &          field%   va(:,:,jb),             &! in, vm1
               &          field%   ta(:,:,jb),             &! in, tm1
               &          field% qtrc(:,:,jb,iqv),         &! in, qm1
               &          field% qtrc(:,:,jb,iqc),         &! in, xlm1
               &          field% qtrc(:,:,jb,iqi),         &! in, xim1
               &          field%   qx(:,:,jb),             &! in, xlm1 + xim1
               &          field% qtrc(:,:,jb,iqt:),        &! in, xtm1
               &          field% mair(:,:,jb),             &! in,     air mass
               &          field% mdry(:,:,jb),             &! in, dry air mass
               &          field% presi_old(:,:,jb),        &! in, aphm1
               &          field% presm_old(:,:,jb),        &! in, apm1
               &          field%   tv(:,:,jb),             &! in, virtual temperaturea
               &          field% aclc(:,:,jb),             &! in, cloud fraction
               &          zxt_emis,                        &! in, zxtems
               &          field% thvvar(:,:,jb),           &! in, variance of theta_v at step t-dt
               &          field%   xvar(:,:,jb),           &! in
               &          field% z0m_tile(:,jb,:),         &! in
               &          field%  tkem1(:,:,jb),           &! in, TKE at step t-dt
               &          field%  ustar(:,  jb),           &! inout
               &          field%  wstar(:,  jb),           &! out, convective velocity scale
               &          field%  wstar_tile(:,jb,:),      &! inout, convective velocity scale (each sfc type)
               &          field% qs_sfc_tile(:,jb,:),      &! out, sfc specific humidity at saturation
               &          field%    ghpbl(:,jb),           &! out, for output
               &          field%      ri (:,:,jb),         &! out, for output
               &          zri_tile (:,:),                  &! out, for nsurf_diag
               &          field%  mixlen (:,:,jb),         &! out, for output
               &          field% cfm     (:,:,jb),         &! out, for output
               &          field% cfm_tile(:,jb,:),         &! out, for output and "vdiff_up"
               &          field% cfh     (:,:,jb),         &! out, for output
               &          field% cfh_tile(:,jb,:),         &! out, for output and "vdiff_up"
               &          field% cfv     (:,:,jb),         &! out, for output
               &          field% cftke   (:,:,jb),         &! out, for output
               &          field% cfthv   (:,:,jb),         &! out, for output
               &          zaa, zaa_btm, zbb, zbb_btm,      &! out, for "vdiff_up"
               &          zfactor_sfc(:),                  &! out, for "vdiff_up"
               &          zcpt_sfc_tile(:,:),              &! out, for "vdiff_up"
               &          zcptgz(:,:),                     &! out, for "vdiff_up"
               &          zthvvar(:,:),                    &! out, for "vdiff_up"
               &          field%   thvsig(:,  jb),         &! out, for "cucall"
               &          ztkevn (:,:),                    &! out, for "vdiff_up"
               &          zch_tile(:,:),                   &! out, for "nsurf_diag"
!!$               &          zchn_tile(:,:),                  &! out, for "nsurf_diag"
!!$               &          zcdn_tile(:,:),                  &! out, for "nsurf_diag"
!!$               &          zcfnc_tile(:,:),                 &! out, for "nsurf_diag"
               &          zbn_tile(:,:),                   &! out, for "nsurf_diag"
               &          zbhn_tile(:,:),                  &! out, for "nsurf_diag"
               &          zbm_tile(:,:),                   &! out, for "nsurf_diag"
               &          zbh_tile(:,:),                   &! out, for "nsurf_diag"
               &          pcsat = field% csat(:,jb),       &! in, optional, area fraction with wet land surface
               &          pcair = field% cair(:,jb),       &! in, optional, area fraction with wet land surface (air)
               &          paz0lh = field% z0h_lnd(:,jb))    ! in, optional, roughness length for heat over land
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 output fields serialization
          !$ser verbatim call serialize_vdiff_down_output(jb, jce, nbdim, nlev,&
          !$ser verbatim   nlevm1, nlevp1, ntrac, nsfc_type, iwtr, iice, ilnd,&
          !$ser verbatim   pdtime, field)
          !
          IF (ltimer) CALL timer_stop(timer_vdiff_down)
          !
          !
          ! Surface processes that provide time-dependent lower boundary
          ! condition for wind, temperature, tracer concentraion, etc.
          !
          field% lhflx_tile(jcs:jce,jb,:) = 0._wp
          field% shflx_tile(jcs:jce,jb,:) = 0._wp
          field% evap_tile (jcs:jce,jb,:) = 0._wp
          !
          IF (ltimer) CALL timer_start(timer_surface)
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 input fields serialization
          !$ser verbatim call serialize_surface_input(jb, nlev, nlevp1, iqv, jg,&
          !$ser verbatim   jce, nbdim, nlev, nsfc_type, iwtr, iice, ilnd, pdtime,&
          !$ser verbatim   field, zfactor_sfc, zaa, zaa_btm, zbb, zbb_btm,&
          !$ser verbatim   zcpt_sfc_tile, nblock = jb, pch_tile = zch_tile)
          !
          CALL update_surface(jg, jce, nbdim, field%kice,                     &! in
               &              nlev, nsfc_type,                                &! in
               &              iwtr, iice, ilnd,                               &! in, indices of surface types
               &              pdtime,                                         &! in, time step
               &              field%frac_tile(:,jb,:),                        &! in, area fraction
               &              field% cfh_tile(:,jb,:),                        &! in, from "vdiff_down"
               &              field% cfm_tile(:,jb,:),                        &! in, from "vdiff_down"
               &              zfactor_sfc(:),                                 &! in, from "vdiff_down"
               &              field% ocu (:,jb),                              &! in, ocean sfc velocity, u-component
               &              field% ocv (:,jb),                              &! in, ocean sfc velocity, v-component
               &              zaa, zaa_btm, zbb, zbb_btm,                     &! inout
               &              zcpt_sfc_tile(:,:),                             &! inout, from "vdiff_down", for "vdiff_up"
               &              field%qs_sfc_tile(:,jb,:),                      &! inout, from "vdiff_down", for "vdiff_up"
               &              field% ts_tile(:,jb,:),                         &! inout
               &              field%u_stress    (:,  jb),                     &! out
               &              field%v_stress    (:,  jb),                     &! out
               &              field% lhflx      (:,  jb),                     &! out
               &              field% shflx      (:,  jb),                     &! out
               &              field%  evap      (:,  jb),                     &! out, for "cucall"
               &              field%u_stress_tile  (:,jb,:),                  &! out
               &              field%v_stress_tile  (:,jb,:),                  &! out
               &              field% lhflx_tile    (:,jb,:),                  &! out
               &              field% shflx_tile    (:,jb,:),                  &! out
               &              field%  evap_tile    (:,jb,:),                  &! out
               &              nblock = jb,                                    &! in
               &              lsm = field%lsmask(:,jb),                       &!< in, land-sea mask
               &              alake = field%alake(:,jb),                      &! in, lake fraction
               &              pu    = field% ua(:,nlev,jb),                   &! in, um1
               &              pv    = field% va(:,nlev,jb),                   &! in, vm1
               &              ptemp = field% ta(:,nlev,jb),                   &! in, tm1
               &              pq = field% qtrc(:,nlev,jb,iqv),                &! in, qm1
               &              prsfl = field% rsfl(:,jb),                      &! in, rain surface large scale (from cloud)
               &              prsfc = field% rsfc(:,jb),                      &! in, rain surface concective (from cucall)
               &              pssfl = field% ssfl(:,jb),                      &! in, snow surface large scale (from cloud)
               &              pssfc = field% ssfc(:,jb),                      &! in, snow surface concective (from cucall)
               &              rlds        = field% rlds (:,jb),               &! in,  downward surface  longwave flux [W/m2]
               &              rlus        = field% rlus (:,jb),               &! inout, upward surface  longwave flux [W/m2]
               &              rsds        = field% rsds (:,jb),               &! in,  downward surface shortwave flux [W/m2]
               &              rsus        = field% rsus (:,jb),               &! in,  upward surface shortwave flux [W/m2]
               !
               &              rvds_dir   = field%rvds_dir   (:,jb),           &! in, all-sky downward direct visible radiation at surface
               &              rpds_dir   = field%rpds_dir   (:,jb),           &! in, all-sky downward direct PAR     radiation at surface
               &              rnds_dir   = field%rnds_dir   (:,jb),           &! in, all-sky downward direct near-IR radiation at surface
               &              rvds_dif   = field%rvds_dif   (:,jb),           &! in, all-sky downward diffuse visible radiation at surface
               &              rpds_dif   = field%rpds_dif   (:,jb),           &! in, all-sky downward diffuse PAR     radiation at surface
               &              rnds_dif   = field%rnds_dif   (:,jb),           &! in, all-sky downward diffuse near-IR radiation at surface
               !
               &              ps = field% presi_old(:,nlevp1,jb),             &! in, paphm1, half level pressure
               &              pcosmu0 = field% cosmu0(:,jb),                  &! in, amu0_x, cos of zenith angle
               &              pch_tile = zch_tile(:,:),                       &! in, from "vdiff_down" for JSBACH
               &              pcsat = field%csat(:,jb),                       &! inout, area fraction with wet land surface
               &              pcair = field%cair(:,jb),                       &! inout, area fraction with wet land surface (air)
               &              q_snocpymlt = zq_snocpymlt(:),                  &! out, heating  by melting snow on the canopy [W/m2]
               &              z0m_tile = field% z0m_tile(:,jb,:),             &! inout, roughness length for momentum over tiles
               &              z0h_lnd  = field% z0h_lnd (:,jb),               &! out, roughness length for heat over land
               &              albvisdir      = field% albvisdir     (:,jb)  , &! inout
               &              albnirdir      = field% albnirdir     (:,jb)  , &! inout
               &              albvisdif      = field% albvisdif     (:,jb)  , &! inout
               &              albnirdif      = field% albnirdif     (:,jb)  , &! inout
               &              albvisdir_tile = field% albvisdir_tile(:,jb,:), &! inout
               &              albnirdir_tile = field% albnirdir_tile(:,jb,:), &! inout
               &              albvisdif_tile = field% albvisdif_tile(:,jb,:), &! inout
               &              albnirdif_tile = field% albnirdif_tile(:,jb,:), &! inout
               &              albedo         = field% albedo        (:,jb)  , &! inout
               &              albedo_tile    = field% albedo_tile(:,jb,:),    &! inout
               &              ptsfc     = field%ts    (:,jb),                 &! out
               &              ptsfc_rad = field%ts_rad(:,jb),                 &! out
               &              rlns_tile = field%lwflxsfc_tile(:,jb,:),        &! out (for coupling)
               &              rsns_tile = field%swflxsfc_tile(:,jb,:),        &! out (for coupling)
               &              lake_ice_frc = field%lake_ice_frc(:,jb),        &! out
               &              Tsurf = field% Tsurf(:,:,jb),                   &! inout, for sea ice
               &              T1    = field% T1   (:,:,jb),                   &! inout, for sea ice
               &              T2    = field% T2   (:,:,jb),                   &! inout, for sea ice
               &              hi    = field% hi   (:,:,jb),                   &! in, for sea ice
               &              hs    = field% hs   (:,:,jb),                   &! in, for sea ice
               &              conc  = field% conc (:,:,jb),                   &! in, for sea ice
               &              Qtop  = field% Qtop (:,:,jb),                   &! out, for sea ice
               &              Qbot  = field% Qbot (:,:,jb),                   &! out, for sea ice
               &              albvisdir_ice = field% albvisdir_ice(:,:,jb),   &! inout ice albedos
               &              albnirdir_ice = field% albnirdir_ice(:,:,jb),   &! inout
               &              albvisdif_ice = field% albvisdif_ice(:,:,jb),   &! inout
               &              albnirdif_ice = field% albnirdif_ice(:,:,jb))    ! inout

          !----------------------------------------------------------------------------------------
          ! Serialbox2 output fields serialization
          !$ser verbatim call serialize_surface_output(jb, jg, jce, nbdim, nlev, nsfc_type,&
          !$ser verbatim   iwtr, iice, ilnd, pdtime, field, zaa, zaa_btm, zbb,&
          !$ser verbatim   zbb_btm, zcpt_sfc_tile, nblock = jb,&
          !$ser verbatim   q_snocpymlt = zq_snocpymlt)
          !
          IF (ltimer) CALL timer_stop(timer_surface)
          !
          !
          IF (mpi_phy_config(jg)%ljsb) THEN
             !
             ! heating accumulated
             ! zq_snocpymlt = heating for melting of snow on canopy
             !              = cooling of atmosphere
             ! --> negative sign
             field% q_phy(jcs:jce,nlev,jb) = field% q_phy(jcs:jce,nlev,jb) - zq_snocpymlt(jcs:jce)
             !
             ! tendency
             tend% ta_sfc(jcs:jce,jb)      = -zq_snocpymlt(jcs:jce) * field% qconv(jcs:jce,nlev,jb)
             !
             ! tendencies accumulated
             tend% ta_phy(jcs:jce,nlev,jb) = tend% ta_phy(jcs:jce,nlev,jb) + tend% ta_sfc(jcs:jce,jb)
             !
          END IF
          !
          !
          ! Turbulent mixing, part II:
          ! - Elimination for the lowest model level using boundary conditions
          !   provided by the surface model(s);
          ! - Back substitution to get solution of the tridiagonal system;
          ! - Compute tendencies and additional diagnostics.
          !
          IF (ltimer) CALL timer_start(timer_vdiff_up)
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 input fields serialization
          !$ser verbatim call serialize_vdiff_up_input(jb, jce, nbdim, nlev,&
          !$ser verbatim   nlevm1, ntrac, nsfc_type, iwtr, pdtime, field)
          !
          CALL vdiff_up(jce, nbdim, nlev, nlevm1,        &! in
               &        ntrac, nsfc_type,                &! in
               &        iwtr,                            &! in, indices of different sfc types
               &        pdtime,                          &! in, time steps
               &        field%frac_tile(:,jb,:),         &! in, area fraction of each sfc type
               &        field% cfm_tile(:,jb,:),         &! in
               &        zaa,                             &! in, from "vdiff_down"
               &         zcptgz(:,:),                    &! in, from "vdiff_down"
               &        field%   ua(:,:,jb),             &! in, um1
               &        field%   va(:,:,jb),             &! in, vm1
               &        field%   ta(:,:,jb),             &! in, tm1
               &        field% mair(:,:,jb),             &! in, moist air mass [kg/m2]
               &        field% mdry(:,:,jb),             &! in, dry   air mass [kg/m2]
               &        field% qtrc(:,:,jb,iqv),         &! in, qm1
               &        field% qtrc(:,:,jb,iqc),         &! in, xlm1
               &        field% qtrc(:,:,jb,iqi),         &! in, xim1
               &        field% qtrc(:,:,jb,iqt:),        &! in, xtm1
               &        field% geom(:,:,jb),             &! in, pgeom1 = geopotential above ground
               &             ztkevn(:,:),                &! in, tke at intermediate time step
               &        zbb,                             &! in
               &        zthvvar(:,:),                    &! inout
               &        field%   xvar(:,:,jb),           &! inout
               &        field% z0m_tile(:,jb,:),         &! inout
               &        field% kedisp(:,  jb),           &! out, vert. integr. diss. kin. energy [W/m2]
               &         tend%   ua_vdf(:,:,jb),         &! out
               &         tend%   va_vdf(:,:,jb),         &! out
               &        field% q_vdf   (:,:,jb),         &! out   heating W/m2
               &         tend% qtrc_vdf(:,:,jb,iqv),     &! out
               &         tend% qtrc_vdf(:,:,jb,iqc),     &! out
               &         tend% qtrc_vdf(:,:,jb,iqi),     &! out
               &         tend% qtrc_vdf(:,:,jb,iqt:),    &! out
               &        field%   z0m   (:,  jb),         &! out, for the next step
               &        field%   thvvar(:,:,jb),         &! out, for the next step
               &        field%      tke(:,:,jb),         &! out
               &        field%   sh_vdiff(:,  jb),       &! out, for energy diagnostic
               &        field%   qv_vdiff(:,  jb)        )! out, for energy diagnostic
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 input fields serialization
          !$ser verbatim call serialize_vdiff_up_output(jb, jce, nbdim, nlev,&
          !$ser verbatim   nlevm1, ntrac, nsfc_type, iwtr, pdtime, field)
          !
          IF (ltimer) CALL timer_stop(timer_vdiff_up)
          !
       END IF
       !
       ! convert    heating
       tend% ta_vdf(jcs:jce,:,jb) = field% q_vdf(jcs:jce,:,jb) * field% qconv(jcs:jce,:,jb)
       !
       ! accumulate heating
       field% q_phy(jcs:jce,:,jb) = field% q_phy(jcs:jce,:,jb) + field% q_vdf(jcs:jce,:,jb)
       !
       ! accumulated tendencies
       tend%   ua_phy(jcs:jce,:,jb)      = tend%   ua_phy(jcs:jce,:,jb)      + tend%   ua_vdf(jcs:jce,:,jb)
       tend%   va_phy(jcs:jce,:,jb)      = tend%   va_phy(jcs:jce,:,jb)      + tend%   va_vdf(jcs:jce,:,jb)
       tend%   ta_phy(jcs:jce,:,jb)      = tend%   ta_phy(jcs:jce,:,jb)      + tend%   ta_vdf(jcs:jce,:,jb)
       tend% qtrc_phy(jcs:jce,:,jb,iqv)  = tend% qtrc_phy(jcs:jce,:,jb,iqv)  + tend% qtrc_vdf(jcs:jce,:,jb,iqv)
       tend% qtrc_phy(jcs:jce,:,jb,iqc)  = tend% qtrc_phy(jcs:jce,:,jb,iqc)  + tend% qtrc_vdf(jcs:jce,:,jb,iqc)
       tend% qtrc_phy(jcs:jce,:,jb,iqi)  = tend% qtrc_phy(jcs:jce,:,jb,iqi)  + tend% qtrc_vdf(jcs:jce,:,jb,iqi)
       tend% qtrc_phy(jcs:jce,:,jb,iqt:) = tend% qtrc_phy(jcs:jce,:,jb,iqt:) + tend% qtrc_vdf(jcs:jce,:,jb,iqt:)

!!$       ! TIME FILTER FOR TURBULENT KINETIC ENERGY
!!$
!!$       IF(.NOT.lstart) THEN
!!$         zeps=eps
!!$       ELSE
!!$         zeps=0._wp
!!$       END IF
!!$       DO 397 jk=ktdia,klev
!!$         DO 396 jl=1,kproma
!!$           ptkem1(jl,jk)=ptkem(jl,jk)                                    &
!!$                     +zeps*(ptkem1(jl,jk)-2._wp*ptkem(jl,jk)+ptke(jl,jk))
!!$           ptkem(jl,jk)=ptke(jl,jk)
!!$396      END DO
!!$397    END DO

       ! 2-tl-scheme
       field% tkem1(jcs:jce,:,jb) = field% tke  (jcs:jce,:,jb)
       !
       ! Turbulent mixing, part III:
       ! - Further diagnostics.
       !
       CALL nsurf_diag(jce, nbdim, nsfc_type,           &! in
            &          ilnd,                            &! in
            &          field%frac_tile(:,jb,:),         &! in
            &          field%  qtrc(:,nlev,jb,iqv),     &! in humidity qm1
            &          field%    ta(:,nlev,jb),         &! in tm1
            &          field% presm_old(:,nlev,jb),     &! in, apm1
            &          field% presi_old(:,nlevp1,jb),   &! in, aphm1
            &          field%   qx(:,nlev,jb),          &! in, xlm1 + xim1
            &          field%   ua(:,nlev,jb),          &! in, um1
            &          field%   va(:,nlev,jb),          &! in, vm1
            &          field% ocu (:,jb),               &! in, ocean sfc velocity, u-component
            &          field% ocv (:,jb),               &! in, ocean sfc velocity, v-component
            &          field% zf  (:,nlev  ,jb),        &! in, height of lowermost full level (m)
            &          field% zh  (:,nlev+1,jb),        &! in, surface height    (m)
            &          zcptgz(:,nlev),                  &! in dry static energy
            &          zcpt_sfc_tile(:,:),              &! in dry static energy
            &          zbn_tile(:,:),                   &! in for diagnostic
            &          zbhn_tile(:,:),                  &! in for diagnostic
            &          zbh_tile(:,:),                   &! in for diagnostic
            &          zbm_tile(:,:),                   &! in for diagnostic
            &          zri_tile(:,:),                   &! in 
            &          field%sfcWind(:,  jb),           &! out 10m windspeed
            &          field%    tas(:,  jb),           &! out temperature in 2m
            &          field%   dew2(:,  jb),           &! out dew point temperature in 2m
            &          field%    uas(:,  jb),           &! out zonal wind in 10m
            &          field%    vas(:,  jb),           &! out meridional wind in 10m
            &          field%tasmax (:,  jb),           &! out max 2m temperature
            &          field%tasmin (:,  jb),           &! out min 2m temperature
            &          field%sfcWind_tile(:,jb,:),      &! out 10m windspeed on tiles
            &          field%    tas_tile(:,jb,:),      &! out temperature in 2m on tiles
            &          field%   dew2_tile(:,jb,:),      &! out dew point temperature in 2m on tiles
            &          field%    uas_tile(:,jb,:),      &! out zonal wind in 10m on tiles
            &          field%    vas_tile(:,jb,:)       )! out meridional wind in 10m on tiles
       !
       !
    ELSE
       !
       !
       ! vdiff_down
       field% ustar          (jcs:jce,  jb  ) = 0.0_wp
       field% wstar          (jcs:jce,  jb  ) = 0.0_wp
       field% wstar_tile     (jcs:jce,  jb,:) = 0.0_wp
       field% qs_sfc_tile    (jcs:jce,  jb,:) = 0.0_wp
       field% ghpbl          (jcs:jce,  jb  ) = 0.0_wp
       field% ri             (jcs:jce,:,jb  ) = 0.0_wp
       field% mixlen         (jcs:jce,:,jb  ) = 0.0_wp
       field% cfm            (jcs:jce,:,jb  ) = 0.0_wp
       field% cfm_tile       (jcs:jce,  jb,:) = 0.0_wp
       field% cfh            (jcs:jce,:,jb  ) = 0.0_wp
       field% cfh_tile       (jcs:jce,  jb,:) = 0.0_wp
       field% cfv            (jcs:jce,:,jb  ) = 0.0_wp
       field% cftke          (jcs:jce,:,jb  ) = 0.0_wp
       field% cfthv          (jcs:jce,:,jb  ) = 0.0_wp
       field% thvsig         (jcs:jce,  jb  ) = 0.0_wp
       !
       ! update_surface
       field% ts_tile        (jcs:jce,  jb,:) = 0.0_wp
       field% u_stress       (jcs:jce,  jb  ) = 0.0_wp
       field% v_stress       (jcs:jce,  jb  ) = 0.0_wp
       field% lhflx          (jcs:jce,  jb  ) = 0.0_wp
       field% shflx          (jcs:jce,  jb  ) = 0.0_wp
       field%  evap          (jcs:jce,  jb  ) = 0.0_wp
       field% u_stress_tile  (jcs:jce,  jb,:) = 0.0_wp
       field% v_stress_tile  (jcs:jce,  jb,:) = 0.0_wp
       field% lhflx_tile     (jcs:jce,  jb,:) = 0.0_wp
       field% shflx_tile     (jcs:jce,  jb,:) = 0.0_wp
       field%  evap_tile     (jcs:jce,  jb,:) = 0.0_wp
       field% rlus           (jcs:jce,  jb  ) = 0.0_wp
       field% rsus           (jcs:jce,  jb  ) = 0.0_wp
       field% csat           (jcs:jce,  jb  ) = 0.0_wp
       field% cair           (jcs:jce,  jb  ) = 0.0_wp
       field% z0h_lnd        (jcs:jce,  jb  ) = 0.0_wp
       field% albvisdir      (jcs:jce,  jb  ) = 0.0_wp
       field% albnirdir      (jcs:jce,  jb  ) = 0.0_wp
       field% albvisdif      (jcs:jce,  jb  ) = 0.0_wp
       field% albnirdif      (jcs:jce,  jb  ) = 0.0_wp
       field% albvisdir_tile (jcs:jce,  jb,:) = 0.0_wp
       field% albnirdir_tile (jcs:jce,  jb,:) = 0.0_wp
       field% albvisdif_tile (jcs:jce,  jb,:) = 0.0_wp
       field% albnirdif_tile (jcs:jce,  jb,:) = 0.0_wp
       field% albedo         (jcs:jce,  jb  ) = 0.0_wp
       field% albedo_tile    (jcs:jce,  jb,:) = 0.0_wp
       field% ts             (jcs:jce,  jb  ) = 0.0_wp
       field% ts_rad         (jcs:jce,  jb  ) = 0.0_wp
       field% lwflxsfc_tile  (jcs:jce,  jb,:) = 0.0_wp
       field% swflxsfc_tile  (jcs:jce,  jb,:) = 0.0_wp
       !
       field% Tsurf          (jcs:jce,:,jb  ) = 0.0_wp
       field% T1             (jcs:jce,:,jb  ) = 0.0_wp
       field% T2             (jcs:jce,:,jb  ) = 0.0_wp
       field% Qtop           (jcs:jce,:,jb  ) = 0.0_wp
       field% Qbot           (jcs:jce,:,jb  ) = 0.0_wp
       field% albvisdir_ice  (jcs:jce,:,jb  ) = 0.0_wp
       field% albnirdir_ice  (jcs:jce,:,jb  ) = 0.0_wp
       field% albvisdif_ice  (jcs:jce,:,jb  ) = 0.0_wp
       field% albnirdif_ice  (jcs:jce,:,jb  ) = 0.0_wp
       !
       tend% ta_sfc          (jcs:jce,  jb  ) = 0.0_wp
       !
       ! vdiff_up
       field% tke            (jcs:jce,:,jb  ) = 0.0_wp
       field% thvvar         (jcs:jce,:,jb  ) = 0.0_wp
       field%   xvar         (jcs:jce,:,jb  ) = 0.0_wp
       field% z0m            (jcs:jce,  jb  ) = 0.0_wp
       field% z0m_tile       (jcs:jce,  jb,:) = 0.0_wp
       field% kedisp         (jcs:jce,  jb  ) = 0.0_wp
       field% sh_vdiff       (jcs:jce,  jb  ) = 0.0_wp
       field% qv_vdiff       (jcs:jce,  jb  ) = 0.0_wp
       !
       tend%   ua_vdf        (jcs:jce,:,jb  ) = 0.0_wp
       tend%   va_vdf        (jcs:jce,:,jb  ) = 0.0_wp
       tend%   ta_vdf        (jcs:jce,:,jb  ) = 0.0_wp
       tend% qtrc_vdf        (jcs:jce,:,jb,:) = 0.0_wp
       !
    END IF

  END SUBROUTINE echam_vdiff_down_surf_up
  !-------------------------------------------------------------------

END MODULE mo_interface_echam_vdiff_surface
