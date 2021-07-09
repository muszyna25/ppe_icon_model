#include "fsel.inc"
!>
!! @brief Subroutines for computing turbulent exchange coefficients of
!! 3D Smagorinsky turbulent scheme.
!!
!! @author Junhong Lee, MPI-M
!!
!! @par Revision History
!! Code originates from ICON-LEM, mo_sgs_turbulence.f90
!! Re-organized by Junhong Lee (2020-02).
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_sma_turbulence_diag

  USE mo_kind              ,ONLY: wp
  USE mo_convect_tables    ,ONLY: compute_qsat
!!$  USE mo_echam_convect_tables, ONLY: prepare_ua_index_spline, lookup_ua_spline
  USE mo_echam_vdf_config  ,ONLY: echam_vdf_config
  USE mo_echam_vdiff_params,ONLY: ckap
  USE mo_physical_constants,ONLY: grav, rd, cpd, cpv, rd_o_cpd,             &
    &                             vtmpc1, p0ref, rgrav
  USE mo_model_domain      ,ONLY: t_patch
  USE mo_nonhydro_types    ,ONLY: t_nh_prog, t_nh_metrics
  USE mo_intp_data_strc    ,ONLY: t_int_state
  USE mo_nonhydro_state    ,ONLY: p_nh_state
  USE mo_dynamics_config   ,ONLY: nnew
  USE mo_intp_data_strc    ,ONLY: p_int_state
  USE mo_fortran_tools     ,ONLY: init
  USE mo_impl_constants    ,ONLY: min_rlcell, min_rledge_int, min_rlcell_int, &
    &                             min_rlvert_int
  USE mo_parallel_config   ,ONLY: p_test_run
  USE mo_loopindices       ,ONLY: get_indices_e, get_indices_c
  USE mo_les_utilities     ,ONLY: brunt_vaisala_freq, vert_intp_full2half_cell_3d
  USE mo_intp              ,ONLY: cells2verts_scalar, cells2edges_scalar
  USE mo_sync              ,ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array, &
    &                             sync_patch_array_mult
  USE mo_intp_rbf          ,ONLY: rbf_vec_interpol_vertex, rbf_vec_interpol_edge
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_nh_testcases_nml  ,ONLY: is_dry_cbl
  USE mo_satad             ,ONLY: spec_humi, sat_pres_water
  USE mo_math_constants    ,ONLY: pi_2, ln2
  USE mo_math_utilities    ,ONLY: tdma_solver
  USE mo_intp_rbf          ,ONLY: rbf_vec_interpol_cell
  USE mo_nh_testcases_nml  ,ONLY: isrfc_type, ufric

  IMPLICIT NONE
  PRIVATE
  REAL(wp),         PARAMETER :: z_1by3  = 1._wp/3._wp

  !Parameter for vertical scheme type
  INTEGER, PARAMETER :: iexplicit = 1
  INTEGER, PARAMETER :: iimplicit = 2

  CHARACTER(len=*), PARAMETER :: inmodule = 'mo_sgs_turbulence:'

  PUBLIC :: atm_exchange_coeff3d, diffuse_hori_velocity, diffuse_vert_velocity, &
          & diffuse_scalar

  !Parameters for surface layer parameterizations: From Zeng_etal 1997 J. Clim
  REAL(wp), PARAMETER :: bsm = 5.0_wp  !Businger Stable Momentum
  REAL(wp), PARAMETER :: bum = 16._wp  !Businger Untable Momentum
  REAL(wp), PARAMETER :: bsh = 5.0_wp  !Businger Stable Heat
  REAL(wp), PARAMETER :: buh = 16._wp  !Businger Untable Heat

CONTAINS
  !>
  !! Compute various thermodynamic variables for all (full) vertical levels;
  !! Diagnose PBL extension;
  !! Diagnose wind shear, buoyancy, Ri-number, mixing length, then compute
  !! the turbulent exchange coefficients of momentum, dry static energy,
  !! tracers, TTE, variance of virtual optential temperation at half levels
  !! [1+1/2, klev-1/2].
  !!
  !! Note that
  !! - for all coeffcient arrays, vertical index k in this subroutine
  !!   correspond to interface (half level) k+1/2;
  !! - the exchange coefficient at model top (level 1/2) is zero, thus does
  !!   not need computing;
  !! Therefore a number of variables are defined on what shall be referred to
  !! as "mid-levels", which are defined on levels with indices xmid(k) corre-
  !! sponding to xh(k+1), where x is some quantity at half-level k+1/2.
  !!
  !! @par Revision History
  !! Separated from vdiff.f90 of ECHAM6 and re-organized by Hui Wan (2010-09).
  !!  updated to echam-6.3.01 by Monika Esch (2014-11)
  !!
  SUBROUTINE atm_exchange_coeff3d( jg, kbdim, nblks_c, nblks_v, nblks_e,  &! in
                               & klev, klevm1, klevp1,                    &! in
                               & ksfc_type, idx_lnd,                      &! in
                               & p_patch,                                 &! in
                               & pz0m, ptsfc, pfrc,                       &! in
                               & ppsfc,                                   &! in
                               & pghf, pghh,                              &! in
                               & pum1, pvm1, pwp1,                        &! in
                               & ptm1, ptvm1,                             &! in
                               & pqm1, pxm1,                              &! in
                               & rho,                                     &! in
                               & papm1, paphm1,                           &! in
                               & pri_tile,                                &! in
                               & pthvsig,                                 &! in
                               & pcfm_tile,                               &! out
                               & pcfh_tile,                               &! out
                               & pqsat_tile, pcpt_tile,                   &! out
                               & pcptgz,                                  &! out
                               & pzthvvar, ptottevn,                      &! out
                               & pmixlen,                                 &! out
                               & pcfm, pcfh, pcfv, pcftotte, pcfthv,      &! out
                               & km_c, km_iv, km_ie, km_ic, kh_ic,        &! out
                               & pprfac,                                  &! out
                               & u_vert, v_vert, div_c,                   &! out
                               & rho_ic, w_vert, w_ie,                    &! out
                               & vn,                                      &! out
                               & pch_tile,                                &! out
                               & pbn_tile, pbhn_tile, pbm_tile, pbh_tile, &! out
                               & paz0lh,                                  &! in, optional
                               & pcsat, pcair                             &! in, optional
                               & )

    ! Arguments

    INTEGER, INTENT(IN) :: jg, nblks_c, nblks_v, nblks_e
    INTEGER, INTENT(IN) :: kbdim
    INTEGER :: nproma
    INTEGER, INTENT(IN) :: klev, klevm1, klevp1
    INTEGER :: nlev, nlevm1, nlevp1
    REAL(wp),INTENT(IN) :: pghf(kbdim,klev,nblks_c)
    REAL(wp),INTENT(IN) :: pghh(kbdim,klevp1,nblks_c)
    REAL(wp),INTENT(IN) :: pxm1(kbdim,klev,nblks_c)
    REAL(wp),INTENT(IN) :: ptvm1(kbdim,klev,nblks_c)
    REAL(wp),INTENT(IN) :: pqm1(kbdim,klev,nblks_c)
    REAL(wp),INTENT(IN) :: pwp1(kbdim,klevp1,nblks_c)
    REAL(wp),INTENT(IN) :: ptm1(kbdim,klev,nblks_c), rho(kbdim,klev,nblks_c)
    REAL(wp),INTENT(IN) :: papm1(kbdim,klev,nblks_c),  paphm1(kbdim,klevp1,nblks_c)

    REAL(wp),INTENT(INOUT) :: pum1(kbdim,klev,nblks_c),  pvm1(kbdim,klev,nblks_c)

    REAL(wp),INTENT(OUT) :: ptottevn(kbdim,klev,nblks_c) !< TTE at intermediate time step
    REAL(wp),INTENT(OUT) :: pcftotte(kbdim,klev,nblks_c) !< exchange coeff. for TTE
    REAL(wp),INTENT(OUT) :: pcfthv  (kbdim,klev,nblks_c) !< exchange coeff. for var. of theta_v
    REAL(wp),INTENT(OUT) :: pcfm    (kbdim,klev,nblks_c) !< exchange coeff. for u, v
    REAL(wp),INTENT(OUT) :: pcfh    (kbdim,klev,nblks_c) !< exchange coeff. for cptgz and tracers
    REAL(wp),INTENT(OUT) :: pcfv    (kbdim,klev,nblks_c) !< exchange coeff. for variance of qx
    REAL(wp),INTENT(OUT) :: pzthvvar(kbdim,klev,nblks_c) !< variance of theta_v at interm. step
    REAL(wp),INTENT(OUT) :: pcptgz  (kbdim,klev,nblks_c)   !< dry static energy
    REAL(wp),INTENT(OUT) :: pprfac  (kbdim,klev,nblks_c) !< prefactor for the exchange coeff.
    REAL(wp),INTENT(OUT) :: pmixlen (kbdim,klev,nblks_c) !< prefactor for the exchange coeff.
    REAL(wp),INTENT(OUT) :: pthvsig (kbdim,nblks_c)

    ! Local variables
    ! - Variables defined at full levels

    REAL(wp) :: ztheta (kbdim,klev,nblks_c)  !< potential temperature

    ! - Variables defined at mid-levels

    REAL(wp) :: zdgmid(kbdim,klevm1,nblks_c) !< geopotential height difference between two full levels
    REAL(wp) :: ztvmid(kbdim,klevm1,nblks_c)

    TYPE(t_patch)   ,TARGET ,INTENT(inout)   :: p_patch
    TYPE(t_nh_metrics) ,POINTER :: p_nh_metrics
    TYPE(t_nh_prog)    ,POINTER :: p_nh_prog     !<the prognostic variables
    TYPE(t_int_state)  ,POINTER :: p_int         !< interpolation state

    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end
    INTEGER :: jc, jb, je

  !Variables for the module
    REAL(wp), INTENT(OUT), DIMENSION(kbdim,klev,nblks_c)   :: km_c
    REAL(wp), INTENT(OUT), DIMENSION(kbdim,klevp1,nblks_v) :: km_iv
    REAL(wp), INTENT(OUT), DIMENSION(kbdim,klevp1,nblks_e) :: km_ie
    REAL(wp), INTENT(OUT), DIMENSION(kbdim,klevp1,nblks_c) :: kh_ic, km_ic
    REAL(wp), INTENT(OUT), DIMENSION(kbdim,klev,nblks_v)   :: u_vert, v_vert
    REAL(wp), INTENT(OUT), DIMENSION(kbdim,klev,nblks_c)   :: div_c
    REAL(wp), INTENT(OUT), DIMENSION(kbdim,klevp1,nblks_v) :: w_vert
    REAL(wp), INTENT(OUT), DIMENSION(kbdim,klevp1,nblks_e) :: w_ie
    REAL(wp), INTENT(OUT), DIMENSION(kbdim,klevp1,nblks_c) :: rho_ic
    REAL(wp), INTENT(OUT), DIMENSION(kbdim,klev,nblks_e)   :: vn !< normal wind vector

    REAL(wp), DIMENSION(kbdim,klev,nblks_c)   :: theta, theta_v
    REAL(wp), DIMENSION(kbdim,klevp1,nblks_c) :: bruvais

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: vn_ie, vt_ie, mech_prod, &
                                                 shear, div_of_stress

    INTEGER, DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk, ieidx, ieblk

    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4
    REAL(wp) :: vt_vert1, vt_vert2, vt_vert3, vt_vert4
    REAL(wp) :: w_full_c1, w_full_c2, w_full_v1, w_full_v2
    REAL(wp) :: D_11, D_12, D_13, D_22, D_23, D_33

    ! - surface variables
    INTEGER, INTENT(IN) :: ksfc_type, idx_lnd

    REAL(wp),INTENT(IN) :: pz0m  (kbdim,nblks_c,ksfc_type) !< aerodynamic roughness length
    REAL(wp),INTENT(IN) :: ptsfc (kbdim,nblks_c,ksfc_type) !< temp. at surface
    REAL(wp),INTENT(IN) :: pfrc  (kbdim,nblks_c,ksfc_type) !< fraction of the grid box occupied
    REAL(wp),INTENT(IN) :: ppsfc (kbdim,nblks_c)  !< surface pressure

    ! optional arguments for use with jsbach
    REAL(wp),OPTIONAL,INTENT(IN) :: paz0lh (kbdim,nblks_c)  !< roughness length for heat over land
    REAL(wp),OPTIONAL,INTENT(IN) :: pcsat  (kbdim,nblks_c)  !< area fraction with wet land surface
    REAL(wp),OPTIONAL,INTENT(IN) :: pcair  (kbdim,nblks_c)  !< area fraction with wet land surface (air)

    REAL(wp),INTENT(OUT) :: pqsat_tile (kbdim,nblks_c,ksfc_type)!< saturation specific humidity
    REAL(wp),INTENT(OUT) :: pcpt_tile (kbdim,nblks_c,ksfc_type) !< dry static energy
    REAL(wp),INTENT(OUT) :: pcfm_tile (kbdim,nblks_c,ksfc_type) !< exchange coeff. of momentum,
                                                                !< for each type of surface
    REAL(wp),INTENT(OUT) :: pcfh_tile (kbdim,nblks_c,ksfc_type) !< exchange coeff. of heat and
                                                                !<  vapor for each surface type
    REAL(wp),INTENT(OUT) :: pbn_tile  (kbdim,nblks_c,ksfc_type) !< for diagnostics
    REAL(wp),INTENT(OUT) :: pbhn_tile (kbdim,nblks_c,ksfc_type) !< for diagnostics
    REAL(wp),INTENT(OUT) :: pbm_tile  (kbdim,nblks_c,ksfc_type) !< for diagnostics
    REAL(wp),INTENT(OUT) :: pbh_tile  (kbdim,nblks_c,ksfc_type) !< for diagnostics
    REAL(wp),INTENT(OUT) :: pch_tile  (kbdim,nblks_c,ksfc_type) !< for TTE boundary condition
    REAL(wp),INTENT(OUT) :: pri_tile  (kbdim,nblks_c,ksfc_type) !< Richardson number for diagnostics

    REAL(wp) :: zqts
    REAL(wp) :: dummy
    REAL(wp) :: zvn1, zvn2
    INTEGER  :: loidx (kbdim,nblks_c,ksfc_type) !< counter for masks
    INTEGER  :: is    (ksfc_type)               !< counter for masks
    INTEGER  :: jsfc, jls, js
    INTEGER  :: jcn,jbn                         !< jc and jb of neighbor cells sharing an edge je
    REAL(wp),parameter :: zcons17 = 1._wp / ckap**2

    ! - surface variables (mo_surface_les.f90)
    REAL(wp) :: umfl_tile(kbdim,nblks_c,ksfc_type)
    REAL(wp) :: vmfl_tile(kbdim,nblks_c,ksfc_type)

    INTEGER  :: itr
    REAL(wp) :: zrough, theta_sfc, qv_s, rhos, mwind, z_mc, RIB, tcn_mom, tcn_heat, &
                shfl, lhfl, bflx1, ustar, obukhov_length, inv_bus_mom
    REAL(wp) :: tch(kbdim,nblks_c,ksfc_type)
    REAL(wp) :: tcm(kbdim,nblks_c,ksfc_type)
    REAL(wp) :: zthetavmid (kbdim,nblks_c,ksfc_type)



    ! - 1D variables and scalars

    INTEGER  :: jk, jl
    REAL(wp) :: zrdp
    REAL(wp) :: zsdep1
    REAL(wp) :: zsdep2

    ! to prevent floating-point arithmetic inconsistencies later in
    ! the interpolation to u 10m and 2m T/T_d: has been 0.01 before
    ! (Cray FP instead of IEEE 754 FP format)
    REAL(wp) :: zepsec = 0.028_wp

    ! Shortcuts to components of echam_vdf_config
    !
    REAL(wp) :: fsl
    !

    p_nh_metrics => p_nh_state(jg)%metrics
    p_nh_prog    => p_nh_state(jg)%prog(nnew(jg))
    p_int        => p_int_state(jg)

    nproma = kbdim
    nlev = klev; nlevm1 = klevm1; nlevp1 = klevp1

    fsl      = echam_vdf_config(jg)% fsl

    rl_start   = 1
    rl_end     = min_rlcell
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        ptottevn(jc,:,jb) = 10._wp  ! for vdiff_up
        pzthvvar(jc,:,jb) = 1._wp   ! for vdiff_up

        pcftotte(jc,:,jb) = 10._wp  ! K TEE
        pcfthv(jc,:,jb)   = 0._wp   ! K theta_v
        pcfv(jc,:,jb)     = 0._wp   ! K qx
        pmixlen(jc,:,jb)  = 0._wp   ! K qx

        pcfm_tile(jc,jb,:) = 0._wp
        pcfh_tile(jc,jb,:) = 0._wp

        pqsat_tile(jc,jb,:)= 0._wp
        pri_tile(jc,jb,:)  = 0._wp
        pch_tile(jc,jb,:)  = 0._wp
        pcpt_tile(jc,jb,:) = 0._wp

        pthvsig(jc,jb)     = 0._wp
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   CALL sync_patch_array(SYNC_C, p_patch, pum1)
   CALL sync_patch_array(SYNC_C, p_patch, pvm1)

    rl_start   = grf_bdywidth_e+1
    rl_end     = min_rledge_int
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,jcn,jbn,zvn1,zvn2)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
          jcn  =   p_patch%edges%cell_idx(je,jb,1)
          jbn  =   p_patch%edges%cell_blk(je,jb,1)
          zvn1 =   pum1(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,1)%v1 &
            &    + pvm1(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,1)%v2
          !
          jcn  =   p_patch%edges%cell_idx(je,jb,2)
          jbn  =   p_patch%edges%cell_blk(je,jb,2)
          zvn2 =   pum1(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,2)%v1 &
            &    + pvm1(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,2)%v2
          !
          vn(je,jk,jb) = p_int%c_lin_e(je,1,jb)*zvn1 &
            &          + p_int%c_lin_e(je,2,jb)*zvn2
        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

   CALL sync_patch_array(SYNC_E, p_patch, vn)

!#########################################################################
!## variables for TTE scheme and JSBACH LSM
!#########################################################################

    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jl,jls,js,i_startidx,i_endidx,is,jsfc,zqts,zrough, theta_sfc,   &
!$OMP            qv_s, mwind, rhos, z_mc, RIB, tcn_mom, tcn_heat, itr, shfl, lhfl,&
!$OMP            bflx1, ustar, obukhov_length, inv_bus_mom)

  DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)


    DO jl = i_startidx, i_endidx
      DO jk = 1, klevm1
        zdgmid(jl,jk,jb) = pghf(jl,jk,jb) - pghf(jl,jk+1,jb)

        ! interpolation coefficients
        zrdp   = 1._wp/(papm1(jl,jk,jb) - papm1(jl,jk+1,jb))
        zsdep1 = (papm1(jl,jk,jb)  - paphm1(jl,jk+1,jb))*zrdp
        zsdep2 = (paphm1(jl,jk+1,jb)- papm1(jl,jk+1,jb))*zrdp

        ztvmid(jl,jk,jb)=zsdep1*ptvm1(jl,jk,jb)+zsdep2*ptvm1(jl,jk+1,jb)

        ! Virtual dry static energy
        pcptgz(jl,jk,jb) = pghf(jl,jk,jb)*grav+ptm1(jl,jk,jb)*(cpd+(cpv-cpd)*pqm1(jl,jk,jb))

        ! Potential temperature
        ztheta(jl,jk,jb) = ptm1(jl,jk,jb)*(p0ref/papm1(jl,jk,jb))**rd_o_cpd

        ! Air density at mid levels, p/(Tv*R)/dz = air density/dz, and the prefactor
        ! that will be multiplied later to the exchange coeffcients to build a linear
        ! algebraic equation set.
        pprfac(jl,jk,jb) = paphm1(jl,jk+1,jb)/(ztvmid(jl,jk,jb)*rd*zdgmid(jl,jk,jb))

      END DO

      ! dry static energy pcpt_tile
      !
      pcptgz(jl,klev,jb) = pghf(jl,klev,jb)*grav+ptm1(jl,klev,jb)*(cpd+(cpv-cpd)*pqm1(jl,klev,jb))

      ! Potential temperature
      ztheta(jl,klev,jb) = ptm1(jl,klev,jb)*(p0ref/papm1(jl,klev,jb))**rd_o_cpd

      ! The prefactor (= air density) that will be multiplied to the exchange
      ! coefficients when building the linear algebraic equations. The density here is
      ! computed using air temperature of the lowest model level at time step n-1.
      pprfac(jl,klev,jb) = ppsfc(jl,jb)                                       &
                           &  /( rd*ptm1(jl,klev,jb)                             &
                           &  *(1._wp+vtmpc1*pqm1(jl,klev,jb)-pxm1(jl,klev,jb)) )
    END DO


    pcfm(:,klev,jb) = 0._wp
    pcfh(:,klev,jb) = 0._wp

!private is, loidx,

    DO jsfc = 1, ksfc_type
     ! check for masks
     !
      is(jsfc) = 0
      DO jl = i_startidx, i_endidx
        IF(pfrc(jl,jb,jsfc).GT.0.0_wp) THEN
          is(jsfc) = is(jsfc) + 1
          loidx(is(jsfc),jb,jsfc) = jl
        ENDIF
      ENDDO

      CALL compute_qsat( kbdim, is(jsfc), loidx(:,jb,jsfc), ppsfc(:,jb), ptsfc(:,jb,jsfc), pqsat_tile(:,jb,jsfc) )

     ! loop over mask only
     !
      DO jls = 1, is(jsfc)
        js=loidx(jls,jb,jsfc)

        ! dry static energy pcpt_tile
        !
        IF(jsfc == idx_lnd ) THEN
          zqts = pcsat(js,jb) * pqsat_tile(js,jb,jsfc) + (1._wp - pcair(js,jb))*pqm1(js,klev,jb) ! q_total at land surface
        ELSE
          zqts = pqsat_tile(js,jb,jsfc)                                              ! q_total at non-land surface
        END IF

        pcpt_tile(js,jb,jsfc) = ptsfc(js,jb,jsfc) * (cpd + (cpv - cpd) * zqts)

        !Surface roughness length
        IF ( pz0m(js,jb,jsfc) .GT. 0.0_wp ) THEN
          zrough = pz0m(js,jb,jsfc)
        ELSE
          zrough = 1.E-3_wp
        END IF

        !Get surface pot. temperature and humidity
        theta_sfc = ptsfc(js,jb,jsfc) / EXP( rd_o_cpd*LOG(ppsfc(js,jb)/p0ref) )
        qv_s    = spec_humi(sat_pres_water(ptsfc(js,jb,jsfc)),ppsfc(js,jb))


        !rho at surface: no qc at suface
        rhos   =  ppsfc(js,jb)/( rd * &
                  ptsfc(js,jb,jsfc)*(1._wp+vtmpc1*qv_s) )

        mwind = MAX( echam_vdf_config(jg)%min_sfc_wind, SQRT(pum1(js,klev,jb)**2+pvm1(js,klev,jb)**2) )

        !Z height to be used as a reference height in surface layer
        z_mc   = p_nh_metrics%z_mc(js,klev,jb) - p_nh_metrics%z_ifc(js,klevp1,jb)

        !First guess for tch and tcm using bulk approach
        RIB = grav * (ztheta(js,klev,jb)-theta_sfc) * (z_mc-zrough) / (theta_sfc * mwind**2)
        tcn_mom = (ckap/LOG(z_mc/zrough))**2
        tcm(js,jb,jsfc) = tcn_mom * stability_function_mom(RIB,z_mc/zrough,tcn_mom)

        tcn_heat        = ckap**2/(LOG(z_mc/zrough)*LOG(z_mc/zrough))
        tch(js,jb,jsfc) = tcn_heat * stability_function_heat(RIB,z_mc/zrough,tcn_heat)

        !now iterate
        DO itr = 1 , 5
           shfl = tch(js,jb,jsfc)*mwind*(theta_sfc-ztheta(js,klev,jb))
           lhfl = tch(js,jb,jsfc)*mwind*(qv_s-pqm1(js,klev,jb))
           bflx1= shfl + vtmpc1 * theta_sfc * lhfl
           ustar= SQRT(tcm(js,jb,jsfc))*mwind

           obukhov_length = -ustar**3 * theta_sfc * rgrav / (ckap * bflx1)

           inv_bus_mom = 1._wp / businger_mom(zrough,z_mc,obukhov_length)
           tch(js,jb,jsfc) = inv_bus_mom / businger_heat(zrough,z_mc,obukhov_length)
           tcm(js,jb,jsfc) = inv_bus_mom * inv_bus_mom
        END DO

        IF ( isrfc_type == 1 ) THEN
          umfl_tile(js,jb,jsfc) = ufric**2 * rhos * pum1(js,klev,jb) / mwind
          vmfl_tile(js,jb,jsfc) = ufric**2 * rhos * pvm1(js,klev,jb) / mwind
        ELSE
          umfl_tile(js,jb,jsfc) = rhos*tcm(js,jb,jsfc)*mwind*pum1(js,klev,jb)
          vmfl_tile(js,jb,jsfc) = rhos*tcm(js,jb,jsfc)*mwind*pvm1(js,klev,jb)
        END IF
        pcfm_tile(js,jb,jsfc) = tcm(js,jb,jsfc)*mwind
        pcfh_tile(js,jb,jsfc) = tch(js,jb,jsfc)*mwind
        pch_tile (js,jb,jsfc) = tch(js,jb,jsfc)

        pcfm(js,klev,jb) = pcfm(js,klev,jb) + pfrc(js,jb,jsfc)*pcfm_tile(js,jb,jsfc)
        pcfh(js,klev,jb) = pcfh(js,klev,jb) + pfrc(js,jb,jsfc)*pcfh_tile(js,jb,jsfc)


        pbn_tile(js,jb,jsfc) = ckap / MAX( zepsec, sqrt(tcn_mom) )
        pbhn_tile(js,jb,jsfc)= ckap / MAX( zepsec, sqrt(tcn_heat) )
        pbm_tile(js,jb,jsfc) = MAX( zepsec, sqrt(pcfm_tile(js,jb,jsfc) * tch(js,jb,jsfc)*zcons17/ (tcn_mom*mwind)) )
        pbh_tile(js,jb,jsfc) = MAX( zepsec, tch(js,jb,jsfc)/pbm_tile(js,jb,jsfc)*zcons17)
        pbm_tile(js,jb,jsfc) = 1._wp / pbm_tile(js,jb,jsfc)
        pbh_tile(js,jb,jsfc) = 1._wp / pbh_tile(js,jb,jsfc)

        zthetavmid(js,jb,jsfc) = fsl*ptvm1(js,nlev,jb)*(p0ref/papm1(js,nlev,jb))**rd_o_cpd + &
                               & (1._wp-fsl) * theta_sfc * (1._wp+vtmpc1*zqts)

        pri_tile(js,jb,jsfc) = pghf(js,klev,jb) * grav *                                   &
                             & ( ptvm1(js,nlev,jb)*(p0ref/papm1(js,nlev,jb))**rd_o_cpd -   &
                             &   theta_sfc * (1._wp+vtmpc1*zqts)                       ) / &
                             & ( zthetavmid(js,jb,jsfc) * mwind )
      END DO !jls
    END DO !jsfc


  END DO !jb
!$OMP END DO
!$OMP END PARALLEL




!#########################################################################
!## initialize
!#########################################################################

!$OMP PARALLEL
    CALL init(km_iv(:,:,:))
    CALL init(km_c(:,:,:))
    CALL init(km_ie(:,:,:))

    IF(p_test_run)THEN
      CALL init(u_vert(:,:,:))
      CALL init(v_vert(:,:,:))
      CALL init(w_vert(:,:,:))
    END IF
!$OMP END PARALLEL


!#########################################################################
!## Convert temperature to potential temperature: all routines within
!## use theta.
!#########################################################################


    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        theta(jc,1:nlev,jb) = ptm1(jc,1:nlev,jb)*(p0ref/papm1(jc,1:nlev,jb))**rd_o_cpd
        theta_v(jc,1:nlev,jb) = ptvm1(jc,1:nlev,jb)*(p0ref/papm1(jc,1:nlev,jb))**rd_o_cpd
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    !Get rho at interfaces to be used later
    CALL vert_intp_full2half_cell_3d(p_patch, p_nh_metrics, rho, rho_ic, &
                                     2, min_rlcell_int-2)

    CALL brunt_vaisala_freq(p_patch, p_nh_metrics, theta_v, bruvais)



!#########################################################################
!## Smagorinsky_model
  !!------------------------------------------------------------------------
  !! Computes the sgs viscosity and diffusivity using Smagorinsky model
  !! \tau_ij = KD_ij where D_ij = du_i/dx_j + du_j/dx_i
  !!
  !! and  K = cs * \Delta**2 * D / sqrt(2), where D = sqrt(D_ijD_ij)
  !!
  !! and D**2 = D_11**2 + D_22**2 + D_33**2 + 2D_12**2 + 2D_13**2 + 2D_23**2
  !!
  !! where, D_11 = 2 * du_1/dx_1
  !!        D_22 = 2 * d_u2/dx_2
  !!        D_33 = 2 * d_u3/dx_3
  !!        D_12 = du_1/dx_2 + du_2/dx_1
  !!        D_13 = du_1/dx_3 + du_3/dx_1
  !!        D_23 = du_2/dx_3 + du_3/dx_2
  !! For triangles: 1=normal, 2=tangential, and 3 = z directions
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-20)
  !! Modify by Junhong Lee, MPI-M (2020-03-03)

    ALLOCATE( vn_ie(nproma,p_patch%nlevp1,p_patch%nblks_e)      &
             ,vt_ie(nproma,p_patch%nlevp1,p_patch%nblks_e)      &
             ,shear(nproma,p_patch%nlev,p_patch%nblks_e)        &
             ,div_of_stress(nproma,p_patch%nlev,p_patch%nblks_e)&
             ,mech_prod(nproma,nlevp1,p_patch%nblks_c)          &
            )

    IF(p_test_run)THEN
      kh_ic(:,:,:) = 0._wp
      km_ic(:,:,:) = 0._wp
    END IF

    !--------------------------------------------------------------------------
    !1) Interpolate velocities at desired locations- mostly around the quadrilateral
    !
    !<It assumes that prog values are all synced at this stage while diag values might not>
    !--------------------------------------------------------------------------


    CALL cells2verts_scalar(pwp1, p_patch, p_int%cells_aw_verts, w_vert,                   &
                            opt_rlend=min_rlvert_int)
    CALL cells2edges_scalar(pwp1, p_patch, p_int%c_lin_e, w_ie, opt_rlend=min_rledge_int-2)

    ! RBF reconstruction of velocity at vertices: include halos
    CALL rbf_vec_interpol_vertex( vn, p_patch, p_int, &
                                  u_vert, v_vert, opt_rlend=min_rlvert_int )

    !sync them
    CALL sync_patch_array_mult(SYNC_V, p_patch, 3, w_vert, u_vert, v_vert)

    !Get vn at interfaces and then get vt at interfaces
    !Boundary values are extrapolated like dynamics although
    !they are not required in current implementation

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)
    rl_start   = 2
    rl_end     = min_rledge_int-3
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO je = i_startidx, i_endidx
#endif
          vn_ie(je,jk,jb) = p_nh_metrics%wgtfac_e(je,jk,jb) * vn(je,jk,jb) +            &
                            ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * vn(je,jk-1,jb)
        END DO
      END DO
      DO je = i_startidx, i_endidx
        vn_ie(je,1,jb)      = p_nh_metrics%wgtfacq1_e(je,1,jb) * vn(je,1,jb) +          &
                              p_nh_metrics%wgtfacq1_e(je,2,jb) * vn(je,2,jb) +          &
                              p_nh_metrics%wgtfacq1_e(je,3,jb) * vn(je,3,jb)

        vn_ie(je,nlevp1,jb) = p_nh_metrics%wgtfacq_e(je,1,jb) * vn(je,nlev,jb)   +      &
                              p_nh_metrics%wgtfacq_e(je,2,jb) * vn(je,nlev-1,jb) +      &
                              p_nh_metrics%wgtfacq_e(je,3,jb) * vn(je,nlev-2,jb)
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    CALL rbf_vec_interpol_edge(vn_ie, p_patch, p_int, vt_ie, opt_rlstart=3, &
                               opt_rlend=min_rledge_int-2)

    !--------------------------------------------------------------------------
    !2) Compute horizontal strain rate tensor at full levels
    !--------------------------------------------------------------------------
    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk


    rl_start   = 4
    rl_end     = min_rledge_int-2
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,vn_vert1,vn_vert2,vn_vert3,vn_vert4,  &
!$OMP            vt_vert1,vt_vert2,vt_vert3,vt_vert4,w_full_c1,w_full_c2,w_full_v1,  &
!$OMP            w_full_v2,D_11,D_12,D_13,D_22,D_23,D_33)
    DO jb = i_startblk,i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif

          vn_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,1)%v1 +   &
                     v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,1)%v2

          vn_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,2)%v1 +   &
                     v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,2)%v2

          vn_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,3)%v1 +   &
                     v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,3)%v2

          vn_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,4)%v1 +   &
                     v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,4)%v2

          vt_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,1)%v1   +   &
                     v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,1)%v2

          vt_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,2)%v1   +   &
                     v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,2)%v2

          vt_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,3)%v1   +   &
                     v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,3)%v2

          vt_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,4)%v1   +   &
                     v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,4)%v2

          ! W at full levels
          w_full_c1  = 0.5_wp *                                              &
                       ( pwp1(iecidx(je,jb,1),jk,iecblk(je,jb,1)) +   &
                         pwp1(iecidx(je,jb,1),jk+1,iecblk(je,jb,1)) )

          w_full_c2  = 0.5_wp *                                              &
                       ( pwp1(iecidx(je,jb,2),jk,iecblk(je,jb,2)) +   &
                         pwp1(iecidx(je,jb,2),jk+1,iecblk(je,jb,2)) )

          ! W at full levels vertices from w at vertices at interface levels
          w_full_v1  = 0.5_wp *                                              &
                       ( w_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) +          &
                         w_vert(ividx(je,jb,1),jk+1,ivblk(je,jb,1)) )

          w_full_v2  = 0.5_wp *                                              &
                       ( w_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) +          &
                         w_vert(ividx(je,jb,2),jk+1,ivblk(je,jb,2)) )


          ! Strain rates at edge center
          D_11       = 2._wp * ( vn_vert4 - vn_vert3 ) *                     &
                       p_patch%edges%inv_vert_vert_length(je,jb)

          D_12       = p_patch%edges%tangent_orientation(je,jb) *            &
                       ( vn_vert2 - vn_vert1 ) *                             &
                       p_patch%edges%inv_primal_edge_length(je,jb) +         &
                       ( vt_vert4-vt_vert3 ) *                               &
                       p_patch%edges%inv_vert_vert_length(je,jb)

          D_13       = ( vn_ie(je,jk,jb) - vn_ie(je,jk+1,jb) ) *             &
                       p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)  +           &
                       ( w_full_c2 - w_full_c1 ) *                           &
                       p_patch%edges%inv_dual_edge_length(je,jb)

          D_22       = 2._wp * ( vt_vert2-vt_vert1 ) *                       &
                       p_patch%edges%tangent_orientation(je,jb) *            &
                       p_patch%edges%inv_primal_edge_length(je,jb)

          D_23       = ( vt_ie(je,jk,jb) - vt_ie(je,jk+1,jb) ) *             &
                       p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)  +           &
                       p_patch%edges%tangent_orientation(je,jb) *            &
                       ( w_full_v2 - w_full_v1 ) *                           &
                       p_patch%edges%inv_primal_edge_length(je,jb)

          D_33       = 2._wp * ( w_ie(je,jk,jb) - w_ie(je,jk+1,jb) ) *       &
                       p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)

          ! Mechanical prod is half of this value divided by km
          shear(je,jk,jb) = D_11**2 + D_22**2 + D_33**2 +                    &
                            2._wp * ( D_12**2 + D_13**2 + D_23**2 )

          ! calculate divergence to get the deviatoric part of stress tensor in
          ! diffusion: D_11 - 1/3 * (D_11 + D_22 + D_33)
          div_of_stress(je,jk,jb) = 0.5_wp * ( D_11 + D_22 + D_33 )
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL


    !Interpolate mech production term from mid level edge to interface level cell
    !except top and bottom boundaries
    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int-1
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
#endif
          div_c(jc,jk,jb) =                                                                       &
                  ( div_of_stress(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb) +  &
                    div_of_stress(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb) +  &
                    div_of_stress(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) )
        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL


    ! Interpolate mech. production term from mid level edge to interface level cell
    ! except top and bottom boundaries
    rl_start   = 3
    rl_end     = min_rlcell_int-1
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
#endif

          mech_prod(jc,jk,jb) = p_nh_metrics%wgtfac_c(jc,jk,jb) * (                      &
                      shear(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb)   +      &
                      shear(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb)   +      &
                      shear(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) ) +      &
                      ( 1._wp - p_nh_metrics%wgtfac_c(jc,jk,jb) ) * (                             &
                      shear(ieidx(jc,jb,1),jk-1,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb) +      &
                      shear(ieidx(jc,jb,2),jk-1,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb) +      &
                      shear(ieidx(jc,jb,3),jk-1,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) )
        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    !--------------------------------------------------------------------------
    ! 3) Classical Smagorinsky model with stability correction due to Lilly 1962
    !    at interface cell centers. At this point mech_prod is twice the actual
    !    mechanical production term.
    !--------------------------------------------------------------------------
    ! MP = Mechanical prod term calculated above
    ! visc = mixing_length_sq * SQRT(MP/2) * SQRT(1-Ri/Pr) where
    ! Ri = (g / theta) * d_theta_dz / (MP/2), where Brunt_vaisala_freq (byncy prod term/kh)
    !    = (g / theta) * d_theta_dz.
    ! After simplification: visc = mixing_length_sq/SQRT(2) * SQRT[MP/2 - (Brunt_vaisala_frq/Pr)]
    ! Note that the factor SQRT(2) with mixing_length_sq is considered into the Smag constant
    rl_start   = 3
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2 , nlev
#else
      DO jk = 2 , nlev
        DO jc = i_startidx, i_endidx
#endif
          kh_ic(jc,jk,jb) = rho_ic(jc,jk,jb) * echam_vdf_config(jg)%rturb_prandtl *                     &
                            p_nh_metrics%mixing_length_sq(jc,jk,jb)         *                     &
                            SQRT( MAX( 0._wp, mech_prod(jc,jk,jb) * 0.5_wp -             &
                            echam_vdf_config(jg)%rturb_prandtl * bruvais(jc,jk,jb) ) )
          km_ic(jc,jk,jb) = kh_ic(jc,jk,jb) * echam_vdf_config(jg)%turb_prandtl
        END DO
      END DO
      DO jc = i_startidx, i_endidx
        kh_ic(jc,1,jb)      = kh_ic(jc,2,jb)
        kh_ic(jc,nlevp1,jb) = kh_ic(jc,nlev,jb)
        km_ic(jc,1,jb)      = km_ic(jc,2,jb)
        km_ic(jc,nlevp1,jb) = km_ic(jc,nlev,jb)
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL


    CALL sync_patch_array(SYNC_C, p_patch, kh_ic)
    CALL sync_patch_array(SYNC_C, p_patch, km_ic)

    !--------------------------------------------------------------------------
    !4) Interpolate difusivity (viscosity) to different locations: calculate them for
    !   halos also because they will be used later in diffusion
    !--------------------------------------------------------------------------

    !4a) visc at cell center
    rl_start = grf_bdywidth_c
    rl_end   = min_rlcell_int-1
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1 , nlev
#else
      DO jk = 1 , nlev
        DO jc = i_startidx, i_endidx
#endif
          km_c(jc,jk,jb) = MAX( echam_vdf_config(jg)%km_min, &
                                ( kh_ic(jc,jk,jb) + kh_ic(jc,jk+1,jb) ) * &
                                0.5_wp * echam_vdf_config(jg)%turb_prandtl )
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !4b) visc at vertices
    CALL cells2verts_scalar(kh_ic, p_patch, p_int%cells_aw_verts, km_iv, &
                            opt_rlstart=5, opt_rlend=min_rlvert_int-1)
    km_iv = MAX( echam_vdf_config(jg)%km_min, km_iv * echam_vdf_config(jg)%turb_prandtl )

    !4c) Now calculate visc at half levels at edge
    CALL cells2edges_scalar(kh_ic, p_patch, p_int%c_lin_e, km_ie, &
                            opt_rlstart=grf_bdywidth_e, opt_rlend=min_rledge_int-1)
    km_ie = MAX( echam_vdf_config(jg)%km_min, km_ie * echam_vdf_config(jg)%turb_prandtl )

    !4d)Get visc at the center on interface level
!    prm_diag%tkvm = MAX( echam_vdf_config(jg)%km_min, prm_diag%tkvh * echam_vdf_config(jg)%turb_prandtl )

    rl_start   = 1
    rl_end     = min_rlcell
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev-1
#else
      DO jk = 1, nlev-1
        DO jc = i_startidx, i_endidx
#endif
          pcfm(jc,jk,jb) = km_ic(jc,jk+1,jb)
          pcfh(jc,jk,jb) = kh_ic(jc,jk+1,jb)

        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

  NULLIFY(p_nh_metrics)
  NULLIFY(p_nh_prog)
  NULLIFY(p_int)

  END SUBROUTINE atm_exchange_coeff3d
  !-------------
  !>
  !!

  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_hori_velocity
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for normal velocity component
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - Option to switch on implicit scheme in vertical
  !!
  !! d_vn/d_t =  d_tau_11/d_x1 + d_tau_12/d_x2 + d_tau_13/d_x3
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-05)
  SUBROUTINE diffuse_hori_velocity( nproma,                &
                                  & p_patch,               &
                                  & km_c, km_iv, km_ie,    &
                                  & u_vert, v_vert, div_c, &
                                  & rho, pum1, pvm1, vn,   &
                                  & ddt_u, ddt_v, dt)

    INTEGER,INTENT(in) :: nproma
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch      !< single patch
    TYPE(t_nh_metrics) ,POINTER :: p_nh_metrics
    TYPE(t_int_state)  ,POINTER :: p_int            !< interpolation state
    REAL(wp),           INTENT(in)       :: dt           !< dt turb
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: km_c
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev+1,p_patch%nblks_v) :: km_iv
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev+1,p_patch%nblks_e) :: km_ie
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_v) :: u_vert, v_vert
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: div_c
    REAL(wp), INTENT(INOUT), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: rho, pum1, pvm1

    REAL(wp), INTENT(out) :: ddt_u(nproma,p_patch%nlev,p_patch%nblks_c) !< u tendency
    REAL(wp), INTENT(out) :: ddt_v(nproma,p_patch%nlev,p_patch%nblks_c) !< v tendency

    REAL(wp) :: flux_up_v, flux_dn_v, flux_up_c, flux_dn_c
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4, dvt, inv_dt
    REAL(wp) :: inv_rhoe(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: vn_new(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp), INTENT(IN) :: vn(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: unew(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: vnew(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: tot_tend(nproma,p_patch%nlev,p_patch%nblks_e)

    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, je, jcn, jbn, jvn, jc
    INTEGER :: nlev, jg, nlevp1

    jg = p_patch%id

    p_nh_metrics => p_nh_state(jg)%metrics
    p_int        => p_int_state(jg)

    ! number of vertical levels
    nlev     = p_patch%nlev
    nlevp1   = nlev+1

    inv_dt   = 1._wp / dt

    ividx  => p_patch%edges%vertex_idx
    ivblk  => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    !total tendency
    tot_tend(:,:,:) = 0._wp

    vn_new(:,:,:) = vn(:,:,:)

    CALL sync_patch_array(SYNC_C, p_patch, rho)

    !density at edge
    CALL cells2edges_scalar(rho, p_patch, p_int%c_lin_e, inv_rhoe,                &
                            opt_rlstart=grf_bdywidth_e+1, opt_rlend=min_rledge_int)

    rl_start   = grf_bdywidth_e+1
    rl_end     = min_rledge_int
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
          inv_rhoe(je,jk,jb) = 1._wp / inv_rhoe(je,jk,jb)
        END DO
      END DO
    END DO
!$OMP END DO

    ! 1) First get the horizontal tendencies

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,vn_vert1,vn_vert2,vn_vert3,vn_vert4,&
!$OMP            dvt,jcn,jbn,flux_up_c,flux_dn_c,jvn,flux_up_v,flux_dn_v)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
          vn_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,1)%v1 +   &
                     v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,1)%v2

          vn_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,2)%v1 +   &
                     v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,2)%v2

          vn_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,3)%v1 +   &
                     v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,3)%v2

          vn_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,4)%v1 +   &
                     v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,4)%v2

          dvt      = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,4)%v1   +   &
                     v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,4)%v2   -   &
                     (u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))    *   &
                     p_patch%edges%dual_normal_vert(je,jb,3)%v1   +   &
                     v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,3)%v2)

          ! tendency in normal direction:
          ! flux = visc*(D_11-2/3DIV) = visc*(2*delta_v/(vert_vert_len/2)-2/3*div_of_stress)

          jcn       = iecidx(je,jb,2)
          jbn       = iecblk(je,jb,2)
          flux_up_c = km_c(jcn,jk,jbn) * ( 4._wp * ( vn_vert4 - vn(je,jk,jb) ) *     &
                      p_patch%edges%inv_vert_vert_length(je,jb) - 2._wp * z_1by3 *   &
                      div_c(jcn,jk,jbn) )


          jcn       = iecidx(je,jb,1)
          jbn       = iecblk(je,jb,1)
          flux_dn_c = km_c(jcn,jk,jbn) * ( 4._wp * ( vn(je,jk,jb) - vn_vert3 ) *     &
                      p_patch%edges%inv_vert_vert_length(je,jb) -  2._wp * z_1by3 *  &
                      div_c(jcn,jk,jbn) )

          ! tendency in tangential direction

          ! D_12 between edge center and the vertex: delta_v/(primal_edge_len/2) +
          ! ((vt4+vt2)/2-(vt3+vt2)/2)/(distance_opp_edges)
          ! flux = D_12 * visc

          ! Note that the tangential velocities at vertices are used in D_12 is an
          ! approximation for speed. Better way is to use vt reconstructed from vn at
          ! each edges. Also, visc at somewhere between edge mid point and the vertex
          ! should be used but this is a good approximation

          jvn       = ividx(je,jb,2)
          jbn       = ivblk(je,jb,2)
          flux_up_v = 0.5_wp * ( km_iv(jvn,jk,jbn) + km_iv(jvn,jk+1,jbn) ) *                      &
                      ( p_patch%edges%tangent_orientation(je,jb) *                                &
                        ( vn_vert2 - vn(je,jk,jb) )    *                                          &
                        p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp +                     &
                        dvt * p_patch%edges%inv_vert_vert_length(je,jb) )

          jvn       = ividx(je,jb,1)
          jbn       = ivblk(je,jb,1)
          flux_dn_v = 0.5_wp * ( km_iv(jvn,jk,jbn) + km_iv(jvn,jk+1,jbn) ) *                      &
                      ( p_patch%edges%tangent_orientation(je,jb) *                                &
                        ( vn(je,jk,jb) - vn_vert1 )    *                                          &
                        p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp +                     &
                        dvt * p_patch%edges%inv_vert_vert_length(je,jb) )

          tot_tend(je,jk,jb) = ( ( flux_up_c - flux_dn_c ) *                                      &
                                 p_patch%edges%inv_dual_edge_length(je,jb) +                      &
                                 p_patch%edges%tangent_orientation(je,jb) *                       &
                                 ( flux_up_v - flux_dn_v ) *                                      &
                                 p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp ) *          &
                               inv_rhoe(je,jk,jb)

        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    CALL sync_patch_array(SYNC_E, p_patch, tot_tend)
    CALL rbf_vec_interpol_cell(tot_tend, p_patch, p_int, ddt_u, ddt_v, opt_rlend=min_rlcell_int)

  NULLIFY(p_nh_metrics)
  NULLIFY(p_int)

  END SUBROUTINE diffuse_hori_velocity
  !-------------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_vert_velocity
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for vertical velocity component
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - Option to switch on implicit scheme in vertical
  !! - only solves for jk=2 to nlev. The bottom and top boundaries are left untouched
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-05)
  SUBROUTINE diffuse_vert_velocity( nproma,                &
                                  & p_patch,               &
                                  & rho_ic, w_vert, w_ie,  &
                                  & km_c, km_iv, km_ic,    &
                                  & u_vert, v_vert, div_c, &
                                  & pum1, pvm1, pwp1, vn,  &
                                  & dt)

    INTEGER,INTENT(in) :: nproma
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch      !< single patch
    TYPE(t_nh_metrics) ,POINTER :: p_nh_metrics
    TYPE(t_nh_prog)    ,POINTER :: p_nh_prog     !<the prognostic variables
    TYPE(t_int_state)  ,POINTER :: p_int         !< interpolation state
    REAL(wp),          INTENT(in)        :: km_ic(nproma,p_patch%nlev+1,p_patch%nblks_c)
    REAL(wp),          INTENT(in)        :: dt

    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev+1,p_patch%nblks_c) :: rho_ic
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev+1,p_patch%nblks_v) :: w_vert
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev+1,p_patch%nblks_e) :: w_ie
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: km_c
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev+1,p_patch%nblks_v) :: km_iv
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_v) :: u_vert, v_vert
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: div_c
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: pum1, pvm1
    REAL(wp), INTENT(INOUT), DIMENSION(nproma,p_patch%nlevp1,p_patch%nblks_c) :: pwp1

    REAL(wp) :: flux_up_c, flux_dn_c, dvn1, dvn2, dvt1, dvt2, flux_up_v, flux_dn_v
    REAL(wp) :: vt_e(nproma,p_patch%nlev,p_patch%nblks_e), inv_dt
    REAL(wp), INTENT(IN) :: vn(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp), DIMENSION(nproma,p_patch%nlev) :: a, b, c, rhs
    REAL(wp)                                 :: var_new(p_patch%nlev)

    !interface level variables but only nlev quantities are needed
    REAL(wp) :: hor_tend(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: tot_tend(nproma,p_patch%nlev+1,p_patch%nblks_c)
    REAL(wp) :: inv_rho_ic(nproma,p_patch%nlev,p_patch%nblks_c)!not necessary to allocate for nlev+1

    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk, ieidx, ieblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end, jg
    INTEGER :: jk, jb, je, jc, jcn, jbn, jvn
    INTEGER  :: nlev

    !patch id
    jg = p_patch%id

    p_nh_metrics => p_nh_state(jg)%metrics
    p_nh_prog    => p_nh_state(jg)%prog(nnew(jg))
    p_int        => p_int_state(jg)

    ! number of vertical levels
    nlev     = p_patch%nlev

    inv_dt  = 1._wp / dt

    ividx  => p_patch%edges%vertex_idx
    ivblk  => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    !Some initializations
    a = 0._wp; c = 0._wp
    !total tendency
    tot_tend(:,:,:) = 0._wp

    IF(p_test_run)THEN
      hor_tend(:,:,:) = 0._wp
    END IF

    CALL rbf_vec_interpol_edge( vn, p_patch, p_int, vt_e, opt_rlend=min_rledge_int-1)

    ! Calculate rho at interface for vertical diffusion
    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
#endif
          inv_rho_ic(jc,jk,jb) = 1._wp / rho_ic(jc,jk,jb)
        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL


    ! 1) Get horizontal tendencies at half level edges
    rl_start   = grf_bdywidth_e
    rl_end     = min_rledge_int-1
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,jcn,jbn,dvn1,dvn2,flux_up_c,flux_dn_c,&
!$OMP            jvn,dvt1,dvt2,flux_up_v,flux_dn_v)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO je = i_startidx, i_endidx
#endif

          ! tendency in normal direction
          ! flux = visc_c * D_31_c where D_31(=D_13) is calculated at half level
          ! cell center

          jcn   = iecidx(je,jb,2)
          jbn   = iecblk(je,jb,2)

          dvn2  = pum1(jcn,jk-1,jbn) * p_patch%edges%primal_normal_cell(je,jb,2)%v1 +  &
                  pvm1(jcn,jk-1,jbn) * p_patch%edges%primal_normal_cell(je,jb,2)%v2 -  &
                  pum1(jcn,jk,jbn) * p_patch%edges%primal_normal_cell(je,jb,2)%v1   -  &
                  pvm1(jcn,jk,jbn) * p_patch%edges%primal_normal_cell(je,jb,2)%v2

          flux_up_c = km_ic(jcn,jk,jbn) * (                                            &
                      dvn2 * p_nh_metrics%inv_ddqz_z_half(jcn,jk,jbn) +                &
                      ( w_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) - w_ie(je,jk,jb) ) *  &
                      p_patch%edges%inv_vert_vert_length(je,jb) * 2.0_wp )

          jcn   = iecidx(je,jb,1)
          jbn   = iecblk(je,jb,1)

          dvn1  = pum1(jcn,jk-1,jbn) * p_patch%edges%primal_normal_cell(je,jb,1)%v1 +  &
                  pvm1(jcn,jk-1,jbn) * p_patch%edges%primal_normal_cell(je,jb,1)%v2 -  &
                  pum1(jcn,jk,jbn) * p_patch%edges%primal_normal_cell(je,jb,1)%v1   -  &
                  pvm1(jcn,jk,jbn) * p_patch%edges%primal_normal_cell(je,jb,1)%v2


          flux_dn_c = km_ic(jcn,jk,jbn) * (                                            &
                      dvn1 * p_nh_metrics%inv_ddqz_z_half(jcn,jk,jbn) +                &
                      ( w_ie(je,jk,jb) - w_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) ) *  &
                      p_patch%edges%inv_vert_vert_length(je,jb) * 2.0_wp )


         ! tendency in tangential direction
         ! flux = visc_v * D_32_v where D_32(= D_23) is calculated at half level
         ! between vertex and edge center

          jvn  = ividx(je,jb,2)
          jbn  = ivblk(je,jb,2)

          dvt2 = ( u_vert(jvn,jk-1,jbn) * p_patch%edges%dual_normal_vert(je,jb,2)%v1 +            &
                   v_vert(jvn,jk-1,jbn) * p_patch%edges%dual_normal_vert(je,jb,2)%v2 +            &
                   vt_e(je,jk-1,jb) ) * 0.5_wp  -                                                 &
                 ( u_vert(jvn,jk,jbn) * p_patch%edges%dual_normal_vert(je,jb,2)%v1 +              &
                   v_vert(jvn,jk,jbn) * p_patch%edges%dual_normal_vert(je,jb,2)%v2 +              &
                   vt_e(je,jk,jb) ) * 0.5_wp

          flux_up_v = km_iv(jvn,jk,jbn) * ( dvt2 * p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) +   &
                      p_patch%edges%tangent_orientation(je,jb) * ( w_vert(jvn,jk,jbn) -           &
                      w_ie(je,jk,jb) ) / p_patch%edges%edge_cell_length(je,jb,2) )


          jvn  = ividx(je,jb,1)
          jbn  = ivblk(je,jb,1)

          dvt1 = ( u_vert(jvn,jk-1,jbn) * p_patch%edges%dual_normal_vert(je,jb,1)%v1 +            &
                   v_vert(jvn,jk-1,jbn) * p_patch%edges%dual_normal_vert(je,jb,1)%v2 +            &
                   vt_e(je,jk-1,jb) ) * 0.5_wp           - &
                 ( u_vert(jvn,jk,jbn) * p_patch%edges%dual_normal_vert(je,jb,1)%v1 +              &
                   v_vert(jvn,jk,jbn) * p_patch%edges%dual_normal_vert(je,jb,1)%v2 +              &
                   vt_e(je,jk,jb) ) * 0.5_wp


          flux_dn_v = km_iv(jvn,jk,jbn) * ( dvt1 * p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) +   &
                      p_patch%edges%tangent_orientation(je,jb) * ( w_ie(je,jk,jb) -               &
                      w_vert(jvn,jk,jbn) ) / p_patch%edges%edge_cell_length(je,jb,1) )

          hor_tend(je,jk,jb) = ( flux_up_c - flux_dn_c ) *                                        &
                                 p_patch%edges%inv_dual_edge_length(je,jb) +                      &
                                 p_patch%edges%tangent_orientation(je,jb)  *                      &
                               ( flux_up_v - flux_dn_v ) *                                        &
                                 p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp

        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    ! Interpolate horizontal tendencies to w point: except top and bottom boundaries
    ! w==0 at these boundaries
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
#endif
          tot_tend(jc,jk,jb) = inv_rho_ic(jc,jk,jb)       *                                       &
                               ( hor_tend(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) *                     &
                                 p_int%e_bln_c_s(jc,1,jb) +                                       &
                                 hor_tend(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) *                     &
                                 p_int%e_bln_c_s(jc,2,jb) +                                       &
                                 hor_tend(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) *                     &
                                 p_int%e_bln_c_s(jc,3,jb) )
         END DO
       END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    ! 2) Vertical tendency: evaluated at w point

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx,a,b,c,rhs,var_new)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
     DO jc = i_startidx, i_endidx
       DO jk = 3, nlev-1
#else
     DO jk = 3, nlev-1
       DO jc = i_startidx, i_endidx
#endif

           a(jc,jk)   = - 2._wp * km_c(jc,jk-1,jb) * p_nh_metrics%inv_ddqz_z_full(jc,jk-1,jb) * &
                          p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb)

           c(jc,jk)   = - 2._wp * km_c(jc,jk,jb) * p_nh_metrics%inv_ddqz_z_full(jc,jk,jb)     * &
                          p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb)

           b(jc,jk)   =  inv_dt - a(jc,jk) - c(jc,jk)

           rhs(jc,jk) =  pwp1(jc,jk,jb) * inv_dt +                                       &
                         2._wp * ( km_c(jc,jk,jb)   * z_1by3 * div_c(jc,jk,jb)     -            &
                                   km_c(jc,jk-1,jb) * z_1by3 * div_c(jc,jk-1,jb) ) *            &
                         p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb)
       END DO
     END DO

        ! Boundary treatment
        !--------------------------------------------------------
        ! jk = 2 (w == 0)
        !--------------------------------------------------------
        DO jc = i_startidx, i_endidx
          c(jc,2)   = - 2._wp * km_c(jc,2,jb) * p_nh_metrics%inv_ddqz_z_full(jc,2,jb) *         &
                        p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * inv_rho_ic(jc,2,jb)

          b(jc,2)   = inv_dt - c(jc,2) + 2._wp * km_c(jc,1,jb) *                                &
                      p_nh_metrics%inv_ddqz_z_full(jc,1,jb) *                                   &
                      p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * inv_rho_ic(jc,2,jb)

          rhs(jc,2) = pwp1(jc,2,jb) * inv_dt +                                           &
                      2._wp * ( km_c(jc,2,jb) * z_1by3 * div_c(jc,2,jb) -                       &
                                km_c(jc,1,jb) * z_1by3 * div_c(jc,1,jb) ) *                     &
                      p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * inv_rho_ic(jc,2,jb)
        END DO
        !--------------------------------------------------------
        ! jk = nlev (w == 0)
        !--------------------------------------------------------
        DO jc = i_startidx, i_endidx
          a(jc,nlev)   = - km_c(jc,nlev-1,jb) * p_nh_metrics%inv_ddqz_z_full(jc,nlev-1,jb) *    &
                           p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * 2._wp *                   &
                           inv_rho_ic(jc,nlev,jb)

          b(jc,nlev)   =   inv_dt - a(jc,nlev) + 2._wp * km_c(jc,nlev,jb) *                     &
                           p_nh_metrics%inv_ddqz_z_full(jc,nlev,jb) *                           &
                           p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * inv_rho_ic(jc,nlev,jb)

          rhs(jc,nlev) =   pwp1(jc,nlev,jb) * inv_dt +                                   &
                           2._wp * ( km_c(jc,nlev,jb) * z_1by3 * div_c(jc,nlev,jb) -            &
                                     km_c(jc,nlev-1,jb) * z_1by3 * div_c(jc,nlev-1,jb) ) *      &
                           p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * inv_rho_ic(jc,nlev,jb)
        END DO

        ! CALL TDMA
        DO jc = i_startidx, i_endidx
          CALL tdma_solver( a(jc,2:nlev), b(jc,2:nlev), c(jc,2:nlev), rhs(jc,2:nlev),           &
                            nlev-1,var_new(2:nlev) )

          tot_tend(jc,2:nlev,jb) = tot_tend(jc,2:nlev,jb) +                                     &
                                   ( var_new(2:nlev) - pwp1(jc,2:nlev,jb) ) * inv_dt
       END DO

    END DO !jb
!$OMP END DO
!$OMP END PARALLEL

    ! update w
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
#endif
          pwp1(jc,jk,jb) = pwp1(jc,jk,jb) + dt * tot_tend(jc,jk,jb)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


   CALL sync_patch_array(SYNC_C, p_patch, pwp1)

  NULLIFY(p_nh_metrics)
  NULLIFY(p_nh_prog)
  NULLIFY(p_int)

  END SUBROUTINE diffuse_vert_velocity
  !-------------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_scalar
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for cell based scalars
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - Option to switch on implicit scheme in vertical
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-05)
  SUBROUTINE diffuse_scalar( nproma, var_temp, &
                           & p_patch,          &
                           & kh_ic, km_ie,     &
                           & hori_tend,        &
                           & rho,              &
                           & scalar_name)

    INTEGER,INTENT(in) :: nproma
    REAL(wp),          INTENT(in)           :: var_temp(:,:,:)      ! input scalar
    TYPE(t_patch),     INTENT(inout),TARGET :: p_patch         !< single patch
    TYPE(t_nh_metrics) ,POINTER :: p_nh_metrics
    TYPE(t_int_state)  ,POINTER :: p_int         !< interpolation state
    REAL(wp),          INTENT(in)           :: rho(:,:,:)      ! density at cell center
    INTEGER,           INTENT(in)           :: scalar_name
    REAL(wp), INTENT(IN)                    :: kh_ic(nproma,p_patch%nlev+1,p_patch%nblks_c)
    REAL(wp), INTENT(IN)                    :: km_ie(nproma,p_patch%nlev+1,p_patch%nblks_e)
    REAL(wp), INTENT(OUT)                   :: hori_tend(nproma,p_patch%nlev,p_patch%nblks_c) !< total tendency
    REAL(wp)                                :: var(nproma,p_patch%nlev,p_patch%nblks_c)      ! input scalar

    !Local variables
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, je, jc, jg
    INTEGER :: nlev, nlevp1

    INTEGER,  DIMENSION(:,:,:), POINTER :: iecidx, iecblk, ieidx, ieblk

    REAL(wp) :: nabla2_e(nproma,p_patch%nlev,p_patch%nblks_e)

    INTEGER, PARAMETER :: tracer_dry_static = 1
    INTEGER, PARAMETER :: tracer_water = 2


    !patch id
    jg = p_patch%id

    p_nh_metrics => p_nh_state(jg)%metrics
    p_int        => p_int_state(jg)

    ! number of vertical levels
    nlev = p_patch%nlev
    nlevp1 = nlev+1

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    hori_tend = 0._wp

    var = var_temp
    CALL sync_patch_array(SYNC_C, p_patch, var)

    !1) First set local vars to 1 for other scalars
    !   Soon get different routines for different scalars

!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)

    !---------------------------------------------------------------
    ! Horizontal diffusion (conservative; following mo_nh_diffusion)
    !---------------------------------------------------------------

    !include halo points and boundary points because these values will be
    !used in next loop
    rl_start   = grf_bdywidth_e
    rl_end     = min_rledge_int-1
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jk,je,jb,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      ! compute kh_ie * grad_horiz(var)
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
          nabla2_e(je,jk,jb) = 0.5_wp * ( km_ie(je,jk,jb) + km_ie(je,jk+1,jb) ) *               &
                               echam_vdf_config(jg)%rturb_prandtl *                                   &
                               p_patch%edges%inv_dual_edge_length(je,jb) *                      &
                               ( var(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -                      &
                                 var(iecidx(je,jb,1),jk,iecblk(je,jb,1)) )
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO

    ! now compute the divergence of the quantity above
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
#endif
          ! horizontal tendency
          hori_tend(jc,jk,jb) = (                                                  &
                        nabla2_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%geofac_div(jc,1,jb) +  &
                        nabla2_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%geofac_div(jc,2,jb) +  &
                        nabla2_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%geofac_div(jc,3,jb) ) / rho(jc,jk,jb)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    IF (is_dry_cbl .AND. scalar_name==tracer_water) THEN
!$OMP PARALLEL
      CALL init(hori_tend(:,:,:))
!$OMP END PARALLEL
    END IF

  NULLIFY(p_nh_metrics)
  NULLIFY(p_int)

  END SUBROUTINE diffuse_scalar
  !-------------------------------------------------------------------------------------


  !
  !! stability_function_mom
  !! Taken from COSMO docs and Holstag & Boville 1992
  !!------------------------------------------------------------------------
  FUNCTION stability_function_mom(RIB, hz0, tc) RESULT(stab_fun)
     REAL(wp), INTENT(IN) :: RIB, hz0, tc

     REAL(wp) :: stab_fun, hz0_fac

     IF(RIB.GE.0._wp)THEN
       !Cosmo
       !stab_fun = 1._wp / ( 1._wp + 10._wp*RIB/SQRT(1._wp+5*RIB) )

       !H&B
       stab_fun = 1._wp / ( 1._wp + 10._wp*RIB*(1._wp+8._wp*RIB) )
    ELSE
       hz0_fac = ( max(hz0, 1._wp)**(1._wp/3._wp) - 1._wp )**1.5_wp ! FLO - hz0 can be < 1 then, the **1.5 is invalid.
       !for water surface (z0/h)**(1/3)<<1 giving hz0_fac=SQRT(h/z0)
       !Generally it is explicitly written for water surface but i don't
       !see any reason to do that.
       stab_fun = 1._wp + 10._wp*ABS(RIB)/(1._wp + 75._wp*tc*hz0_fac*SQRT(ABS(RIB)))
     END IF

  END FUNCTION stability_function_mom
  !>
  !! stability_function_heat
  !!------------------------------------------------------------------------
  FUNCTION stability_function_heat(RIB, hzh, tc) RESULT(stab_fun)
     REAL(wp), INTENT(IN) :: RIB, hzh, tc

     REAL(wp) :: stab_fun, hzh_fac

     IF(RIB.GE.0._wp)THEN
       !Cosmo
       !stab_fun = 1._wp / ( 1._wp + 15._wp*RIB*SQRT(1._wp+5*RIB) )

       !H&B
       stab_fun = 1._wp / ( 1._wp + 10._wp*RIB*(1._wp+8._wp*RIB) )
     ELSE
       hzh_fac = ( max(hzh, 1._wp)**(1._wp/3._wp) - 1._wp )**1.5_wp
       !for water surface (zh/h)**(1/3)<<1 giving hzh_fac=SQRT(h/zh)
       !Generally it is explicitly written for water surface but i don't
       !see any reason to do that.
       stab_fun = 1._wp + 15._wp*ABS(RIB)/(1._wp + 75._wp*tc*hzh_fac*SQRT(ABS(RIB)))
     END IF
  END FUNCTION stability_function_heat

  !>
  !! factor_mom
  !!------------------------------------------------------------------------
  !! Businger Dyer similarity profile:
  !! Louis (1979) A Parametirc model of vertical eddy fluxes in the atmosphere
  !! and R. B. Stull's book
  !!------------------------------------------------------------------------
  FUNCTION businger_mom(z0, z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z0, z1, L
     REAL(wp) :: factor, zeta, psi, lamda
     REAL(wp) :: zeta0, psi0, lamda0

     IF(L > 0._wp)THEN !Stable
       zeta  = z1/L
       zeta0 = z0/L
       IF(zeta > 1._wp)THEN !Zeng etal 1997 J. Clim
         psi    = -bsm*LOG(zeta) - zeta + 1
         psi0   = -bsm*LOG(zeta0) - zeta0 + 1
         factor = (LOG(L/z0) + bsh - psi + psi0  ) / ckap
       ELSE
         psi  = -bsm*zeta
         psi0 = -bsm*zeta0
         factor = ( LOG(z1/z0) - psi + psi0 ) / ckap
       END IF
     ELSEIF(L < 0._wp)THEN !unstable
       zeta   = z1/L
       zeta0  = z0/L
       lamda  = SQRT(SQRT(1._wp - bum*zeta))
       lamda0 = SQRT(SQRT(1._wp - bum*zeta0))

       psi    = 2._wp * LOG(1._wp+lamda) + LOG(1._wp+lamda*lamda) - &
                2._wp * ATAN(lamda) + pi_2 - 3._wp*ln2

       psi0   = 2._wp * LOG(1._wp+lamda0) + LOG(1._wp+lamda0*lamda0) - &
                2._wp * ATAN(lamda0) + pi_2 - 3._wp*ln2

       factor = ( LOG(z1/z0) - psi + psi0 ) / ckap
     ELSE !neutral
       factor = LOG(z1/z0) / ckap
     END IF

  END FUNCTION businger_mom
  !>
  !! factor_heat
  !!------------------------------------------------------------------------
  !! Businger Dyer similarity profile:
  !! Louis (1979) A Parametirc model of vertical eddy fluxes in the atmosphere
  !! and R. B. Stull's book
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-06)
  FUNCTION businger_heat(z0, z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z0, z1, L
     REAL(wp) :: factor, zeta, lamda, psi
     REAL(wp) :: zeta0, lamda0, psi0

     IF(L > 0._wp)THEN !Stable
       zeta   = z1/L
       zeta0  = z0/L
       IF(zeta > 1._wp)THEN !Zeng etal 1997 J. Clim
         psi    = -bsh*LOG(zeta) - zeta + 1
         psi0   = -bsh*LOG(zeta0) - zeta0 + 1
         factor = (LOG(L/z0) + bsh - psi + psi0  ) / ckap
       ELSE
         psi    = -bsh*zeta
         psi0   = -bsh*zeta0
         factor = (LOG(z1/z0) - psi + psi0) / ckap
       END IF
     ELSEIF(L < 0._wp)THEN !unstable
       zeta   = z1/L
       zeta0  = z0/L
       lamda  = SQRT(1._wp - buh*zeta)
       lamda0 = SQRT(1._wp - buh*zeta0)
       psi    = 2._wp * ( LOG(1._wp+lamda) - ln2 )
       psi0   = 2._wp * ( LOG(1._wp+lamda0) - ln2 )
       factor = (LOG(z1/z0) - psi + psi0) / ckap
     ELSE !Neutral
       factor = LOG(z1/z0) / ckap
     END IF

  END FUNCTION businger_heat


  !-------------
END MODULE mo_sma_turbulence_diag
