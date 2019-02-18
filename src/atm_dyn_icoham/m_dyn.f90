!>
!!  contains the vertical discretization.
!!
!!
!! @par Revision History
!!  Original version: SUBROUTIN dyn in ECHAM5
!!  Modified for the ICON hydrostatic dynamical core by Hui Wan, MPI-M (2006-02)
!!  Modified for the ICON hydrostatic dynamical core by Hui Wan, MPI-M (2006-08)
!!  Alternative formulations for the pressure gradient force and
!!  vorticity flux term added by Hui Wan, MPI-M (2006-08-18).
!!  Vertical advection of normal wind corrected by Hui Wan, MPI-M (2006-08-25)
!!  <i>SUBROUTINE vortflux_ecnsv</i> added by Hui Wan, MPI-M (2006-09-15)
!! @par
!!  Modifications by Hui Wan, MPI-M (2008-04-05)
!!  - ldivavg, divavg_c0 and divavg_cj renamed ldiv_avg, div_avg_c0
!!    and div_avg_c2, respectively.
!!  - ldivdamp and kdivdamp renamed ldiv_damp and kdiv_damp, respectively.
!!  - lhnormavg and lknormavg renamed lgrad_avg_h and lgrad_avg_kin,
!!    respctively.
!! @par
!!  Modifications by Hui Wan, MPI-M (2008-04-24)
!!  In order to avoid if- and case-blocks in the do-loops:
!!   * The additional diagnositcs for the tendency terms were temporarily
!!     removed.
!!   * Some parts of the subroutine dyn in which smoothing of the divergence
!!     operator was applied were restructured.
!! @par
!!  Modifications by Almut Gassmann, MPI-M, (2008-09-11)
!!  - Cleaning and optimizing
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
MODULE m_dyn
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
  USE mo_physical_constants, ONLY: cpd, rd, grav, vtmpc1, p0ref
  USE mo_model_domain,       ONLY: t_patch
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_math_gradients,     ONLY: grad_fd_norm
  USE mo_math_divrot,        ONLY: div, div_avg
  USE mo_dynamics_config,    ONLY: idiv_method
  USE mo_ha_dyn_config,      ONLY: ha_dyn_config
  USE mo_io_config,          ONLY: l_outputtime
  USE mo_parallel_config,    ONLY: nproma
  USE mo_run_config,         ONLY: nlev, nlevm1, nlevp1,iqv, iforcing, &
                                   iqm_max, output_mode
  USE mo_icoham_dyn_types,   ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_intp,               ONLY: cell_avg,                    &
                                 & cells2edges_scalar, edges2cells_scalar, &
                                 & verts2edges_scalar, cells2verts_scalar
  USE mo_interpol_config,    ONLY: i_cori_method, sick_a, sick_o
  USE mo_nonlinear_adv,      ONLY: kin_vel_rot, lamb_rot
  USE mo_eta_coord_diag,     ONLY: half_level_pressure, full_level_pressure, &
                                   auxhyb, geopot
  USE mo_impl_constants,     ONLY: iecham,ildf_echam, inwp
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_vertical_coord_table, ONLY: delpr, nplev, nplvp1,             &
    &                                vct_b, alrrdic, rdlnp0i, rdt0ral, t0icao

  USE mo_sync,               ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
  
  USE mo_timer,              ONLY: ltimer, timer_start, timer_stop,&
    & timer_dyn_theta

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: dyn_theta

  CONTAINS

!-------------------------------------------------------------------------
  !>
  !!               Computes adiabatic tendencies.
  !!
  !!               Version using potential temperature as prognostic variable
  !!
  !! Method:
  !!
  !! The primary purpose is to compute adiabatic tendencies.
  !! Auxilliary variables connected with the vertical difference
  !! scheme are saved as they are subsequently required to compute
  !! input to physical parameterizations.
  !!
  !! Externals:
  !! *pres* and *auxhyb* are called to calculate auxilliary variables.
  !! *geopot* is called to calculate geopotential deviations.
  !!
  !! @par Revision History
  !!   This subroutine originated from <i>SUBROUTINE dyn</i>
  !!   Initial version by Guenther Zaengl, DWD, 2009-01-13
  !!
  SUBROUTINE dyn_theta( pt_patch, pt_int_state, pt_ext_data, & ! input
    &                   pt_prog, pt_diag , pt_tend_dyn       ) ! output
!

  IMPLICIT NONE


  TYPE(t_patch), TARGET, INTENT(INOUT) :: pt_patch     !< grid/patch info.
  TYPE(t_int_state),TARGET,INTENT(IN):: pt_int_state !< horizontal interpolation coeff.
  TYPE(t_external_data), INTENT(IN) :: pt_ext_data  !< external data

  TYPE(t_hydro_atm_prog),INTENT(INOUT) :: pt_prog     !< prognostic variables
  TYPE(t_hydro_atm_diag),INTENT(INOUT) :: pt_diag     !< diagnostic variables
  TYPE(t_hydro_atm_prog),INTENT(INOUT) :: pt_tend_dyn !< tendency of the prognostic variables

! Local array bounds:

  INTEGER :: nlen                    ! arb. field length
  INTEGER :: nblks_c, nblks_e        ! number of blocks for cells / edges
  INTEGER :: npromz_e, npromz_c      ! length of last block line

! index ranges needed for grid refinement
  INTEGER :: i_startblk, i_startidx, i_endidx

! Local scalars:

  INTEGER  :: jb, jk, je, jc, jt ! loop indices
  INTEGER  :: jkp                ! vertical level indices
  REAL(wp) :: rovcp              ! R/cp

! Local arrays:
! ( *_c* and *_e* denote cells and edges, respectively.
!   *_i* and *_m* denote half and full vertical levels,
!                   i.e. layer interfaces and layers.)

! arrays associated to cells

  REAL(wp), DIMENSION ( nproma, nlev, pt_patch%nblks_c ) :: &
            & z_tv_c,      & ! virtual temperature
            & z_tvp_c,     & ! virtual temperature perturbation
            & z_tbar_mc,   & ! reference temperature
            & z_zbar_mc,   & ! reference geopotential (full level)
            & z_geop_mc,   & ! geopot. perturb.(full levels)
            & z_mdiv,      & ! mass divergence
            & z_aux,       & ! auxiliary field for smoothing operations
            & z_tmp_mc,    & ! help array on cells
            & z_tmdiv,     & ! divergence ( V.T)
            & z_rdlnp_mc,  & ! Rd.log pressure on cells
            & z_dpdt,      & ! adiaba conversion term
            & z_e_mech_c     ! mechanical energy on cells

  REAL(wp), DIMENSION ( nproma, nlevp1, pt_patch%nblks_c ):: &
            & z_lnp_ic,    & ! log(pressure) on half levels
            & z_zbar_ic,   & ! reference geopotential (half level)
            & z_theta_ic,  & ! potential temperature at interface levels
            & z_geop_ic,   & ! geopot. perturb.(half levels)
            & z_mdiv_int     ! accummulated mass divergence

  REAL(wp), DIMENSION ( nproma, pt_patch%nblks_c ) :: &
            & z_geos_c       ! modified surface geopotential

  REAL(wp), DIMENSION ( nproma ) :: z_help

! arrays associated to edges
  REAL(wp), DIMENSION ( nproma, nlev, pt_patch%nblks_e ) :: &
            & z_tvp_e,       & ! virtual temperature perturb.
            & z_tv_e,        & ! virtual temperature
            & z_tmp_me,      & ! temporary array
            & z_temp_e,      & ! temperature at edges
            & z_temp_flux_e, & ! V.T at edges
            & z_pres_grad1,  & ! part of horizontal pressure gradient
            & z_pres_grad2,  & ! part of horizontal pressure gradient
            & z_pres_grad4,  & ! part of horizontal pressure gradient
            & z_rdlnp_me       ! Rd.log pressure on cells

  REAL(wp), DIMENSION ( nproma, nlevp1, pt_patch%nblks_e ) :: &
            & z_weta_e,      & ! weta averaged to edges
            & z_tmp_e          ! edge help value

  REAL(wp), DIMENSION ( nproma, nlev, pt_patch%nblks_v ) :: &
            & z_tmp_v          ! vertex help value

  REAL(wp):: z_aux_tracer(nproma,nlev,pt_patch%nblks_c)! needed for virt. increment

!  Executable statements
  IF (ltimer) CALL timer_start(timer_dyn_theta)

!-------------------------------------------------------------------------------
! 0. Preliminary calculations
!-------------------------------------------------------------------------------

   nblks_c   = pt_patch%nblks_c
   npromz_c  = pt_patch%npromz_c
   nblks_e   = pt_patch%nblks_e
   npromz_e  = pt_patch%npromz_e

   SELECT CASE (iforcing)
   CASE (iecham,ildf_echam, inwp)  ! real physics
     z_aux_tracer(:,:,:) = 0.0_wp
   CASE DEFAULT
     ! nothing to be done
   END SELECT


  ! 0.1 Compute half-level/full-level pressure and auxiliary variables.

   rovcp = rd/cpd

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = 1,nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      CALL half_level_pressure( pt_prog%pres_sfc(:,jb), nproma, nlen, &! in
                                pt_diag%pres_ic(:,:,jb)              ) ! out

      CALL full_level_pressure( pt_diag%pres_ic(:,:,jb),nproma, nlen, &! in
                                pt_diag%pres_mc(:,:,jb)              ) ! out

      CALL auxhyb( pt_diag%pres_ic(:,:,jb), nproma, nlen,           &! in
                   pt_diag%delp_c(:,:,jb), pt_diag%rdelp_c(:,:,jb), &! out
                   z_lnp_ic(:,:,jb), pt_diag%rdlnpr_c(:,:,jb),      &! out
                   pt_diag%rdalpha_c(:,:,jb)                       ) ! out

   ! Compute exner function
      pt_diag%exner(1:nlen,1:nlev,jb) = &
        EXP(rovcp*LOG(pt_diag%pres_mc(1:nlen,1:nlev,jb)/p0ref))
   ! Diagnose physical temperature
      pt_prog%temp(1:nlen,1:nlev,jb) = pt_prog%theta(1:nlen,1:nlev,jb)* &
        pt_diag%rdelp_c(1:nlen,1:nlev,jb)* pt_diag%exner(1:nlen,1:nlev,jb)

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  ! 0.2 layer thickness at edge centres
  !     (All the loops are locate inside the "cells2edges_scalar" subroutine.)

   CALL cells2edges_scalar (pt_diag%delp_c, pt_patch, & ! input
                            pt_int_state%c_lin_e, &
                            pt_diag%delp_e,     & ! output
                            nplvp1, nlev )        ! optional input

   IF(i_cori_method >= 2) THEN

     z_tmp_v=0.0_wp
     CALL cells2verts_scalar (pt_diag%delp_c, pt_patch,    & ! input
                              pt_int_state%cells_aw_verts, & ! area weight
                              z_tmp_v,                     & ! output
                              nplvp1, nlev )                 ! optional input
     CALL sync_patch_array(SYNC_V,pt_patch,z_tmp_v)
     CALL verts2edges_scalar (z_tmp_v, pt_patch,    & ! input
                              pt_int_state%v_1o2_e, & ! direct distance ave.
                              z_tmp_e,       & ! output
                              nplvp1, nlev )          ! optional input

!$OMP PARALLEL
     i_startblk = pt_patch%edges%start_blk(2,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk) ICON_OMP_DEFAULT_SCHEDULE
     DO jb = i_startblk,nblks_e
       CALL get_indices_e(pt_patch, jb, i_startblk, nblks_e, &
                          i_startidx, i_endidx, 2)
       DO jk = nplvp1,nlev
         pt_diag%delp_e(i_startidx:i_endidx,jk,jb)  =  &
         &   sick_a*z_tmp_e(i_startidx:i_endidx,jk,jb) &
         &  +sick_o*pt_diag%delp_e(i_startidx:i_endidx,jk,jb)
       ENDDO

     ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   ENDIF

  ! 0.3 Mass flux at edges

!$OMP PARALLEL
   i_startblk = pt_patch%edges%start_blk(2,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk,nblks_e

     CALL get_indices_e(pt_patch, jb, i_startblk, nblks_e, &
                        i_startidx, i_endidx, 2)

      DO jk = 1,nplev   ! only constant layer thickness
        pt_diag%delp_e(i_startidx:i_endidx,jk,jb)  = delpr(jk)
      ENDDO

      pt_diag%mass_flux_e(i_startidx:i_endidx,:,jb) = &
              pt_diag%delp_e(i_startidx:i_endidx,:,jb) *    &
              pt_prog%vn(i_startidx:i_endidx,:,jb)

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   CALL sync_patch_array(SYNC_E,pt_patch, pt_diag%mass_flux_e)

!-------------------------------------------------------------------------------
! 1. The continuity equation
!-------------------------------------------------------------------------------
!(Sum mass divergence to compute surface pressure tendency;
! Compute the vertical velocity eta-dot at layer interfaces.)
!-------------------------------------------------------------------------------

  !-1.1 Compute mass divergence
  SELECT CASE(idiv_method)

  CASE(1)

    CALL div(pt_diag%mass_flux_e, pt_patch, pt_int_state, &  ! input
    &        z_mdiv, opt_rlstart=2 )           ! output

  CASE(2)

    CALL div_avg(pt_diag%mass_flux_e, pt_patch, pt_int_state, &
    &            pt_int_state%c_bln_avg, z_mdiv, opt_rlstart=2 )

  END SELECT

!$OMP PARALLEL PRIVATE(i_startblk)
   i_startblk = pt_patch%cells%start_blk(2,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jkp) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk,nblks_c

     CALL get_indices_c(pt_patch, jb, i_startblk, nblks_c, &
                        i_startidx, i_endidx, 2)

     !-1.3 vertically integrated mass divergence

      z_mdiv_int(i_startidx:i_endidx,1,jb) = 0.0_wp
      DO jk = 1,nlev
         jkp = jk + 1
         z_mdiv_int(i_startidx:i_endidx,jkp,jb) =       &
                z_mdiv_int(i_startidx:i_endidx,jk,jb) + &
                z_mdiv    (i_startidx:i_endidx,jk,jb)
      ENDDO

     !-1.4 Vertical velocity at the interfaces (half-levels)

      pt_diag%weta(i_startidx:i_endidx,1,jb)      = 0._wp   ! upper boundary
      pt_diag%weta(i_startidx:i_endidx,nlevp1,jb) = 0._wp   ! lower boundary

      DO jk = 2,nlev
        pt_diag%weta(i_startidx:i_endidx,jk,jb)  = &
  &     -z_mdiv_int(i_startidx:i_endidx,jk,jb) + &
  &     z_mdiv_int(i_startidx:i_endidx,nlevp1,jb)*vct_b(jk)
      ENDDO

   ENDDO
!$OMP END DO

   ! For nested domains, tendencies are interpolated from the parent domain
   ! on a boundary zone with a width of grf_bdywidth_c for cells and
   ! grf_bdywidth_e for edges, respectively. These tendencies must not be
   ! overwritten.
   i_startblk = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk,nblks_c

     CALL get_indices_c(pt_patch, jb, i_startblk, nblks_c, &
                        i_startidx, i_endidx, grf_bdywidth_c+1)

     !-1.5 Tendency of surface pressure

      pt_tend_dyn%pres_sfc(i_startidx:i_endidx,jb) = &
 &             -z_mdiv_int(i_startidx:i_endidx,nlevp1,jb)
   ENDDO
!$OMP END DO


!-------------------------------------------------------------------------------
! 2. Compute reference temperature/geopotential and the deviation of
!    virtual temprature
!-------------------------------------------------------------------------------
  ! 2.1 the reference state for temperature and geopotential
  !     ( ICAO(1964) standard atmosphere.
  !       Ref. MPI-M Report 349, p20 )

   i_startblk = pt_patch%cells%start_blk(2,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,z_help,jkp) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk,nblks_c

     CALL get_indices_c(pt_patch, jb, i_startblk, nblks_c, &
                        i_startidx, i_endidx, 2)

      IF (ha_dyn_config%lref_temp ) THEN

         DO jk = 1,nlev
            z_help(i_startidx:i_endidx)           = EXP(alrrdic* &
 &            (rd*LOG(pt_diag%pres_mc(i_startidx:i_endidx,jk,jb))-rdlnp0i))
            z_zbar_mc(i_startidx:i_endidx,jk,jb)  = &
               -rdt0ral*z_help(i_startidx:i_endidx)
            z_tbar_mc(i_startidx:i_endidx,jk,jb)  = &
               t0icao*z_help(i_startidx:i_endidx)

            jkp = jk + 1
            z_zbar_ic(i_startidx:i_endidx,jkp,jb) = -rdt0ral*EXP( alrrdic* &
 &            (rd*LOG(pt_diag%pres_ic(i_startidx:i_endidx,jkp,jb))-rdlnp0i))
         ENDDO

      ELSE

         z_tbar_mc(:,:,jb) = 0.0_wp
         z_zbar_mc(:,:,jb) = 0.0_wp
         z_zbar_ic(:,:,jb) = 0.0_wp

      ENDIF

  ! 2.2 (The deviation of) virtual temperature (from the reference state).

   SELECT CASE (iforcing)
   CASE (iecham, ildf_echam, inwp)  ! real physics with moist atmosphere

     z_aux_tracer(:,:,jb) = 0._wp

     DO jt=2,iqm_max
       DO jk = 1, nlev
         DO jc = i_startidx,i_endidx
           z_aux_tracer(jc,jk,jb) = z_aux_tracer(jc,jk,jb) + pt_prog%tracer(jc,jk,jb,jt)
         END DO
       END DO
     END DO

     pt_diag%virt_incr(i_startidx:i_endidx,:,jb) =                     &
          &  vtmpc1*pt_prog%tracer(i_startidx:i_endidx,:,jb,iqv)       &
          &  - z_aux_tracer(i_startidx:i_endidx,:,jb)
     z_tv_c(i_startidx:i_endidx,:,jb) =                                &
          &  pt_prog%temp(i_startidx:i_endidx,:,jb)                    &
          &  * ( 1.0_wp + pt_diag%virt_incr(i_startidx:i_endidx,:,jb))
     z_tvp_c(i_startidx:i_endidx,:,jb) =                               &
          &  z_tv_c(i_startidx:i_endidx,:,jb)                          &
          &  - z_tbar_mc(i_startidx:i_endidx,:,jb)

   CASE DEFAULT ! (Held Suarez or no forcing)

     z_tv_c(i_startidx:i_endidx,:,jb) =              &
          &  pt_prog%temp(i_startidx:i_endidx,:,jb)
     z_tvp_c(i_startidx:i_endidx,:,jb) =             &
          &  pt_prog%temp(i_startidx:i_endidx,:,jb)  &
          &  - z_tbar_mc(i_startidx:i_endidx,:,jb)

   END SELECT

   ENDDO
!$OMP END DO

!-------------------------------------------------------------------------------
! 3. Compute vertical advection
!-------------------------------------------------------------------------------
  ! 3.1 vertical advection of temperature


   ! For nested domains, tendencies are interpolated from the parent domain
   ! on a boundary zone with a width of grf_bdywidth_c for cells and
   ! grf_bdywidth_e for edges, respectively. These tendencies must not be
   ! overwritten.
   i_startblk = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jkp) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk,nblks_c

     CALL get_indices_c(pt_patch, jb, i_startblk, nblks_c, &
                        i_startidx, i_endidx, grf_bdywidth_c+1)

     ! Uppermost interface contribution vanishes
      pt_tend_dyn%temp(i_startidx:i_endidx,1,jb) = 0._wp
      z_theta_ic(i_startidx:i_endidx,1,jb) = 0._wp
      z_theta_ic(i_startidx:i_endidx,nlevp1,jb) = 0._wp

      DO jk = 1, nlevm1

        jkp = jk + 1

        z_theta_ic(i_startidx:i_endidx,jkp,jb) = 0.5_wp*( &
          pt_prog%theta(i_startidx:i_endidx,jk,jb)*       &
          pt_diag%rdelp_c(i_startidx:i_endidx,jk,jb) +    &
          pt_prog%theta(i_startidx:i_endidx,jkp,jb)*       &
          pt_diag%rdelp_c(i_startidx:i_endidx,jkp,jb) )

      ENDDO

      pt_tend_dyn%temp(i_startidx:i_endidx,1:nlev,jb) = &
        pt_diag%weta(i_startidx:i_endidx,1:nlev,jb)*            &
        z_theta_ic(i_startidx:i_endidx,1:nlev,jb) - &
        pt_diag%weta(i_startidx:i_endidx,2:nlev+1,jb)*            &
        z_theta_ic(i_startidx:i_endidx,2:nlev+1,jb)
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  ! 3.2 vertical advection of normal wind

  ! First average the vertical velocity (m*eta-dot) from cell centers to edges
  ! (The block loop locates inside the "cells2edges_scalar" subroutine.)

   CALL cells2edges_scalar (pt_diag%weta, pt_patch, & ! input
                            pt_int_state%c_lin_e, &
                            z_weta_e, opt_rlstart=4)  ! output

  ! Then use the same formula as for temperature, but at edge centers.

   ! For nested domains, tendencies are interpolated from the parent domain
   ! on a boundary zone with a width of grf_bdywidth_c for cells and
   ! grf_bdywidth_e for edges, respectively. These tendencies must not be
   ! overwritten.
!$OMP PARALLEL PRIVATE(i_startblk)
   i_startblk = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jkp) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk,nblks_e

     CALL get_indices_e(pt_patch, jb, i_startblk, nblks_e, &
                        i_startidx, i_endidx, grf_bdywidth_e+1)

      pt_tend_dyn%vn(i_startidx:i_endidx,1,jb) = 0._wp

      DO jk = 1, nlevm1

         jkp = jk + 1

         pt_tend_dyn%vn(i_startidx:i_endidx,jkp,jb) = &
 &             z_weta_e(i_startidx:i_endidx,jkp,jb)   &
 &        *( pt_prog%vn(i_startidx:i_endidx,jk ,jb)   &
 &          -pt_prog%vn(i_startidx:i_endidx,jkp,jb) )

         pt_tend_dyn%vn(i_startidx:i_endidx,jk ,jb) = &
 &       pt_tend_dyn%vn(i_startidx:i_endidx,jkp,jb)   &
 &      +pt_tend_dyn%vn(i_startidx:i_endidx,jk ,jb)

      ENDDO

      pt_tend_dyn%vn(i_startidx:i_endidx,1:nlev,jb) =     &
 &    pt_tend_dyn%vn(i_startidx:i_endidx,1:nlev,jb)       &
 &   /pt_diag%delp_e(i_startidx:i_endidx,1:nlev,jb)*0.5_wp

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!-------------------------------------------------------------------------------
! 4. Compute horizontal advection of potential temperature
!-------------------------------------------------------------------------------

  ! The only term is here div(theta*delp*V)

  !interpolate temperature from cells to edges
   CALL cells2edges_scalar( pt_prog%theta, pt_patch, & !input
                            pt_int_state%c_lin_e, &
                            z_temp_e )   !output

  !calculate heat flux divergence at cell centres.

  i_startblk = pt_patch%edges%start_blk(2,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk,nblks_e

    CALL get_indices_e(pt_patch, jb, i_startblk, nblks_e, &
                       i_startidx, i_endidx, 2)

    z_temp_flux_e(i_startidx:i_endidx,:,jb) =  &
    &  z_temp_e(i_startidx:i_endidx,:,jb) * pt_prog%vn(i_startidx:i_endidx,:,jb)
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  SELECT CASE (idiv_method)
  CASE(1)

    CALL div(z_temp_flux_e, pt_patch, pt_int_state, & ! input
             z_tmdiv, opt_rlstart=2 )     ! output

  CASE (2)

    CALL div_avg(z_temp_flux_e, pt_patch, pt_int_state, pt_int_state%c_bln_avg, &
                 z_tmdiv, opt_rlstart=2 )

  END SELECT


   ! For nested domains, tendencies are interpolated from the parent domain
   ! on a boundary zone with a width of grf_bdywidth_c for cells and
   ! grf_bdywidth_e for edges, respectively. These tendencies must not be
   ! overwritten.

   i_startblk = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk,nblks_c

     CALL get_indices_c(pt_patch, jb, i_startblk, nblks_c, &
                        i_startidx, i_endidx, grf_bdywidth_c+1)

     pt_tend_dyn%temp(i_startidx:i_endidx,:,jb) = &
 &   pt_tend_dyn%temp(i_startidx:i_endidx,:,jb)   &
 &          - z_tmdiv(i_startidx:i_endidx,:,jb)

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!-------------------------------------------------------------------------------
! 5. horizontal advection of momentum and the Coriolis force
!-------------------------------------------------------------------------------

    ! get kinetic energy and tangential wind
    CALL kin_vel_rot(pt_prog%vn, pt_patch, pt_int_state, pt_diag )

    ! get complete rotation term (abs. vort. times v)
    CALL lamb_rot( pt_patch, pt_int_state, pt_diag, &! in
                   pt_tend_dyn%vn,                  &! inout
                   opt_rlstart=grf_bdywidth_e+1)     ! optional

!-------------------------------------------------------------------------------
! 6. Compute the geopotential + kinetic energy and its horizontal gradient
!-------------------------------------------------------------------------------

!$OMP PARALLEL PRIVATE(i_startblk)
  i_startblk = pt_patch%cells%start_blk(2,1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk,nblks_c

     CALL get_indices_c(pt_patch, jb, i_startblk, nblks_c, &
                        i_startidx, i_endidx, 2)

! 6.1 If a temperature reference state is employed, compute the modified
!     boudary condition for geopotential at the earth's surface;
!     Otherwise get the values of surface geopotential at cell centers.

      z_geos_c(i_startidx:i_endidx,jb) =  &
        pt_ext_data%atm%topography_c(i_startidx:i_endidx,jb) * grav &
        - z_zbar_ic(i_startidx:i_endidx,nlevp1,jb)

! 6.2 Compute geopotential and its perturbation at other vertical levels.
!     First the perturbation:

      CALL geopot( z_tvp_c(:,:,jb), pt_diag%rdlnpr_c(:,:,jb),    &! in
                   pt_diag%rdalpha_c(:,:,jb),                    &! in
                   z_geos_c(:,jb), nproma, i_startidx, i_endidx, &! in
                   z_geop_mc(:,:,jb), z_geop_ic(:,:,jb)         ) ! inout

!  Then add the reference state to get the geopotential

!      IF (l_outputtime ) THEN
        pt_diag%geo_ic(i_startidx:i_endidx,2:nlevp1,jb) =  &
          z_geop_ic(i_startidx:i_endidx,2:nlevp1,jb) &
          + z_zbar_ic(i_startidx:i_endidx,2:nlevp1,jb)

        pt_diag%geo_mc(i_startidx:i_endidx,1:nlev,jb) =  &
          z_geop_mc(i_startidx:i_endidx,1:nlev,jb) &
          + z_zbar_mc(i_startidx:i_endidx,1:nlev,jb)
!     ENDIF

   ENDDO
!$OMP END DO

! 6.3 the gradient of geopotential and its perturbation at edges on full levels
!     add the kinetic energy (call only once the gradient)

   i_startblk = pt_patch%cells%start_blk(2,1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk,nblks_c

     CALL get_indices_c(pt_patch, jb, i_startblk, nblks_c, &
                        i_startidx, i_endidx, 2)

      z_e_mech_c(i_startidx:i_endidx,:,jb) = z_geop_mc(i_startidx:i_endidx,:,jb) &
                              + pt_diag%e_kin(i_startidx:i_endidx,:,jb)

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!-------------------------------------------------------------------------------
! 7. Compute the pressure gradient force
!-------------------------------------------------------------------------------
! 7.0 horizontal gradient of surface pressure at edges
!    The adiabatic heating term in the thermodynamic equation
!    is discretized accordingly to conserve the total energy.
!    (EM/DM scheme)

  ! 7-1.1 the first part of pressure gradient force: geop+kin gradient

    CALL grad_fd_norm( z_e_mech_c, pt_patch, z_pres_grad1, opt_rlstart=4 )

  ! average virtual temperature perturbation from cells to edges
  ! Note: by giving the vertical start index optionally, the comp. are
  ! done only for non pressure levels

   CALL cells2edges_scalar( z_tvp_c, pt_patch, pt_int_state%c_lin_e, &
                            z_tvp_e, nplvp1, opt_rlstart=4)

  ! The counterpart to be used in the thermodynamic equation:
  ! Even if the reference state is involved, the adiabatic heating
  ! is still calculated using the original variables.

   IF (l_outputtime) THEN
     CALL cells2edges_scalar( z_tv_c, pt_patch,  pt_int_state%c_lin_e, &
                              z_tv_e, nplvp1, opt_rlstart=4)
   ENDIF

  ! pressure gradient at edges

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jkp) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jk = nplvp1, nlev
         jkp = jk+1

         z_rdlnp_mc(1:nlen,jk,jb) = rd*z_lnp_ic(1:nlen,jkp,jb) - &
                                    pt_diag%rdalpha_c(1:nlen,jk ,jb)
      ENDDO
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   CALL grad_fd_norm( z_rdlnp_mc, pt_patch, z_rdlnp_me, nplvp1 )

  ! pressure gradient force, part 2, for momentum equation

!$OMP PARALLEL
   i_startblk = pt_patch%edges%start_blk(4,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk,nblks_e

     CALL get_indices_e(pt_patch, jb, i_startblk, nblks_e, &
                        i_startidx, i_endidx, 4)

     DO jk = nplvp1, nlev

     ! pressure gradient force, part 2, for momentum equation
      z_pres_grad2(i_startidx:i_endidx,jk,jb) = &
       -z_tvp_e(i_startidx:i_endidx,jk,jb)*z_rdlnp_me(i_startidx:i_endidx,jk,jb)
     ENDDO

      z_pres_grad2(i_startidx:i_endidx,1:nplev,jb) = 0._wp

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  ! pressure gradient force, part 2, for temperature equation
  ! (Only needed for omega output)
   IF (l_outputtime) THEN
!$OMP PARALLEL
     i_startblk = pt_patch%edges%start_blk(4,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk) ICON_OMP_DEFAULT_SCHEDULE
     DO jb = i_startblk,nblks_e

     CALL get_indices_e(pt_patch, jb, i_startblk, nblks_e, &
                        i_startidx, i_endidx, 4)

       DO jk = nplvp1, nlev

        ! pressure gradient force, part 2, for temperature equation
        z_pres_grad4(i_startidx:i_endidx,jk,jb) = &
         -z_tv_e (i_startidx:i_endidx,jk,jb)*z_rdlnp_me(i_startidx:i_endidx,jk,jb)
       ENDDO

      z_pres_grad4(i_startidx:i_endidx,1:nplev,jb) = 0._wp

     ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
   ENDIF

  ! 7-1.7 accumulate momentum tendency

   ! For nested domains, tendencies are interpolated from the parent domain
   ! on a boundary zone with a width of grf_bdywidth_c for cells and
   ! grf_bdywidth_e for edges, respectively. These tendencies must not be
   ! overwritten.
!$OMP PARALLEL PRIVATE(i_startblk)
   i_startblk = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk,nblks_e

     CALL get_indices_e(pt_patch, jb, i_startblk, nblks_e, &
                        i_startidx, i_endidx, grf_bdywidth_e+1)

     DO jk = 1, nlev
       DO je = i_startidx, i_endidx
         pt_tend_dyn%vn(je,jk,jb) =  pt_tend_dyn%vn(je,jk,jb) &
                                  &  - z_pres_grad1(je,jk,jb) &
                                  &  + z_pres_grad2(je,jk,jb)
       ENDDO
     ENDDO
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!-------------------------------------------------------------------------------
! 8. Compute the vertical velocity dp/dt
!-------------------------------------------------------------------------------
  ! 8.1 vn*[Rd*T/p*grad(p)]

! For theta advection, the following computations are needed only
! if omega is requested as an output field

IF (l_outputtime .AND. .NOT. output_mode%l_none) THEN

!$OMP PARALLEL PRIVATE(i_startblk)
   i_startblk = pt_patch%edges%start_blk(4,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk,nblks_e

     CALL get_indices_e(pt_patch, jb, i_startblk, nblks_e, &
                        i_startidx, i_endidx, 4)

      z_tmp_me(i_startidx:i_endidx,:,jb) = &
       - pt_diag%mass_flux_e(i_startidx:i_endidx,:,jb) * &
         z_pres_grad4(i_startidx:i_endidx,:,jb)

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   CALL edges2cells_scalar( z_tmp_me, pt_patch, pt_int_state%e_inn_c, &
                            z_tmp_mc, opt_rlstart=3 )

   IF ( idiv_method == 2 ) THEN

     z_aux = z_tmp_mc

     CALL cell_avg(z_aux, pt_patch, pt_int_state%c_bln_avg, z_tmp_mc, opt_rlstart=4)

   END IF

  ! 8.2 Rd*T/p *[ p-tendency + vertical-adv ]

!$OMP PARALLEL PRIVATE(i_startblk)
   i_startblk = pt_patch%cells%start_blk(3,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk,nblks_c

     CALL get_indices_c(pt_patch, jb, i_startblk, nblks_c, &
                        i_startidx, i_endidx, 3)

     DO jk = 1, nlev
       DO jc = i_startidx, i_endidx
      z_dpdt(jc,jk,jb) = pt_diag%rdelp_c(jc,jk,jb) * &
            ( z_tmp_mc(jc,jk,jb) - z_tv_c(jc,jk,jb) *  &
            ( pt_diag%rdlnpr_c (jc,jk,jb) * z_mdiv_int(jc,jk,jb) &
             +pt_diag%rdalpha_c(jc,jk,jb) * z_mdiv    (jc,jk,jb) ) )
       ENDDO
     ENDDO

   ENDDO
!$OMP END DO

   i_startblk = pt_patch%cells%start_blk(4,1)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = i_startblk,nblks_c

     CALL get_indices_c(pt_patch, jb, i_startblk, nblks_c, &
                        i_startidx, i_endidx, 4)

     DO jk = 1, nlev
       DO jc = i_startidx, i_endidx
         pt_diag%wpres_mc(jc,jk,jb) = z_dpdt(jc,jk,jb)*pt_diag%pres_mc(jc,jk,jb)/&
                                   rd/z_tv_c(jc,jk,jb)

       ENDDO
     ENDDO

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

ENDIF

 CALL sync_patch_array( SYNC_C, pt_patch, pt_tend_dyn%pres_sfc )
 CALL sync_patch_array( SYNC_C, pt_patch, pt_tend_dyn%temp     )
 CALL sync_patch_array( SYNC_E, pt_patch, pt_tend_dyn%vn       )
  
  IF (ltimer) CALL timer_stop(timer_dyn_theta)

END SUBROUTINE dyn_theta

!-------------------------------------------------------------------------

END MODULE m_dyn

