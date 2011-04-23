!>
!! This module contains subroutines for evaluating the right-hand side
!! of the primitive equations
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan (MPI-M, 2009-11)
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
MODULE mo_ha_dynamics


  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: rcpd, rd
  USE mo_model_domain,       ONLY: t_patch
  USE mo_ext_data,           ONLY: t_external_data
  USE mo_math_operators,     ONLY: grad_fd_norm, div, div_avg
  USE mo_dynamics_nml,       ONLY: idiv_method, lref_temp
  USE mo_run_nml,            ONLY: nproma, lshallow_water, nlev, nlevp1
  USE mo_icoham_dyn_types,   ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_interpolation,      ONLY: t_int_state, cell_avg, cells2edges_scalar, &
                                   edges2cells_scalar

  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_ha_dynamics_adv,    ONLY: temp_adv_vertical, temp_adv_horizontal, &
                                     vn_adv_vertical,   vn_adv_horizontal
  USE mo_ha_diag_util,       ONLY: update_diag_state
  USE mo_vertical_coord_table, ONLY: nplev, nplvp1, vct_b, alrrdic, rdlnp0i,    &
    &                                rdt0ral, t0icao

  USE mo_sync,                 ONLY: SYNC_C, SYNC_E, sync_patch_array

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  PUBLIC :: dyn_temp
  PUBLIC :: continuity
  PUBLIC :: energy_conversion_terms

CONTAINS
  !>
  !!
  !! @par Revision History
  !!
  !! First version by Hui Wan (MPI-M, 2009-11)
  !!
  SUBROUTINE dyn_temp( pt_patch, pt_int_state, pt_prog,    & ! input
                       pt_ext_data, pt_diag, pt_tend_dyn   ) ! input, output

  !! Arguments

  TYPE(t_patch),TARGET,    INTENT(IN)    :: pt_patch
  TYPE(t_int_state),TARGET,INTENT(IN)    :: pt_int_state
  TYPE(t_hydro_atm_prog),INTENT(IN)    :: pt_prog
  TYPE(t_external_data),   INTENT(IN)    :: pt_ext_data   !< external data

  TYPE(t_hydro_atm_diag),INTENT(INOUT) :: pt_diag
  TYPE(t_hydro_atm_prog),INTENT(INOUT) :: pt_tend_dyn

  !! Local variables

  REAL(wp) :: z_mdiv_int (nproma,nlevp1,pt_patch%nblks_c)
  REAL(wp) :: z_mdiv     (nproma,nlev  ,pt_patch%nblks_c)

  INTEGER :: nblks_c, nblks_e
  INTEGER :: jb, jbs, is, ie

! Dimension parameters related to refinement and MPI parallelisation

   nblks_c = pt_patch%nblks_int_c
   nblks_e = pt_patch%nblks_int_e

! Update the diagnostic state vector. This includes the calculation of
! pressure and related quantities, vorticity and divergence, u- and v-wind,
! virtual temperature and geopotential.

   CALL update_diag_state( pt_prog, pt_patch, pt_int_state, pt_ext_data, &
     &                     pt_diag )

! Diagnose the mass flux, eta-coordinate vertical velocity,
! and calculate surface pressure tendency

   CALL continuity( pt_prog%vn, pt_diag%delp_e,         &! in
                    pt_patch, pt_int_state, .TRUE.,     &! in
                    z_mdiv, z_mdiv_int,                 &! inout
                    pt_diag%mass_flux_e,                &! inout
                    pt_tend_dyn%pres_sfc,               &! inout
                    pt_diag%weta                    )    ! inout

! Initialize velocity and temperature tendencies in the interior of each patch

   jbs = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1)
      pt_tend_dyn%vn(is:ie,:,jb) = 0._wp
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

   IF (.NOT.lshallow_water) THEN
      jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
      DO jb = jbs,nblks_c
         CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)
         pt_tend_dyn%temp(is:ie,:,jb) = 0._wp
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
   ENDIF

! From now on, after calling each subroutine the individual contribution
! will be added to the tendency state.

   IF (.NOT.lshallow_water) THEN

! Vertical advection of momentum

     CALL vn_adv_vertical( pt_prog%vn, pt_diag%weta, pt_diag%delp_e, &! in
                           pt_patch,   pt_int_state,                 &! in
                           pt_tend_dyn%vn                          )  ! inout

! Vertical and horizontal advection of temperature

     CALL temp_adv_vertical( pt_prog%temp, pt_diag%weta,    &! in
                             pt_diag%rdelp_c, pt_patch,     &! in
                             pt_tend_dyn%temp            )   ! inout

     CALL temp_adv_horizontal( pt_prog%temp, pt_diag%mass_flux_e,  &! in
                               pt_diag%rdelp_c, z_mdiv,            &! in
                               pt_patch, pt_int_state,             &! in
                               pt_tend_dyn%temp                  )  ! inout
   ENDIF

! Horizontal advection and Coriolis force

   CALL sync_patch_array( SYNC_E, pt_patch, pt_diag%mass_flux_e )
   CALL vn_adv_horizontal( pt_prog%vn,                          &! in
                           pt_diag%rel_vort, pt_diag%rel_vort_e,&! in
                           pt_diag%mass_flux_e, pt_diag%delp_c, &! in
                           pt_patch, pt_int_state,              &! in
                           pt_tend_dyn%vn,                      &! inout
                           pt_diag%e_kin, pt_diag%vt,           &! inout
                           pt_diag%delp_v, pt_diag%delp_e,      &! inout
                           opt_rlstart=grf_bdywidth_e+1 )        ! for nesting

! Pressure gradient force and adiabatic heating

   CALL energy_conversion_terms( pt_diag, pt_patch, pt_int_state, &! inout,in,in
                                 z_mdiv, z_mdiv_int,              &! in
                                 pt_tend_dyn%temp,                &! inout
                                 pt_tend_dyn%vn                  ) ! inout

! Synchronize tendencies


   CALL sync_patch_array( SYNC_C, pt_patch, pt_tend_dyn%pres_sfc )
   CALL sync_patch_array( SYNC_C, pt_patch, pt_tend_dyn%temp     )
   CALL sync_patch_array( SYNC_E, pt_patch, pt_tend_dyn%vn       )

  END SUBROUTINE dyn_temp
  !----------------------

  !>
  !! Spatially discretized continuity equation
  !!
  !! Purpose:
  !! Vertically integrate mass divergence to compute surface pressure tendency;
  !! Compute the vertical velocity eta-dot at layer interfaces.)
  !!
  !! @par Revision History
  !! Separated from subroutine dyn by Hui Wan (MPI-M, 2009-11-17)
  !!
  SUBROUTINE continuity( p_vn, p_delp_e, pt_patch, pt_int_state, ldiag_weta, &
                         p_mdiv, p_mdiv_int, p_mflux, p_ddt_psfc, p_weta )

  IMPLICIT NONE

  !! Arguments

  REAL(wp),       INTENT(in) :: p_vn    (:,:,:) !< normal velocity
  REAL(wp),       INTENT(in) :: p_delp_e(:,:,:) !< layer thickness at edges
  TYPE(t_patch), TARGET, INTENT(in) :: pt_patch   !< grid information
  TYPE(t_int_state),INTENT(in) :: pt_int_state    !< interpolation coefficients

  LOGICAL,INTENT(in) :: ldiag_weta              !< if .true., diagnoise weta

  REAL(wp),INTENT(inout) :: p_mdiv    (:,:,:)   !< mass divergence
  REAL(wp),INTENT(inout) :: p_mdiv_int(:,:,:)   !< mass divergence
                                                !< vertically integrated

  REAL(wp),INTENT(inout) :: p_mflux(:,:,:)  !< mass flux at edges

  REAL(wp),INTENT(inout) :: p_ddt_psfc(:,:) !< tendency of surface pressure

  REAL(wp),INTENT(inout),OPTIONAL :: p_weta (SIZE(p_mdiv_int,1), &!< vertical
                                             SIZE(p_mdiv_int,2), &!< velocity
                                             SIZE(p_mdiv_int,3) ) !< rho*eta-dot

  !! Local variables

  INTEGER  :: nblks_e, nblks_c
  INTEGER  :: jb, jbs, is,ie, jk,jkp

! Dimension parameters

  nblks_e  = pt_patch%nblks_int_e
  nblks_c  = pt_patch%nblks_int_c

! Divergence of mass flux at full levels

!$OMP PARALLEL
   jbs = pt_patch%edges%start_blk(2,1)
!$OMP DO PRIVATE(jb,is,ie)
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, 2)
      p_mflux(is:ie,:,jb) = p_delp_e(is:ie,:,jb)*p_vn(is:ie,:,jb)
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

   CALL sync_patch_array(SYNC_E, pt_patch, p_mflux)

   SELECT CASE(idiv_method)
   CASE(1)

       CALL div( p_mflux, pt_patch, pt_int_state, p_mdiv, opt_rlstart=2 )

   CASE(2)

       CALL div_avg( p_mflux, pt_patch, pt_int_state, &
                     pt_int_state%c_bln_avg, p_mdiv, opt_rlstart=2 )

   END SELECT

! Vertically integrated mass divergence at layer interfaces

   jbs = pt_patch%cells%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk,jkp)
   DO jb = jbs,nblks_c
      CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, 2)

      p_mdiv_int(is:ie,1,jb) = 0.0_wp
      DO jk = 1,nlev
         jkp = jk + 1
         p_mdiv_int(is:ie,jkp,jb) = p_mdiv_int(is:ie,jk,jb)+p_mdiv(is:ie,jk,jb)
      ENDDO
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

! Vertical velocity at the interfaces (half-levels)

   IF ((.NOT.lshallow_water).AND.ldiag_weta) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk)
      DO jb = jbs,nblks_c
        CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, 2)

        p_weta(is:ie,1,jb)      = 0._wp   ! upper boundary
        p_weta(is:ie,nlevp1,jb) = 0._wp   ! lower boundary

        DO jk = 2,nlev
           p_weta(is:ie,jk,jb)  = -p_mdiv_int(is:ie,jk,jb) + &
                                   p_mdiv_int(is:ie,nlevp1,jb)*vct_b(jk)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
   ENDIF

! Tendency of surface pressure
! (For nested domains, tendencies are interpolated from the parent domain
! on a boundary zone with a width of grf_bdywidth_c for cells and
! grf_bdywidth_e for edges, respectively. These tendencies must not be
! overwritten.)

   jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
   DO jb = jbs,nblks_c
      CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)
      p_ddt_psfc(is:ie,jb) = -p_mdiv_int(is:ie,nlevp1,jb)
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE continuity


  !>
  !! energy_conversion_terms
  !!
  !! Purpose:
  !!
  !! @par Revision History
  !! Separated from subroutine dyn and rewritten by Hui Wan (MPI-M, 2009-11-18)
  !!
  SUBROUTINE energy_conversion_terms( pt_diag, pt_patch, pt_int_state, &
                                      p_mdiv, p_mdiv_int,              &
                                      p_ddt_temp, p_ddt_vn,            &
                                      opt_lseparate, opt_ddt_temp_fast )

  !! Arguments

  TYPE(t_hydro_atm_diag),INTENT(INOUT) :: pt_diag
  TYPE(t_patch),TARGET,    INTENT(IN)    :: pt_patch
  TYPE(t_int_state),       INTENT(IN)    :: pt_int_state

  REAL(wp),INTENT(IN)    :: p_mdiv     (:,:,:)
  REAL(wp),INTENT(IN)    :: p_mdiv_int (:,:,:)
  REAL(wp),INTENT(INOUT) :: p_ddt_temp (:,:,:)
  REAL(wp),INTENT(INOUT) :: p_ddt_vn   (:,:,:)

  LOGICAL, INTENT(IN),   OPTIONAL :: opt_lseparate
  REAL(wp),INTENT(INOUT),OPTIONAL :: opt_ddt_temp_fast(:,:,:)

  !! Local variables

  LOGICAL  :: lseparate
  INTEGER  :: jk,jkp
  INTEGER  :: jb,jbs,is,ie
  INTEGER  :: nblks_e,nblks_c
  REAL(wp) :: z2d (nproma,nlev)

  REAL(wp) :: z_geo_mc( nproma,nlev,pt_patch%nblks_c )
  REAL(wp) :: z_tv_c  ( nproma,nlev,pt_patch%nblks_c )
  REAL(wp) :: z_tv_e  ( nproma,nlev,pt_patch%nblks_e )
  REAL(wp) :: z_tvp_c ( nproma,nlev,pt_patch%nblks_c )
  REAL(wp) :: z_tvp_e ( nproma,nlev,pt_patch%nblks_e )
  REAL(wp) :: z_tmp_e ( nproma,nlev,pt_patch%nblks_e )
  REAL(wp) :: z_tmp_c ( nproma,nlev,pt_patch%nblks_c )
  REAL(wp) :: z_rlp_c ( nproma,nlev,pt_patch%nblks_c )
  REAL(wp) :: z_fast  ( nproma,nlev,pt_patch%nblks_c )

! Check optional input

  IF (PRESENT(opt_lseparate)) THEN
     lseparate = opt_lseparate
  ELSE
     lseparate = .FALSE.
  ENDIF

! Dimension parameters

  nblks_c = pt_patch%nblks_int_c
  nblks_e = pt_patch%nblks_int_e

!=====================================================================
! Horizontal gradient of geopotential
!=====================================================================
IF (.NOT.lshallow_water) THEN
! For the hydrostatic model,
! if the use of a reference state is desired (to reduce the numerical
! error near steep topography), we need to first construct the reference
! state and calculate the temperature and geopotential perturbation.
! The reference state we use here is the ICAO(1964) standard atmosphere.
! See MPI-M Report 349, p20.

  IF (lref_temp) THEN

    jbs = pt_patch%cells%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,z2d)
    DO jb = jbs,nblks_c
       CALL get_indices_c( pt_patch, jb,jbs,nblks_c, is,ie, 2)

       z2d(is:ie,:) = rd*LOG(pt_diag%pres_mc(is:ie,:,jb)) - rdlnp0i
       z2d(is:ie,:) = EXP( alrrdic*z2d(is:ie,:) )

       z_tv_c  (is:ie,:,jb) = pt_diag%tempv (is:ie,:,jb)
       z_tvp_c (is:ie,:,jb) = pt_diag%tempv (is:ie,:,jb) -  t0icao*z2d(is:ie,:)
       z_geo_mc(is:ie,:,jb) = pt_diag%geo_mc(is:ie,:,jb) + rdt0ral*z2d(is:ie,:)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  ELSE !Do not use reference state

    jbs = pt_patch%cells%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
    DO jb = jbs,nblks_c
       CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, 2)
       z_tv_c  (is:ie,:,jb) = pt_diag%tempv (is:ie,:,jb)
       z_tvp_c (is:ie,:,jb) = pt_diag%tempv (is:ie,:,jb)
       z_geo_mc(is:ie,:,jb) = pt_diag%geo_mc(is:ie,:,jb)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDIF

ELSE !shallow water model
!$OMP PARALLEL
!$OMP WORKSHARE
      z_geo_mc(:,:,:) = pt_diag%geo_mc(:,:,:)
!$OMP END WORKSHARE
!$OMP END PARALLEL
ENDIF

!--------------------------------------------------------------------------
! Calculate the gradient of geopotential, and accumulate velocity tendency
!--------------------------------------------------------------------------

  CALL grad_fd_norm( z_geo_mc, pt_patch, z_tmp_e, opt_rlstart=4 )

  jbs = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1)
      p_ddt_vn(is:ie,:,jb) = p_ddt_vn(is:ie,:,jb) - z_tmp_e(is:ie,:,jb)
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

!=====================================================================
! The other part of the pressure gradient force and the thermodynamic
! equation only exist in the hydrostatic model
!=====================================================================
IF (.NOT.lshallow_water) THEN

! Average virtual temperature perturbation from cells to edges
! Note: by giving the vertical start index optionally, the computation
! is done only for non-pressure levels

   CALL cells2edges_scalar( z_tvp_c, pt_patch, pt_int_state%c_lin_e, &
                            z_tvp_e, nplvp1, opt_rlstart=4)

! The counterpart to be used in the thermodynamic equation:
! Even if the reference state is involved, the adiabatic heating
! is still calculated using the original variables.

   CALL cells2edges_scalar( z_tv_c, pt_patch, pt_int_state%c_lin_e, &
                            z_tv_e, nplvp1, opt_rlstart=4)

! Pressure gradient at edges

   jbs = pt_patch%cells%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk,jkp)
   DO jb = jbs, nblks_c
      CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, 2)
      DO jk = nplvp1, nlev
         jkp = jk+1
         z_rlp_c(is:ie,jk,jb) = rd*pt_diag%lnp_ic(is:ie,jkp,jb) &
                              - pt_diag%rdalpha_c(is:ie,jk ,jb)
      ENDDO
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

   CALL grad_fd_norm( z_rlp_c, pt_patch, z_tmp_e, nplvp1, opt_rlstart=4 )

!--------------------------------------------------------------------
! Accumulate velocity tendency
!--------------------------------------------------------------------

!$OMP PARALLEL
   jbs = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP DO PRIVATE(jb,is,ie,jk)
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1)

      DO jk = nplvp1, nlev
         p_ddt_vn(is:ie,jk,jb) = p_ddt_vn(is:ie,jk,jb) &
                                -z_tvp_e(is:ie,jk,jb)*z_tmp_e(is:ie,jk,jb)
      ENDDO
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

!--------------------------------------------------------------------
! Adiabatic heating terms in the thermodynamic equation
!--------------------------------------------------------------------
! Part I: vn*[Rd*T/p*grad(p)]. Use the grad(lnp) calculated above.

   jbs = pt_patch%edges%start_blk(4,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk)
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, 4)

      DO jk = 1, nplev
         z_tmp_e(is:ie,jk,jb) = 0._wp  !Pressure gradient vanishes on p-levels
      ENDDO

      DO jk = nplvp1, nlev
         z_tmp_e(is:ie,jk,jb) =  pt_diag%mass_flux_e(is:ie,jk,jb) &
                                *z_tv_e (is:ie,jk,jb)             &
                                *z_tmp_e(is:ie,jk,jb)
      ENDDO
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

! Interpolation from edges to cell centers

   CALL edges2cells_scalar( z_tmp_e, pt_patch, pt_int_state%e_inn_c, &
                            z_tmp_c, opt_rlstart=3 )

! Add smoothing if desired. (Only when using the original gradient operator.
! Otherwise the model will become unstable.

   SELECT CASE (idiv_method)

   CASE(2)

!$OMP PARALLEL
!$OMP WORKSHARE
       z_rlp_c = z_tmp_c
!$OMP END WORKSHARE
!$OMP END PARALLEL
       CALL cell_avg( z_rlp_c, pt_patch, pt_int_state%c_bln_avg, &
                      z_tmp_c, opt_rlstart=4)

   END SELECT

! Divide by the pseudo-density

   jbs = pt_patch%cells%start_blk(3,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
   DO jb = jbs,nblks_c
      CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, 3)
      z_tmp_c(is:ie,:,jb) = z_tmp_c(is:ie,:,jb) * pt_diag%rdelp_c(is:ie,:,jb)
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

!-------------------------------------------------------------------
! Part II: Rd*T/p *[ p-tendency + vertical-adv ]

   jbs = pt_patch%cells%start_blk(3,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
   DO jb = jbs,nblks_c
      CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, 3)

      z_fast(is:ie,:,jb) = -pt_diag%rdelp_c(is:ie,:,jb)*z_tv_c(is:ie,:,jb)      &
                 *( pt_diag%rdlnpr_c (is:ie,:,jb) * p_mdiv_int(is:ie,1:nlev,jb) &
                   +pt_diag%rdalpha_c(is:ie,:,jb) * p_mdiv    (is:ie,:,jb)     )
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

!-------------------------------------------------------------------
! Accumulate temperature tendency

   jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP PARALLEL
  IF (lseparate) THEN !Provide the fast and slow components separately

!$OMP DO PRIVATE(jb,is,ie)
   DO jb = jbs,nblks_c
     CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)
     p_ddt_temp(is:ie,:,jb) = p_ddt_temp(is:ie,:,jb) + z_tmp_c(is:ie,:,jb)*rcpd
     opt_ddt_temp_fast(is:ie,:,jb) = z_fast(is:ie,:,jb)*rcpd
   ENDDO
!$OMP END DO

  ELSE !Add both fast and slow components to the tendency state

!$OMP DO PRIVATE(jb,is,ie)
   DO jb = jbs,nblks_c
     CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)
     p_ddt_temp(is:ie,:,jb) = p_ddt_temp(is:ie,:,jb) &
                            + ( z_tmp_c(is:ie,:,jb) + z_fast(is:ie,:,jb) )*rcpd
   ENDDO
!$OMP END DO

  ENDIF
!$OMP END PARALLEL

ENDIF !shallow_water vs hydrostatic
!------------------------------------------
  END SUBROUTINE energy_conversion_terms

END MODULE mo_ha_dynamics
