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
MODULE mo_ha_dynamics_adv

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message,finish
  USE mo_model_domain,       ONLY: t_patch
  USE mo_math_operators,     ONLY: grad_fd_norm, div, div_avg
  USE mo_dynamics_nml,       ONLY: idiv_method
  USE mo_io_nml,             ONLY: l_outputtime, lwrite_omega, l_diagtime
  USE mo_run_nml,            ONLY: nproma, nlev, nlevp1, lshallow_water,   &
                                   i_cell_type
  USE mo_interpolation,      ONLY: t_int_state, rbf_vec_interpol_edge,     &
                                   cells2edges_scalar, edges2cells_scalar, &
                                   verts2edges_scalar, cells2verts_scalar, &
                                   edges2verts_scalar, verts2cells_scalar, &
                                   edges2edges_scalar, i_cori_method,      &
                                   sick_a, sick_o, l_corner_vort
  USE mo_vertical_coord_table, ONLY: rdelpr, nlevm1, nplev, nplvp1
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_impl_constants,     ONLY: min_rledge
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_sync,               ONLY: SYNC_C, SYNC_V, sync_patch_array

  IMPLICIT NONE

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  PUBLIC

CONTAINS

  !>
  !! Purpose: calculate the tendency of normal wind induced by vertical
  !! advection
  !!
  !! General comment:
  !! For nested domains, tendencies are interpolated from the parent domain
  !! on a boundary zone with a width of grf_bdywidth_c for cells and
  !! grf_bdywidth_e for edges, respectively. These tendencies must not be
  !! overwritten.
  !!
  !! Reference:
  !! Simmons and Burridge (1981, Mon. Wea. Rev.)
  !! Roeckner et al. (2005, MPI Report 349)
  !!
  !! @par Revision History
  !! Separated from subroutine dyn and re-written by Hui Wan (MPI-M, 2009-11-17)
  !!
  SUBROUTINE vn_adv_vertical( p_vn, p_weta_c, p_delp_e, &
                              pt_patch, pt_int_state,   &
                              p_ddt_vn )
  !! Arguments

  REAL(wp),INTENT(IN)    :: p_vn    (:,:,:)   !<  normal velocity
  REAL(wp),INTENT(INOUT) :: p_weta_c(:,:,:)   !<  cell-based vertical velocity
  REAL(wp),INTENT(IN)    :: p_delp_e(:,:,:)   !<  layer thickness at edges

  TYPE(t_patch),TARGET,INTENT(IN) :: pt_patch   !<  grid information
  TYPE(t_int_state),INTENT(IN) :: pt_int_state  !<  interpolation coefficients

  REAL(wp), INTENT(INOUT) :: p_ddt_vn(:,:,:)  !<  tendency of normal velocity

  !! Local variables

  INTEGER  :: nblks_e
  INTEGER  :: jb,jbs,is,ie,jk,jkp

  REAL(wp) :: z_weta_e( SIZE(p_vn,1),SIZE(p_weta_c,2),SIZE(p_vn,3) )
  REAL(wp) :: z_tmp_e ( SIZE(p_vn,1),SIZE(p_vn,    2),SIZE(p_vn,3) )

! Dimension parameter

   nblks_e  = pt_patch%nblks_int_e

! Interpolate vertical velocity (rho*eta-dot) from cell centers to edges

   CALL sync_patch_array( SYNC_C, pt_patch, p_weta_c )

   CALL cells2edges_scalar( p_weta_c, pt_patch, pt_int_state%c_lin_e, &! in
                            z_weta_e, opt_rlstart=4 )    ! out, optional in

! Tendency of velocity due to vertical advection.

!$OMP PARALLEL PRIVATE(jbs)
   jbs = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP DO PRIVATE(jb,is,ie,jk,jkp)
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1)

      z_tmp_e(is:ie,1,jb) = 0._wp

      DO jk = 1, nlevm1
         jkp = jk + 1
         z_tmp_e(is:ie,jkp,jb) = z_weta_e(is:ie,jkp,jb)  &
                                 *( p_vn(is:ie,jk,jb)- p_vn(is:ie,jkp,jb) )
         z_tmp_e(is:ie,jk ,jb) = z_tmp_e(is:ie,jk,jb) + z_tmp_e(is:ie,jkp,jb)
      ENDDO

     !Pressure levels
      DO jk = 1,nplev
        p_ddt_vn(is:ie,jk,jb) = p_ddt_vn(is:ie,jk,jb) &
                              + z_tmp_e(is:ie,jk,jb)*rdelpr(jk)*0.5_wp
      ENDDO

     !Sigma and transition levels
      DO jk = nplvp1,nlev
        p_ddt_vn(is:ie,jk,jb) = p_ddt_vn(is:ie,jk,jb) &
                              + z_tmp_e(is:ie,jk,jb)/p_delp_e(is:ie,jk,jb) &
                               *0.5_wp
      ENDDO
   ENDDO
!$OMP END DO
!$OMP END PARALLEL
  END SUBROUTINE vn_adv_vertical

  !>
  !! Vertical advection of temperature
  !!
  !! For nested domains, tendencies are interpolated from the parent domain
  !! on a boundary zone with a width of grf_bdywidth_c for cells and
  !! grf_bdywidth_e for edges, respectively. These tendencies must not be
  !! overwritten.
  !!
  !! Reference:
  !! Simmons and Burridge (1981, Mon. Wea. Rev.)
  !! Roeckner et al. (2005, MPI Report 349)
  !!
  !! @par Revision History
  !! Separation from subroutine dyn and rewriting by Hui Wan (MPI-M, 2009-11-18)
  !!
  SUBROUTINE temp_adv_vertical( p_temp, p_weta, p_rdelp, pt_patch, p_ddt_temp )

  !! Arguments

  REAL(wp),INTENT(IN)    :: p_temp    (:,:,:) !< temperature
  REAL(wp),INTENT(IN)    :: p_weta    (:,:,:) !< vertical velocity rho*eta-dot
  REAL(wp),INTENT(IN)    :: p_rdelp   (:,:,:) !< 1./layer_thickness (delp)
  REAL(wp),INTENT(INOUT) :: p_ddt_temp(:,:,:) !< temperature tendency
  TYPE(t_patch),INTENT(IN) :: pt_patch          !< grid information

  !! Local variables

  INTEGER  :: jb,jbs,is,ie,jk,jkp
  INTEGER  :: nblks_c
  REAL(wp) :: z_tmp_c (SIZE(p_temp,1),SIZE(p_temp,2),SIZE(p_temp,3))

! Dimension parameter

  nblks_c = pt_patch%nblks_int_c

! Vertical advection

!$OMP PARALLEL  PRIVATE(jbs)
  jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP DO PRIVATE(jb,is,ie,jk,jkp)
  DO jb = jbs,nblks_c
     CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)

     z_tmp_c(is:ie,1,jb) = 0._wp

     DO jk = 1, nlevm1
        jkp = jk + 1
        z_tmp_c(is:ie,jkp,jb) = p_weta(is:ie,jkp,jb)  &
                               *( p_temp(is:ie,jk, jb) - p_temp(is:ie,jkp,jb) )
        z_tmp_c(is:ie,jk ,jb) = z_tmp_c(is:ie,jk,jb) + z_tmp_c(is:ie,jkp,jb)
     ENDDO

     p_ddt_temp(is:ie,1:nlev,jb) = p_ddt_temp(is:ie,1:nlev,jb)  &
                                 + z_tmp_c(is:ie,1:nlev,jb)     &
                                  *p_rdelp(is:ie,1:nlev,jb)*0.5_wp
   ENDDO
!$OMP END DO
!$OMP END PARALLEL
  END SUBROUTINE temp_adv_vertical

  !>
  !! Horizontal advection of temperature
  !!
  !! Energy conserving formula [ div(T*delp*V) - T*div(delp*V) ]/delp.
  !! For nested domains, tendencies are interpolated from the parent domain
  !! on a boundary zone with a width of grf_bdywidth_c for cells and
  !! grf_bdywidth_e for edges, respectively. These tendencies must not be
  !! overwritten.
  !!
  !! @par Revision History
  !! Separation from subroutine dyn and rewriting by Hui Wan (MPI-M, 2009-11-18)
  !!
  SUBROUTINE temp_adv_horizontal( p_temp, p_mflux_e, p_rdelp_c, p_mdiv, &
                                  pt_patch, pt_int_state, p_ddt_temp )
  !! Arguments

  REAL(wp),INTENT(IN)    :: p_temp     (:,:,:)
  REAL(wp),INTENT(IN)    :: p_mflux_e  (:,:,:)
  REAL(wp),INTENT(IN)    :: p_rdelp_c  (:,:,:)
  REAL(wp),INTENT(IN)    :: p_mdiv     (:,:,:)
  REAL(wp),INTENT(INOUT) :: p_ddt_temp (:,:,:)

  TYPE(t_patch),TARGET,INTENT(IN) :: pt_patch
  TYPE(t_int_state),INTENT(IN) :: pt_int_state

  !! Local variables

  INTEGER :: jb,jbs,is,ie,nblks_e,nblks_c

  REAL(wp),DIMENSION(SIZE(p_mflux_e,1),SIZE(p_mflux_e,2),SIZE(p_mflux_e,3)) &
          :: z_temp_e, z_flux_e, z_tmdiv

! Dimension parameters
  nblks_e = pt_patch%nblks_int_e
  nblks_c = pt_patch%nblks_int_c

! Interpolate temperature from cells to edges
  CALL cells2edges_scalar( p_temp, pt_patch, pt_int_state%c_lin_e, z_temp_e )

! Calculate heat divergence at cell centres

   jbs = pt_patch%edges%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, 2)
      z_flux_e(is:ie,:,jb) = z_temp_e(is:ie,:,jb)*p_mflux_e(is:ie,:,jb)
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

  SELECT CASE(idiv_method)
  CASE(1)
    CALL div( z_flux_e, pt_patch, pt_int_state, z_tmdiv, opt_rlstart=2 )
  CASE(2)
    CALL div_avg( z_flux_e, pt_patch, pt_int_state, pt_int_state%c_bln_avg, &
                  z_tmdiv, opt_rlstart=2 )
  END SELECT

! Temperature tendency due to horizontal advection

   jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
   DO jb = jbs,nblks_c
      CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)

      p_ddt_temp(is:ie,:,jb) = p_ddt_temp(is:ie,:,jb)     &
                             - p_rdelp_c(is:ie,:,jb)      &
                              *( z_tmdiv(is:ie,:,jb)      &
                                -p_mdiv(is:ie,:,jb)*p_temp(is:ie,:,jb) )

   ENDDO
!$OMP END DO
!$OMP END PARALLEL
  END SUBROUTINE temp_adv_horizontal


  !>
  !! Horizontal momentum advection and Coriolis force
  !! The algorithm on triangular grids is different from that on the
  !! hexagonal/pentagonal grids.
  !!
  !! @par Revision History
  !! Separated from mo_nonlinear_advection and re-written
  !! by Hui Wan (MPI-M, 2009-11-20)
  !!
  SUBROUTINE vn_adv_horizontal( p_vn, p_vort, p_vort_e, p_mflux_e, p_delp_c,   &
                                pt_patch, pt_int,                    &
                                p_ddt_vn, p_kin, p_vt,               &
                                p_delp_v, p_delp_e,                  &
                                opt_rlstart, opt_rlend              )
  !! Arguments
  IMPLICIT NONE

  REAL(wp),INTENT(IN) :: p_vn     (:,:,:)  !< normal velocity at edges
  REAL(wp),INTENT(IN) :: p_vort   (:,:,:)  !< relative vorticity at dual centers
  REAL(wp),INTENT(IN) :: p_vort_e (:,:,:)  !< relative vorticity at edges
  REAL(wp),INTENT(IN) :: p_mflux_e(:,:,:)  !< mass flux at edges
  REAL(wp),INTENT(IN) :: p_delp_c (:,:,:)  !< pseudo-dencity at cell centers

  TYPE(t_patch),    TARGET,INTENT(IN) :: pt_patch
  TYPE(t_int_state),TARGET,INTENT(IN) :: pt_int

  REAL(wp),INTENT(INOUT) :: p_ddt_vn (:,:,:) !< tendency of normal wind
  REAL(wp),INTENT(INOUT) :: p_kin    (:,:,:) !< kinetic energy at cell centers
  REAL(wp),INTENT(INOUT) :: p_vt     (:,:,:) !< tangential wind at edges
  REAL(wp),INTENT(INOUT) :: p_delp_v (:,:,:) !< pseudo-dencity at dual grid
  REAL(wp),INTENT(INOUT) :: p_delp_e (:,:,:) !< pseudo-dencity at edges

  INTEGER,INTENT(IN),OPTIONAL :: opt_rlstart, opt_rlend
                                 !< start and end values of refin_ctrl flag
  !! Local variables

  INTEGER  :: jb,jbs,jbe,is,ie,je,jk,ji
  INTEGER  :: nblks_e, npromz_e, i_nchdom, nlen, nincr, nblks_c, npromz_c
  INTEGER  :: rl_start, rl_end

  REAL(wp) :: z_tmp_e    ( nproma, nlev, pt_patch%nblks_e )
  REAL(wp) :: z_tmp_c    ( nproma, nlev, pt_patch%nblks_c )
  REAL(wp) :: z_tmp_v    ( nproma, nlev, pt_patch%nblks_v )
  REAL(wp) :: z_vort_e   ( nproma, nlev, pt_patch%nblks_e )
  REAL(wp) :: z_potvort_r( nproma, nlev, pt_patch%nblks_e )
  REAL(wp) :: z_potvort_c( nproma, nlev, pt_patch%nblks_c )
  REAL(wp) :: z_potvort_v( nproma, nlev, pt_patch%nblks_v )
  REAL(wp) :: z_u_e      ( nproma, nlev, pt_patch%nblks_e )
  REAL(wp) :: z_v_e      ( nproma, nlev, pt_patch%nblks_e )
  REAL(wp) :: z_u_c      ( nproma, nlev, pt_patch%nblks_c )
  REAL(wp) :: z_v_c      ( nproma, nlev, pt_patch%nblks_c )

  INTEGER, DIMENSION(:,:,:),   POINTER :: ieidx, ieblk, icidx, icblk, ividx, ivblk

! Dimension parameters

  nblks_e  = pt_patch%nblks_int_e
  npromz_e = pt_patch%npromz_int_e
  nblks_c  = pt_patch%nblks_int_c
  npromz_c = pt_patch%npromz_int_c

  i_nchdom = MAX(1,pt_patch%n_childdom)

! For diagnosing potential enstrophy in the
! shallow water model, interpolate the pseudo-dencity from triangle centers
! to vertices.
  IF (l_diagtime .AND. lshallow_water) THEN
    CALL cells2verts_scalar(p_delp_c, pt_patch, &
                            pt_int%cells_aw_verts, p_delp_v)
  ENDIF

SELECT CASE (i_cell_type)
CASE(3)
!=====================================================================
! Triangular model
!=====================================================================
! Kinetic energy and tangential velocity
!-----------------------------------------

! Reconstruct the tangential wind at edges from the normal velocity;
! Define the specific kinetic energy at edges as ke = 0.5*( vn*vn + vt*vt )
! The cell-based value needed for computing ke gradient is then obtained
! by interpolation.

  CALL rbf_vec_interpol_edge( p_vn, pt_patch, pt_int, p_vt)

  jbs = pt_patch%edges%start_blk(2,1) !for rbf_vec_dim_edge==4
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
  DO jb = 1,nblks_e
     CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, 2)
     z_tmp_e(is:ie,:,jb) = 0.5_wp*( p_vt(is:ie,:,jb)*p_vt(is:ie,:,jb) &
                                   +p_vn(is:ie,:,jb)*p_vn(is:ie,:,jb) )
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  CALL edges2cells_scalar( z_tmp_e, pt_patch, pt_int%e_bln_c_s,  &
                           p_kin, opt_rlstart=2 )

!--------------------------
! Absolute vorticity flux
!--------------------------
! Dimension parameters

  IF ( PRESENT(opt_rlstart) ) THEN
    IF ((opt_rlstart >= 1) .AND. (opt_rlstart <= 2)) THEN
      CALL finish ('mo_ha_dynamics_adv:vn_adv_horizontal',  &
                   'opt_rlstart must not be 1 or 2')
    ENDIF
    rl_start = opt_rlstart
  ELSE
    rl_start = 3
  END IF

  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rledge
  END IF

! Compute the absolute vorticity flux

  CALL verts2edges_scalar( p_vort, pt_patch, pt_int%v_1o2_e, &
                           z_vort_e, opt_rlstart=3)

  jbs = pt_patch%edges%start_blk(rl_start,1)
  jbe = pt_patch%edges%end_blk(rl_end,i_nchdom)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk)
  DO jb = jbs,jbe
     CALL get_indices_e(pt_patch, jb,jbs,jbe, is,ie, rl_start,rl_end)

     DO jk = 1,nlev
        p_ddt_vn(is:ie,jk,jb) = p_ddt_vn(is:ie,jk,jb)             &
                              - p_vt(is:ie,jk,jb)                 &
                               *( z_vort_e(is:ie,jk,jb)           &
                                 +pt_patch%edges%f_e(is:ie,jb) )
     ENDDO
     !!HW: minus sign here because the tangential direction in the code is
     !!    in fact the opposite from what is written in BR05.
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

CASE(6)
!=====================================================================
! Hexagonal model
!=====================================================================
! Tangential velocity is not needed in the hexagonal model.
! Kinetic energy is defined from the normal component at edges,
! and then interpolated to cell centers.
!------------------------------------------------------------------
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen)
  DO jb = 1, nblks_e

     IF (jb /= nblks_e) THEN
       nlen = nproma
     ELSE
       nlen = npromz_e
     ENDIF

     z_tmp_e(1:nlen,:,jb) = 0.5_wp * p_vn(1:nlen,:,jb) * p_vn(1:nlen,:,jb)
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  CALL edges2cells_scalar( z_tmp_e, pt_patch, pt_int%e_inn_c, &! in
                           p_kin, opt_rlstart=2              ) ! out,in

  IF (i_cori_method >= 2 )THEN

    CALL edges2verts_scalar( z_tmp_e, pt_patch, pt_int%e_inn_v, &! in
                             z_tmp_v                             ) ! out,in
    CALL sync_patch_array(SYNC_V,pt_patch,z_tmp_v)
    CALL verts2cells_scalar( z_tmp_v, pt_patch, pt_int%verts_aw_cells, &! in
                             z_tmp_c                             ) ! out,in
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      DO jk = 1,nlev
        p_kin(1:nlen,jk,jb) = sick_o*p_kin(1:nlen,jk,jb)+sick_a*z_tmp_c(1:nlen,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  ENDIF

!----------------------------------
! Absolute vorticity flux
!----------------------------------
! Horizontal interpolation of thickness to vertex

  ! Absolute potential vorticity at rhombi

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk)
  DO jb = 1, nblks_e
    IF (jb /= nblks_e) THEN
      nlen = nproma
    ELSE
      nlen = npromz_e
    ENDIF
    DO jk = 1,nlev
      z_potvort_r(1:nlen,jk,jb) = ( p_vort_e(1:nlen,jk,jb)             &
                                   +pt_patch%edges%f_e(1:nlen,jb) )    &
                                   /p_delp_e(1:nlen,jk,jb)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  SELECT CASE (i_cori_method)
  CASE (1,2)

    ! Tendency of vn caused by vorticity flux

    ieidx => pt_int%heli_vn_idx
    ieblk => pt_int%heli_vn_blk
    SELECT CASE (i_cori_method)
    CASE (1)
      nincr = 14
    CASE (2)
      nincr = 10
    END SELECT
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,ji,je)
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
#ifdef __LOOP_EXCHANGE
      DO je = 1, nlen
        DO ji = 1, nincr
          DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO ji = 1, nincr
          DO je = 1, nlen
#endif
            p_ddt_vn(je,jk,jb) =  p_ddt_vn(je,jk,jb)             &
            & +pt_int%heli_coeff(ji,je,jb)                       &
            & *p_mflux_e(ieidx(ji,je,jb),jk,ieblk(ji,je,jb))     &
            & *( z_potvort_r(ieidx(ji,je,jb),jk,ieblk(ji,je,jb)) &
            &   +z_potvort_r(je             ,jk,jb             ))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  CASE (3) ! i_cori_method

    ieidx => pt_patch%edges%quad_idx
    ieblk => pt_patch%edges%quad_blk
    icidx => pt_patch%edges%cell_idx
    icblk => pt_patch%edges%cell_blk

    ! first, reconstruct mass flux vectors at centers of rhombi and hexagons
    CALL edges2cells_scalar(p_mflux_e,pt_patch,pt_int%hex_east  ,z_u_c)
    CALL edges2cells_scalar(p_mflux_e,pt_patch,pt_int%hex_north ,z_v_c)
    CALL edges2edges_scalar(p_mflux_e,pt_patch,pt_int%quad_east ,z_u_e)
    CALL edges2edges_scalar(p_mflux_e,pt_patch,pt_int%quad_north,z_v_e)

    ! second, average absolute potential vorticity from rhombi to centers
    IF (l_corner_vort) THEN
      CALL edges2verts_scalar(z_potvort_r,pt_patch,pt_int%e_1o3_v       ,z_tmp_v)
      CALL sync_patch_array(SYNC_V,pt_patch,z_tmp_v)
      CALL verts2cells_scalar(z_tmp_v    ,pt_patch,pt_int%verts_aw_cells,z_potvort_c)
    ELSE
      CALL edges2cells_scalar(z_potvort_r,pt_patch,pt_int%e_aw_c,z_potvort_c)
    ENDIF

    ! third, multiply the absolute vorticities with the velocities,
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk)
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO jk = 1,nlev
        z_u_e(1:nlen,jk,jb) = z_u_e(1:nlen,jk,jb)*z_potvort_r(1:nlen,jk,jb)
        z_v_e(1:nlen,jk,jb) = z_v_e(1:nlen,jk,jb)*z_potvort_r(1:nlen,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP DO PRIVATE(jb,nlen,jk)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      DO jk = 1,nlev
        z_u_c(1:nlen,jk,jb) = z_u_c(1:nlen,jk,jb)*z_potvort_c(1:nlen,jk,jb)
        z_v_c(1:nlen,jk,jb) = z_v_c(1:nlen,jk,jb)*z_potvort_c(1:nlen,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
    ! fourth, compute vorticity flux term
!$OMP DO PRIVATE(jb,nlen,jk,je)
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
#ifdef __LOOP_EXCHANGE
      DO je = 1, nlen
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = 1, nlen
#endif
          p_ddt_vn(je,jk,jb) =  p_ddt_vn(je,jk,jb)             &
          & + pt_int%heli_coeff( 1,je,jb)*z_v_c(icidx(je,jb,1),jk,icblk(je,jb,1))&
          & + pt_int%heli_coeff( 2,je,jb)*z_u_c(icidx(je,jb,1),jk,icblk(je,jb,1))&
          & + pt_int%heli_coeff( 3,je,jb)*z_v_c(icidx(je,jb,2),jk,icblk(je,jb,2))&
          & + pt_int%heli_coeff( 4,je,jb)*z_u_c(icidx(je,jb,2),jk,icblk(je,jb,2))&
          & + pt_int%heli_coeff( 5,je,jb)*z_v_e(je            ,jk,jb            )&
          & + pt_int%heli_coeff( 6,je,jb)*z_u_e(je            ,jk,jb            )&
          & + pt_int%heli_coeff( 7,je,jb)*z_v_e(ieidx(je,jb,1),jk,ieblk(je,jb,1))&
          & + pt_int%heli_coeff( 8,je,jb)*z_u_e(ieidx(je,jb,1),jk,ieblk(je,jb,1))&
          & + pt_int%heli_coeff( 9,je,jb)*z_v_e(ieidx(je,jb,2),jk,ieblk(je,jb,2))&
          & + pt_int%heli_coeff(10,je,jb)*z_u_e(ieidx(je,jb,2),jk,ieblk(je,jb,2))&
          & + pt_int%heli_coeff(11,je,jb)*z_v_e(ieidx(je,jb,3),jk,ieblk(je,jb,3))&
          & + pt_int%heli_coeff(12,je,jb)*z_u_e(ieidx(je,jb,3),jk,ieblk(je,jb,3))&
          & + pt_int%heli_coeff(13,je,jb)*z_v_e(ieidx(je,jb,4),jk,ieblk(je,jb,4))&
          & + pt_int%heli_coeff(14,je,jb)*z_u_e(ieidx(je,jb,4),jk,ieblk(je,jb,4))
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  CASE(4) ! the same as i_cori_method=3 but instead edge PVs use corner PVs

    ieidx => pt_patch%edges%quad_idx
    ieblk => pt_patch%edges%quad_blk
    icidx => pt_patch%edges%cell_idx
    icblk => pt_patch%edges%cell_blk
    ividx => pt_patch%edges%vertex_idx
    ivblk => pt_patch%edges%vertex_blk

    ! first, reconstruct mass flux vectors at centers of rhombi and hexagons
    CALL edges2cells_scalar(p_mflux_e,pt_patch,pt_int%hex_east  ,z_u_c)
    CALL edges2cells_scalar(p_mflux_e,pt_patch,pt_int%hex_north ,z_v_c)
    CALL edges2edges_scalar(p_mflux_e,pt_patch,pt_int%quad_east ,z_u_e)
    CALL edges2edges_scalar(p_mflux_e,pt_patch,pt_int%quad_north,z_v_e)

    ! second, average absolute potential vorticity from rhombi to vertices and centers
    CALL edges2verts_scalar(z_potvort_r,pt_patch,pt_int%e_1o3_v,z_potvort_v)
    CALL sync_patch_array(SYNC_V,pt_patch,z_potvort_v)
    IF (l_corner_vort) THEN
      CALL verts2cells_scalar(z_potvort_v,pt_patch,pt_int%verts_aw_cells,z_potvort_c)
    ELSE
      CALL edges2cells_scalar(z_potvort_r,pt_patch,pt_int%e_aw_c,z_potvort_c)
    ENDIF

!$OMP PARALLEL
    ! third, compute vorticity flux term
!$OMP DO PRIVATE(jb,nlen,jk,je)
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
#ifdef __LOOP_EXCHANGE
      DO je = 1, nlen
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = 1, nlen
#endif
          p_ddt_vn(je,jk,jb) =  p_ddt_vn(je,jk,jb)             &
          & + (pt_int%heli_coeff( 1,je,jb)*z_v_c(icidx(je,jb,1),jk,icblk(je,jb,1)) &
          &   +pt_int%heli_coeff( 2,je,jb)*z_u_c(icidx(je,jb,1),jk,icblk(je,jb,1)))&
          &                         *z_potvort_c(icidx(je,jb,1),jk,icblk(je,jb,1)) &
          & + (pt_int%heli_coeff( 3,je,jb)*z_v_c(icidx(je,jb,2),jk,icblk(je,jb,2)) &
          &   +pt_int%heli_coeff( 4,je,jb)*z_u_c(icidx(je,jb,2),jk,icblk(je,jb,2)))&
          &                         *z_potvort_c(icidx(je,jb,2),jk,icblk(je,jb,2)) &
          & + (pt_int%heli_coeff( 5,je,jb)*z_v_e(je            ,jk,jb            ) &
          &   +pt_int%heli_coeff( 6,je,jb)*z_u_e(je            ,jk,jb            ))&
          &                 *0.5_wp*(z_potvort_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) &
          &                         +z_potvort_v(ividx(je,jb,2),jk,ivblk(je,jb,2)))&
          & + (pt_int%heli_coeff( 7,je,jb)*z_v_e(ieidx(je,jb,1),jk,ieblk(je,jb,1)) &
          &   +pt_int%heli_coeff( 8,je,jb)*z_u_e(ieidx(je,jb,1),jk,ieblk(je,jb,1)) &
          &   +pt_int%heli_coeff( 9,je,jb)*z_v_e(ieidx(je,jb,2),jk,ieblk(je,jb,2)) &
          &   +pt_int%heli_coeff(10,je,jb)*z_u_e(ieidx(je,jb,2),jk,ieblk(je,jb,2)))&
          &                         *z_potvort_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) &
          & + (pt_int%heli_coeff(11,je,jb)*z_v_e(ieidx(je,jb,3),jk,ieblk(je,jb,3)) &
          &   +pt_int%heli_coeff(12,je,jb)*z_u_e(ieidx(je,jb,3),jk,ieblk(je,jb,3)) &
          &   +pt_int%heli_coeff(13,je,jb)*z_v_e(ieidx(je,jb,4),jk,ieblk(je,jb,4)) &
          &   +pt_int%heli_coeff(14,je,jb)*z_u_e(ieidx(je,jb,4),jk,ieblk(je,jb,4)))&
          &                         *z_potvort_v(ividx(je,jb,2),jk,ivblk(je,jb,2))
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SELECT ! i_cori_method

END SELECT ! i_cell_type

!=============================================================================
! Calculate the gradient of kinetic energy, and accumulate velocity tendency
!=============================================================================

  CALL grad_fd_norm( p_kin, pt_patch, z_tmp_e, opt_rlstart=4 )

  jbs = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1)
      p_ddt_vn(is:ie,:,jb) = p_ddt_vn(is:ie,:,jb) - z_tmp_e(is:ie,:,jb)
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

   END SUBROUTINE vn_adv_horizontal

END MODULE mo_ha_dynamics_adv
