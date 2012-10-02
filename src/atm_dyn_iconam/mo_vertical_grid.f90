!>
!! mo_vertical_grid provides all subroutines for handling of the vertical
!! grid in the nonhydrostatic model branch.
!! Routines are provided for
!!  1) Creating the hybrid height based coordinate system.
!!  2) (De)Allocations of related metric fields. Note, that they are defined
!!     in mo_model_domain
!!
!! @author Almut Gassmann, MPI-M
!!
!! @par Revision History
!! Initial release by Almut Gassmann (2009-04-14)
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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


MODULE mo_vertical_grid

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_model_domain,          ONLY: t_patch
  USE mo_ext_data_types,        ONLY: t_external_data
  USE mo_grid_config,           ONLY: n_dom
  USE mo_nonhydrostatic_config, ONLY: rayleigh_type, rayleigh_coeff, damp_height, &
    &                                 igradp_method, ivctype, vwind_offctr,       &
    &                                 exner_expol, l_zdiffu_t, thslp_zdiffu,      &
    &                                 thhgtd_zdiffu
  USE mo_diffusion_config,      ONLY: diffusion_config
  USE mo_parallel_config,       ONLY: nproma, p_test_run
  USE mo_run_config,            ONLY: msg_level
  USE mo_sleve_config,          ONLY: lread_smt
  USE mo_vertical_coord_table,  ONLY: vct_a
  USE mo_impl_constants,        ONLY: MAX_CHAR_LENGTH, max_dom, RAYLEIGH_CLASSIC, &
    &                                 RAYLEIGH_KLEMP, min_rlcell_int, min_rlcell
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c, grf_bdywidth_e, grf_fbk_start_c
  USE mo_physical_constants,    ONLY: grav, p0ref, rd, rd_o_cpd, cpd, p0sl_bg
  USE mo_math_gradients,        ONLY: grad_fd_norm, grad_fd_tang
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_intp,                  ONLY: cells2edges_scalar, cells2verts_scalar, verts2edges_scalar
  USE mo_math_constants,        ONLY: pi_2
  USE mo_loopindices,           ONLY: get_indices_e, get_indices_c
  USE mo_nonhydro_types,        ONLY: t_nh_state
  USE mo_nh_init_utils,         ONLY: nflat, nflatlev, compute_smooth_topo, init_vert_coord
  USE mo_sync,                  ONLY: SYNC_E, SYNC_V, sync_patch_array, global_sum_array, &
                                      sync_patch_array_mult, global_min, global_max

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  ! Constants used for the computation of the background reference atmosphere of the nh-model
  !
  REAL(wp), PARAMETER :: h_scal_bg = 10000._wp    ! [m]      scale height
  REAL(wp), PARAMETER :: t0sl_bg   = 288.15_wp    ! [K]      sea level temperature
  REAL(wp), PARAMETER :: del_t_bg  = 75._wp       ! [K]      difference between sea level
  !                                                          temperature and asymptotic
  !                                                          stratospheric temperature
  REAL(wp), PARAMETER :: grav_o_cpd      = grav/cpd

  INTEGER:: nrdmax(max_dom), nflat_gradp(max_dom)

  PUBLIC :: set_nh_metrics, nrdmax, nflat_gradp

  CONTAINS


  !----------------------------------------------------------------------------
  !>
  !! Initializes values for the fields in type nh_metrics of the NH state.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-04-14)
  !! New formulation of the tangential slope by A. Gassmann (2009-11-17)
  !! Modification by Almut Gassmann (2009-11-17)
  !! - Adding Rayleigh damping coeff. at upper boundary
  !!
  SUBROUTINE set_nh_metrics(p_patch, p_nh, p_int, ext_data)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_vertical_grid:set_nh_metrics'

    TYPE(t_patch), TARGET, INTENT(INOUT) :: p_patch(n_dom)  !< patch
    TYPE(t_nh_state),      INTENT(INOUT) :: p_nh(n_dom)
    TYPE(t_int_state),     INTENT(   IN) :: p_int(n_dom)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(n_dom)

    INTEGER :: jg, jk, jk1, jk_start, jb, jc, je, jv, jn, jgc, nlen, &
               nblks_c, npromz_c, nblks_e, npromz_e, nblks_v, npromz_v, ic
    INTEGER :: nlev, nlevp1              !< number of full levels

    ! Note: the present way of setting up coordinate surfaces will not work for vertical refinement
    INTEGER :: i_startidx, i_endidx, i_startblk, i_endblk, i_nchdom, icount_total
    INTEGER :: ica(max_dom)

    REAL(wp) :: z_diff, z_sin_diff, z_sin_diff_full, z_tanh_diff,    &
      &         z1, z2, z3, z_help(nproma),                          &
      &         z_temp(nproma), z_aux1(nproma), z_aux2(nproma),      &
      &         z0, coef1, coef2, coef3, dn1, dn2, dn3, dn4, dn5, dn6
    REAL(wp) :: z_maxslope, z_maxhdiff, z_offctr
    REAL(wp), ALLOCATABLE :: z_ifv(:,:,:), z_mfv(:,:,:)
    REAL(wp), ALLOCATABLE :: z_me(:,:,:),z_maxslp(:,:,:),z_maxhgtd(:,:,:),z_shift(:,:,:)
    REAL(wp) :: extrapol_dist
    INTEGER,  ALLOCATABLE :: flat_idx(:,:), imask(:,:,:),icount(:)
    INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
    LOGICAL :: l_half_lev_centr, l_found(nproma), lfound_all

    !------------------------------------------------------------------------

    DO jg = 1,n_dom

      SELECT CASE (p_patch(jg)%cell_type)
      CASE (6)
        l_half_lev_centr = .TRUE.
        ! The HALF LEVEL where the model layer are flat, moves one layer upward.
        ! there could also be a zero there
        nflat = nflat-1
      CASE DEFAULT
        l_half_lev_centr = .FALSE.
      END SELECT

      nblks_c   = p_patch(jg)%nblks_int_c
      npromz_c  = p_patch(jg)%npromz_int_c

      i_nchdom   = MAX(1,p_patch(jg)%n_childdom)

      ! number of vertical levels
      nlev   = p_patch(jg)%nlev
      nlevp1 = p_patch(jg)%nlevp1

      ! total shift of model top with respect to global domain
      IF (jg > 1) THEN
        nflatlev(jg)     = nflatlev(1) - p_patch(jg)%nshift_total
      ENDIF
      IF (p_patch(jg)%cell_type == 6) nflatlev(jg) = nflat
      IF (jg > 1 .AND. p_patch(jg)%nshift_total > 0 .AND. nflatlev(jg) <= 1) THEN
        CALL finish (TRIM(routine), &
                     'flat_height too close to the top of the innermost nested domain')
      ENDIF

      ! Compute smooth topography when SLEVE coordinate is used
      IF ( ivctype == 2 .AND. .NOT. lread_smt ) THEN
        CALL compute_smooth_topo(p_patch(jg), p_int(jg), ext_data(jg)%atm%topography_c, &
                      ext_data(jg)%atm%topography_smt_c, ext_data(jg)%atm%topography_v, &
                      ext_data(jg)%atm%topography_smt_v                                 )
      ENDIF

      ! Initialize vertical coordinate for cell points
      CALL init_vert_coord(ext_data(jg)%atm%topography_c, ext_data(jg)%atm%topography_smt_c, &
                           p_nh(jg)%metrics%z_ifc, p_nh(jg)%metrics%z_mc, nlev,              &
                           nblks_c, npromz_c, p_patch(jg)%nshift_total, nflatlev(jg),        &
                           l_half_lev_centr                                                  )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1,nblks_c

        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
           ! To avoid uninitialized field elements going in the output routine
           p_nh(jg)%metrics%geopot(nlen+1:nproma,:,jb) = 0._wp
        ENDIF

        DO jk = 1, nlev
          ! geopotential on full levels
          p_nh(jg)%metrics%geopot(1:nlen,jk,jb) = grav*p_nh(jg)%metrics%z_mc(1:nlen,jk,jb)
        ENDDO

        DO jk = 1, nlev
          ! functional determinant of the metrics (is positive), full levels
          p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk,jb) = &
          & p_nh(jg)%metrics%z_ifc(1:nlen,jk  ,jb)- &
          & p_nh(jg)%metrics%z_ifc(1:nlen,jk+1,jb)
          ! difference of geopotential between the levels
          p_nh(jg)%metrics%dgeopot_mc (1:nlen,jk,jb) = grav * &
          &         ( p_nh(jg)%metrics%z_ifc(1:nlen,jk  ,jb) &
          &         - p_nh(jg)%metrics%z_ifc(1:nlen,jk+1,jb) )
          ! inverse layer thickness (for runtime optimization)
          p_nh(jg)%metrics%inv_ddqz_z_full(1:nlen,jk,jb) = &
            1._wp/p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk,jb)
        ENDDO
        ! functional determinant of the metrics (is positive), half levels
        DO jk = 2, nlev
          p_nh(jg)%metrics%ddqz_z_half(1:nlen,jk,jb) = &
          & (p_nh(jg)%metrics%z_mc(1:nlen,jk-1,jb)     &
          & -p_nh(jg)%metrics%z_mc(1:nlen,jk  ,jb))
          IF (p_patch(jg)%cell_type == 6) THEN
            p_nh(jg)%metrics%diff_1_o_dz(1:nlen,jk,jb) = &
            &  1.0_wp/p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk-1,jb) & 
            & -1.0_wp/p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk  ,jb)
            p_nh(jg)%metrics%mult_1_o_dz(1:nlen,jk,jb) = &
            &  1.0_wp/(p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk-1,jb) & 
            &         *p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk  ,jb))
          ENDIF
        ENDDO
        IF (p_patch(jg)%cell_type == 6) THEN
          p_nh(jg)%metrics%diff_1_o_dz(1:nlen,1,jb) = 0.0_wp
          p_nh(jg)%metrics%mult_1_o_dz(1:nlen,1,jb) = 0.0_wp
          p_nh(jg)%metrics%diff_1_o_dz(1:nlen,nlevp1,jb) = 0.0_wp
          p_nh(jg)%metrics%mult_1_o_dz(1:nlen,nlevp1,jb) = 0.0_wp
        ENDIF
        p_nh(jg)%metrics%ddqz_z_half(1:nlen,1,jb) =   &
        & 2.0_wp*(p_nh(jg)%metrics%z_ifc(1:nlen,1,jb) &
        &       - p_nh(jg)%metrics%z_mc (1:nlen,1,jb))
        p_nh(jg)%metrics%ddqz_z_half(1:nlen,nlevp1,jb) =    &
        & 2.0_wp*(p_nh(jg)%metrics%z_mc (1:nlen,nlev  ,jb)  &
        &       - p_nh(jg)%metrics%z_ifc(1:nlen,nlevp1,jb))
        IF (p_patch(jg)%cell_type==3) THEN
          ! layer distance between jk+1 and jk-1
          p_nh(jg)%metrics%inv_ddqz_z_half2(:,1,jb) = 0._wp
          DO jk = 2, nlev
            p_nh(jg)%metrics%inv_ddqz_z_half2(1:nlen,jk,jb) = 1._wp / &
              ( p_nh(jg)%metrics%z_ifc(1:nlen,jk-1,jb) -              &
                p_nh(jg)%metrics%z_ifc(1:nlen,jk+1,jb) )
          ENDDO
        ENDIF

        !-------------------------------------------------------------
        ! geopot above ground  - because physics needs positive values
        !-------------------------------------------------------------

          p_nh(jg)%metrics%geopot_agl_ifc(1:nlen,nlevp1,jb) = 0._wp     

          DO jk = nlev,1,-1
            ! geopotential (interfaces)
            p_nh(jg)%metrics%geopot_agl_ifc(1:nlen,jk,jb)=             &
          &            p_nh(jg)%metrics%geopot_agl_ifc(1:nlen,jk+1,jb) &
          &         +  grav *                                      &
          &          ( p_nh(jg)%metrics%z_ifc(1:nlen,jk  ,jb)      &
          &          - p_nh(jg)%metrics%z_ifc(1:nlen,jk+1,jb))
          ENDDO

          ! geopotential above ground (at full levels)
          DO jk = 2,nlevp1
            p_nh(jg)%metrics%geopot_agl(1:nlen,jk-1,jb) =  0.5_wp *           &
           &              ( p_nh(jg)%metrics%geopot_agl_ifc(1:nlen,jk-1  ,jb )&
           &              + p_nh(jg)%metrics%geopot_agl_ifc(1:nlen,jk    ,jb))
          ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (msg_level >= 15) THEN
   !  WRITE(message_text,'(a,3E15.5)') 'vertical height and GEOPOT  = ',&
   !       &  MINVAL(p_nh(jg)%metrics%geopot_agl(:,:,:)),&
   !       &  MINVAL( p_nh(jg)%metrics%geopot_agl_ifc(:,:,:)),&
   !       &  MINVAL(p_nh(jg)%metrics%dgeopot_mc(:,:,:))
   !     CALL message('', TRIM(message_text))

      DO jk = 1, nlev
       WRITE(message_text,'(a,i4,3E15.5)') 'GEOPOT full/half,dgeopot  = ',jk,&
          &  p_nh(jg)%metrics%geopot_agl(1,jk,2), p_nh(jg)%metrics%geopot_agl_ifc(1,jk,2),&
          &  p_nh(jg)%metrics%dgeopot_mc(1,jk,2)
        CALL message(TRIM(routine), TRIM(message_text))
      ENDDO
    ENDIF

      ! For the tangential slope we need temporarily also the height
      ! at the vertices
      ALLOCATE(z_ifv(nproma,nlevp1,p_patch(jg)%nblks_v))
      ALLOCATE(z_mfv(nproma,nlev  ,p_patch(jg)%nblks_v))
      nblks_v   = p_patch(jg)%nblks_int_v
      npromz_v  = p_patch(jg)%npromz_int_v

      ! Initialize vertical coordinate for vertex points
      CALL init_vert_coord(ext_data(jg)%atm%topography_v, ext_data(jg)%atm%topography_smt_v, &
                           z_ifv, z_mfv, nlev, nblks_v, npromz_v, p_patch(jg)%nshift_total,  &
                           nflatlev(jg), l_half_lev_centr                                    )

      ! Start index for slope computations
      i_startblk = p_patch(jg)%edges%start_blk(2,1)

      p_nh(jg)%metrics%ddxt_z_half = 0._wp
      p_nh(jg)%metrics%ddxt_z_full = 0._wp
      p_nh(jg)%metrics%ddxn_z_half = 0._wp
      p_nh(jg)%metrics%ddxn_z_full = 0._wp

      ! slope of the terrain (tangential direction)
      CALL grad_fd_tang ( z_ifv, &
           &              p_patch(jg), &
           &              p_nh(jg)%metrics%ddxt_z_half, &
           &              1, nlevp1 )
      DEALLOCATE(z_ifv)
      DEALLOCATE(z_mfv)

      ! slope of the terrain (normal direction)
      CALL grad_fd_norm ( p_nh(jg)%metrics%z_ifc, &
           &              p_patch(jg), &
           &              p_nh(jg)%metrics%ddxn_z_half, &
           &              1, nlevp1 )
      IF (l_half_lev_centr) THEN
        CALL grad_fd_norm ( p_nh(jg)%metrics%z_mc, &
             &              p_patch(jg), &
             &              p_nh(jg)%metrics%ddxn_z_full, &
             &              1, nlev )
      ENDIF

      IF (p_patch(jg)%cell_type == 6) THEN
        IF (p_test_run) p_nh(jg)%metrics%z_mc_e = 0._wp

        CALL cells2edges_scalar(p_nh(jg)%metrics%z_mc, &
                                p_patch(jg),p_int(jg)%c_lin_e, &
                                p_nh(jg)%metrics%z_mc_e )

        CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%z_mc_e)
      ENDIF

      CALL sync_patch_array_mult(SYNC_E,p_patch(jg),3,p_nh(jg)%metrics%ddxt_z_half,       &
                                 p_nh(jg)%metrics%ddxn_z_half,p_nh(jg)%metrics%ddxn_z_full)

      ! vertically averaged metrics
      nblks_e   = p_patch(jg)%nblks_int_e
      npromz_e  = p_patch(jg)%npromz_int_e
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,nblks_e

        CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, 2)

        IF (.NOT.l_half_lev_centr) THEN
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
              p_nh(jg)%metrics%ddxn_z_full(je,jk,jb) = 0.5_wp * &
              & (p_nh(jg)%metrics%ddxn_z_half(je,jk  ,jb) + &
              &  p_nh(jg)%metrics%ddxn_z_half(je,jk+1,jb))

              p_nh(jg)%metrics%ddxt_z_full(je,jk,jb) = 0.5_wp * &
              & (p_nh(jg)%metrics%ddxt_z_half(je,jk  ,jb) + &
              &  p_nh(jg)%metrics%ddxt_z_half(je,jk+1,jb))
            ENDDO
          ENDDO
        ENDIF

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (p_test_run) p_nh(jg)%metrics%ddqz_z_full_e = 0._wp

      ! functional determinant at full level edges
      CALL cells2edges_scalar(p_nh(jg)%metrics%ddqz_z_full, &
           &                  p_patch(jg), p_int(jg)%c_lin_e, &
           &                  p_nh(jg)%metrics%ddqz_z_full_e )

      CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%ddqz_z_full_e)

      IF(p_patch(jg)%cell_type==6)THEN

        IF (p_test_run) THEN
          p_nh(jg)%metrics%ddqz_z_half_e = 0._wp
          p_nh(jg)%metrics%ddqz_z_full_v = 0._wp
          p_nh(jg)%metrics%ddqz_z_full_r = 0._wp
        ENDIF
        ! functional determinant at half level the edges
        CALL cells2edges_scalar(p_nh(jg)%metrics%ddqz_z_half, &
                                p_patch(jg),p_int(jg)%c_lin_e, &
                                p_nh(jg)%metrics%ddqz_z_half_e, 1, nlevp1 )

        CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%ddqz_z_half_e)

        ! functional determinant at full level dual centers
        CALL cells2verts_scalar(p_nh(jg)%metrics%ddqz_z_full, &
                                p_patch(jg),p_int(jg)%cells_aw_verts, &
                                p_nh(jg)%metrics%ddqz_z_full_v, 1, nlev)

        CALL sync_patch_array(SYNC_V,p_patch(jg), p_nh(jg)%metrics%ddqz_z_full_v)

        CALL verts2edges_scalar(p_nh(jg)%metrics%ddqz_z_full_v, &
                                p_patch(jg), p_int(jg)%tria_aw_rhom, &
                                p_nh(jg)%metrics%ddqz_z_full_r)

        CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%ddqz_z_full_r)

        iidx => p_patch(jg)%verts%edge_idx
        iblk => p_patch(jg)%verts%edge_blk
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, jv) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1,nblks_v
          IF (jb /= nblks_v) THEN
             nlen = nproma
          ELSE
             nlen = npromz_v
          ENDIF
          DO jk = 1, nlev
            DO jv = 1, nlen
              p_nh(jg)%metrics%ddnorth_z(jv,jk,jb) = &
              &   p_nh(jg)%metrics%ddxn_z_full(iidx(jv,jb,1),jk,iblk(jv,jb,1)) &
              & * p_int(jg)%tria_north(1,jv,jb) &
              & + p_nh(jg)%metrics%ddxn_z_full(iidx(jv,jb,2),jk,iblk(jv,jb,2)) &
              & * p_int(jg)%tria_north(2,jv,jb) &
              & + p_nh(jg)%metrics%ddxn_z_full(iidx(jv,jb,3),jk,iblk(jv,jb,3)) &
              & * p_int(jg)%tria_north(3,jv,jb)
              p_nh(jg)%metrics%ddeast_z(jv,jk,jb) = &
              &   p_nh(jg)%metrics%ddxn_z_full(iidx(jv,jb,1),jk,iblk(jv,jb,1)) &
              & * p_int(jg)%tria_east(1,jv,jb) &
              & + p_nh(jg)%metrics%ddxn_z_full(iidx(jv,jb,2),jk,iblk(jv,jb,2)) &
              & * p_int(jg)%tria_east(2,jv,jb) &
              & + p_nh(jg)%metrics%ddxn_z_full(iidx(jv,jb,3),jk,iblk(jv,jb,3)) &
              & * p_int(jg)%tria_east(3,jv,jb)
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO NOWAIT

        ! The reference Exner pressure field is needed for physics initialization
!$OMP DO PRIVATE(jb, nlen, jk, z_aux1) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1,nblks_c
          IF (jb /= nblks_c) THEN
             nlen = nproma
          ELSE
             nlen = npromz_c
          ENDIF

          DO jk = 1, nlev
            ! Reference pressure, full level mass points
            z_aux1(1:nlen) = p0sl_bg*EXP(-grav/rd*h_scal_bg/(t0sl_bg-del_t_bg) &
              &  *LOG((EXP(p_nh(jg)%metrics%z_mc(1:nlen,jk,jb)/h_scal_bg)      &
              &  *(t0sl_bg-del_t_bg) +del_t_bg)/t0sl_bg))

            ! Reference Exner pressure, full level mass points
            p_nh(jg)%metrics%exner_ref_mc(1:nlen,jk,jb) = (z_aux1(1:nlen)/p0ref)**rd_o_cpd
          ENDDO
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      ENDIF

      ! offcentering in vertical mass flux 
      p_nh(jg)%metrics%vwind_expl_wgt(:,:) = 0.5_wp - vwind_offctr
      p_nh(jg)%metrics%vwind_impl_wgt(:,:) = 0.5_wp + vwind_offctr

      ! Rayleigh damping properties

      p_nh(jg)%metrics%rayleigh_w(:)   = 0.0_wp
      p_nh(jg)%metrics%enhfac_diffu(:) = 1.0_wp

      ! Determine end index of damping layer
      nrdmax(jg) = 1
      DO jk = 2, nlevp1
        jk1 = jk + p_patch(jg)%nshift_total
        IF (vct_a(jk1) >= damp_height(jg))  nrdmax(jg) = jk
      ENDDO


      ! Rayleigh damping coefficient for vn and/or w
      DO jk = 1, nrdmax(jg)
        jk1 = jk + p_patch(jg)%nshift_total

        ! z - z_h  with z at half levels
        z_sin_diff  = MAX(0.0_wp,vct_a(jk1)-damp_height(jg))

        ! z_top - z    
        z_tanh_diff = vct_a(1) - vct_a(jk1)

        ! z - z_h  with z at full levels
        z_sin_diff_full = MAX(0.0_wp,0.5_wp*(vct_a(jk1)+vct_a(jk1+1))-damp_height(jg))

        IF ( rayleigh_type == RAYLEIGH_CLASSIC ) THEN
          ! classical Rayleigh damping as known from idealized limited area test cases.
          ! Requires knowledge of a velocity reference state!
          !
          p_nh(jg)%metrics%rayleigh_vn(jk)= rayleigh_coeff(jg)*(SIN(pi_2*z_sin_diff_full/ &
              MAX(1.e-3_wp,vct_a(p_patch(jg)%nshift_total+1)-damp_height(jg))))**2

          p_nh(jg)%metrics%rayleigh_w(jk)= rayleigh_coeff(jg)*(SIN(pi_2*z_sin_diff/ &
              MAX(1.e-3_wp,vct_a(p_patch(jg)%nshift_total+1)-damp_height(jg))))**2

        ELSE IF ( rayleigh_type == RAYLEIGH_KLEMP ) THEN
          ! Rayleigh damping based on Klemp et al. (2008), MRW 136, pp 3987-4004
          ! No reference state needed, thus applicable to real cases!
          !
          p_nh(jg)%metrics%rayleigh_w(jk)= rayleigh_coeff(jg)*&
          (1._wp-TANH(3.8_wp*z_tanh_diff/MAX(1.e-6_wp,vct_a(1)-damp_height(jg))))

          p_nh(jg)%metrics%rayleigh_vn(jk)= 0._wp
        ENDIF
      ENDDO


      ! Enhancement coefficient for nabla4 background diffusion near model top
      DO jk = 1, nrdmax(jg)
        jk1 = jk + p_patch(jg)%nshift_total
!        z_diff = MAX(0.0_wp,0.5_wp*(vct_a(jk1)+vct_a(jk1+1))-damp_height(jg))
        z_diff = 0.5_wp*(vct_a(1)+vct_a(2))-0.5_wp*(vct_a(jk1)+vct_a(jk1+1))
        p_nh(jg)%metrics%enhfac_diffu(jk) = 1._wp + &
 !         (hdiff_efdt_ratio/hdiff_min_efdt_ratio-1._wp)*(SIN(pi_2*z_diff/ &
 !         MAX(1.e-3_wp,0.5_wp*(vct_a(1)+vct_a(2))-damp_height(jg))))**2
          (diffusion_config(jg)%hdiff_efdt_ratio/  &
           diffusion_config(jg)%hdiff_min_efdt_ratio-1._wp)* &
          (1._wp-TANH(3.8_wp*z_diff                          &
          /MAX(1.e-6_wp,0.5_wp*(vct_a(1)+vct_a(2))-damp_height(jg))))
      ENDDO

      IF (p_patch(jg)%cell_type == 3 .AND. msg_level >= 10) THEN
        WRITE(message_text,'(a,i4,a,i4)') 'Domain', jg, &
          '; end index of Rayleigh damping layer for w: ', nrdmax(jg)
        CALL message(TRIM(routine),message_text)
        WRITE(message_text,'(a)') &
          'Damping coefficient for w; diffusion enhancement coefficient:'
        CALL message('mo_vertical_grid',message_text)
        DO jk = 1, nrdmax(jg)
          jk1 = jk + p_patch(jg)%nshift_total
          WRITE(message_text,'(a,i5,a,f8.1,2e13.5)') 'level',jk,', half-level height',vct_a(jk1),&
            p_nh(jg)%metrics%rayleigh_w(jk),p_nh(jg)%metrics%enhfac_diffu(jk)
          CALL message('mo_vertical_grid',message_text)
        ENDDO
      ENDIF


      ! Compute variable Exner extrapolation factors and offcentering coefficients for the
      ! vertical wind solver to optimize numerical stability over steep orography
      iidx => p_patch(jg)%cells%edge_idx
      iblk => p_patch(jg)%cells%edge_blk

      i_startblk = p_patch(jg)%cells%start_blk(2,1)

      IF (p_patch(jg)%cell_type == 3) THEN
        ALLOCATE (z_maxslp(nproma,nlev,p_patch(jg)%nblks_c), &
                  z_maxhgtd(nproma,nlev,p_patch(jg)%nblks_c) )

!$OMP PARALLEL
!$OMP WORKSHARE
        ! Initialization to ensure that values are properly set at lateral boundaries
        p_nh(jg)%metrics%exner_exfac(:,:,:) = exner_expol

        z_maxslp(:,:,:)  = 0._wp
        z_maxhgtd(:,:,:) = 0._wp
!$OMP END WORKSHARE

!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jc, z_maxslope, z_offctr) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,nblks_c

          CALL get_indices_c(p_patch(jg), jb, i_startblk, nblks_c, &
                             i_startidx, i_endidx, 2)

          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx

              z_maxslp(jc,jk,jb) =                                                     &
                MAX(ABS(p_nh(jg)%metrics%ddxn_z_full(iidx(jc,jb,1),jk,iblk(jc,jb,1))), &
                    ABS(p_nh(jg)%metrics%ddxn_z_full(iidx(jc,jb,2),jk,iblk(jc,jb,2))), &
                    ABS(p_nh(jg)%metrics%ddxn_z_full(iidx(jc,jb,3),jk,iblk(jc,jb,3))) )

              z_maxhgtd(jc,jk,jb) =                                                  &
                MAX(ABS(p_nh(jg)%metrics%ddxn_z_full(iidx(jc,jb,1),jk,iblk(jc,jb,1)) &
                 * p_patch(jg)%edges%dual_edge_length(iidx(jc,jb,1),iblk(jc,jb,1))), &
                    ABS(p_nh(jg)%metrics%ddxn_z_full(iidx(jc,jb,2),jk,iblk(jc,jb,2)) &
                 * p_patch(jg)%edges%dual_edge_length(iidx(jc,jb,2),iblk(jc,jb,2))), &
                    ABS(p_nh(jg)%metrics%ddxn_z_full(iidx(jc,jb,3),jk,iblk(jc,jb,3)) &
                 * p_patch(jg)%edges%dual_edge_length(iidx(jc,jb,3),iblk(jc,jb,3)) ) )

              ! Exner extrapolation reaches zero for a slope of 1/4 or a height difference of 500 m
              ! between adjacent grid points (empirically determined values)
              p_nh(jg)%metrics%exner_exfac(jc,jk,jb) = exner_expol * &
                MIN(1._wp-(4._wp*z_maxslp(jc,jk,jb))**2._wp,         &
                    1._wp-(0.002_wp*z_maxhgtd(jc,jk,jb))**2._wp)
              p_nh(jg)%metrics%exner_exfac(jc,jk,jb) =  &
                MAX(0._wp,p_nh(jg)%metrics%exner_exfac(jc,jk,jb))
              ! For extremely steep slopes, going a bit behind time level nnow turned out
              ! to further improve stability
              IF (z_maxslp(jc,jk,jb) > 1.5_wp) &
                p_nh(jg)%metrics%exner_exfac(jc,jk,jb) = &
                  MAX(-1._wp/6._wp,1._wp/9._wp*(1.5_wp-z_maxslp(jc,jk,jb)))
            ENDDO
          ENDDO

          jk = nlevp1
          DO jc = i_startidx, i_endidx

            z_maxslope = MAX(ABS(p_nh(jg)%metrics%ddxn_z_half(iidx(jc,jb,1),jk,iblk(jc,jb,1))),&
                             ABS(p_nh(jg)%metrics%ddxn_z_half(iidx(jc,jb,2),jk,iblk(jc,jb,2))),&
                             ABS(p_nh(jg)%metrics%ddxn_z_half(iidx(jc,jb,3),jk,iblk(jc,jb,3))),&
                             ABS(p_nh(jg)%metrics%ddxt_z_half(iidx(jc,jb,1),jk,iblk(jc,jb,1))),&
                             ABS(p_nh(jg)%metrics%ddxt_z_half(iidx(jc,jb,2),jk,iblk(jc,jb,2))),&
                             ABS(p_nh(jg)%metrics%ddxt_z_half(iidx(jc,jb,3),jk,iblk(jc,jb,3))) )

            ! Empirically determined values to ensure stability over steep slopes
            z_offctr =   MAX(vwind_offctr, 0.425_wp*z_maxslope**0.75_wp)
            z_offctr =   MIN(MAX(vwind_offctr,0.5_wp),z_offctr)

            p_nh(jg)%metrics%vwind_expl_wgt(jc,jb) = 0.5_wp - z_offctr
            p_nh(jg)%metrics%vwind_impl_wgt(jc,jb) = 0.5_wp + z_offctr

          ENDDO
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      ENDIF

      ! The remaining computations are needed for the triangular grid only
      ! once the initialization in mo_nh_testcases is properly rewritten for hexagons
      IF (p_patch(jg)%cell_type == 6) CYCLE

      ! Index lists for boundary nudging (including halo cells so that no 
      ! sync is needed afterwards; halo edges are excluded, however, because
      ! a sync follows afterwards anyway
      ic = 0
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        DO jc = 1, nlen
          IF (p_int(jg)%nudgecoeff_c(jc,jb) > 1.e-10_wp) THEN
            ic = ic+1
          ENDIF
        ENDDO
      ENDDO
      p_nh(jg)%metrics%nudge_c_dim = ic

      ALLOCATE(p_nh(jg)%metrics%nudge_c_idx(ic),p_nh(jg)%metrics%nudge_c_blk(ic))
      ic = 0
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        DO jc = 1, nlen
          IF (p_int(jg)%nudgecoeff_c(jc,jb) > 1.e-10_wp) THEN
            ic = ic+1
            p_nh(jg)%metrics%nudge_c_idx(ic) = jc
            p_nh(jg)%metrics%nudge_c_blk(ic) = jb
          ENDIF
        ENDDO
      ENDDO

      ic = 0
      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
        DO je = 1, nlen
          IF (p_patch(jg)%edges%owner_mask(je,jb) .AND.     &
              p_int(jg)%nudgecoeff_e(je,jb) > 1.e-10_wp) THEN
            ic = ic+1
          ENDIF
        ENDDO
      ENDDO
      p_nh(jg)%metrics%nudge_e_dim = ic

      ALLOCATE(p_nh(jg)%metrics%nudge_e_idx(ic),p_nh(jg)%metrics%nudge_e_blk(ic))

      ic = 0
      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
        DO je = 1, nlen
          IF (p_patch(jg)%edges%owner_mask(je,jb) .AND.     &
              p_int(jg)%nudgecoeff_e(je,jb) > 1.e-10_wp) THEN
            ic = ic+1
            p_nh(jg)%metrics%nudge_e_idx(ic) = je
            p_nh(jg)%metrics%nudge_e_blk(ic) = jb
          ENDIF
        ENDDO
      ENDDO

      ! Index lists needed to minimize the number of halo communications in solve_nh and feedback
      ! (Remark: setting i_startblk for i_nchdom is necessary to get an empty index list
      !  for non-MPI-parallelized setups with multiple nests per nest level)

      p_nh(jg)%metrics%mask_prog_halo_c(:,:) = .FALSE.

      i_startblk = p_patch(jg)%cells%start_blk(min_rlcell_int-1,i_nchdom)
      i_endblk   = p_patch(jg)%cells%end_blk(min_rlcell,i_nchdom)

      ic = 0
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, min_rlcell_int-1, min_rlcell, i_nchdom)

        DO jc = i_startidx, i_endidx
          IF (p_patch(jg)%cells%refin_ctrl(jc,jb)>=1 .AND. &
              p_patch(jg)%cells%refin_ctrl(jc,jb)<=4) THEN
            ic = ic+1
          ENDIF
        ENDDO
      ENDDO
      p_nh(jg)%metrics%bdy_halo_c_dim = ic

      ! Index list for halo points belonging to the lateral boundary interpolation zone
      ALLOCATE(p_nh(jg)%metrics%bdy_halo_c_idx(ic),p_nh(jg)%metrics%bdy_halo_c_blk(ic))

      ic = 0
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, min_rlcell_int-1, min_rlcell, i_nchdom)

        DO jc = i_startidx, i_endidx
          IF (p_patch(jg)%cells%refin_ctrl(jc,jb)>=1 .AND. &
              p_patch(jg)%cells%refin_ctrl(jc,jb)<=4) THEN
            ic = ic+1
            p_nh(jg)%metrics%bdy_halo_c_idx(ic) = jc
            p_nh(jg)%metrics%bdy_halo_c_blk(ic) = jb
          ELSE
            p_nh(jg)%metrics%mask_prog_halo_c(jc,jb) = .TRUE.
          ENDIF
        ENDDO
      ENDDO

      ! Index list for halo points belonging to the nest overlap zone

      i_startblk = p_patch(jg)%cells%start_blk(min_rlcell_int-1,i_nchdom)
      i_endblk   = p_patch(jg)%cells%end_blk(min_rlcell,i_nchdom)
      ica(:)     = 0

      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, min_rlcell_int-1, min_rlcell, i_nchdom)

        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)

          DO jc = i_startidx, i_endidx
            IF (p_patch(jg)%cells%refin_ctrl(jc,jb)<=grf_fbk_start_c .AND. &
                p_patch(jg)%cells%child_id(jc,jb)==jgc) THEN
              ica(jn) = ica(jn)+1
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ic = MAXVAL(ica)
      jn = MAX(1,p_patch(jg)%n_childdom)
      ALLOCATE(p_nh(jg)%metrics%ovlp_halo_c_dim(jn),   &
               p_nh(jg)%metrics%ovlp_halo_c_idx(ic,jn),&
               p_nh(jg)%metrics%ovlp_halo_c_blk(ic,jn))

      p_nh(jg)%metrics%ovlp_halo_c_dim(1:jn) = ica(1:jn)

      ica(:) = 0
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, min_rlcell_int-1, min_rlcell, i_nchdom)

        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)

          DO jc = i_startidx, i_endidx
            IF (p_patch(jg)%cells%refin_ctrl(jc,jb)<=grf_fbk_start_c .AND. &
                p_patch(jg)%cells%child_id(jc,jb)==jgc) THEN
              ica(jn) = ica(jn)+1
              p_nh(jg)%metrics%ovlp_halo_c_idx(ica(jn),jn) = jc
              p_nh(jg)%metrics%ovlp_halo_c_blk(ica(jn),jn) = jb

            ENDIF
          ENDDO
        ENDDO
      ENDDO

      IF (l_zdiffu_t) THEN
        CALL prepare_zdiffu(p_patch(jg), p_nh(jg), p_int(jg), z_maxslp, z_maxhgtd)
      ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, jc, z1, z2, z3) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        DO jk = 2, nlev
          p_nh(jg)%metrics%wgtfac_c(1:nlen,jk,jb) =  &
           (p_nh(jg)%metrics%z_ifc(1:nlen,jk-1,jb) - &
            p_nh(jg)%metrics%z_ifc(1:nlen,jk,jb)) /  &
           (p_nh(jg)%metrics%z_ifc(1:nlen,jk-1,jb) - &
            p_nh(jg)%metrics%z_ifc(1:nlen,jk+1,jb))
        ENDDO
          p_nh(jg)%metrics%wgtfac_c(1:nlen,1,jb) = &
           (p_nh(jg)%metrics%z_ifc(1:nlen,2,jb) -  &
            p_nh(jg)%metrics%z_ifc(1:nlen,1,jb)) / &
           (p_nh(jg)%metrics%z_ifc(1:nlen,3,jb) -  &
            p_nh(jg)%metrics%z_ifc(1:nlen,1,jb))
          p_nh(jg)%metrics%wgtfac_c(1:nlen,nlevp1,jb) = &
           (p_nh(jg)%metrics%z_ifc(1:nlen,nlev,jb) -    &
            p_nh(jg)%metrics%z_ifc(1:nlen,nlevp1,jb)) / &
           (p_nh(jg)%metrics%z_ifc(1:nlen,nlev-1,jb) -  &
            p_nh(jg)%metrics%z_ifc(1:nlen,nlevp1,jb))
        DO jc = 1, nlen
          z1 = 0.5_wp * ( p_nh(jg)%metrics%z_ifc(jc,nlev,jb) - &
                          p_nh(jg)%metrics%z_ifc(jc,nlevp1,jb) )
          z2 = 0.5_wp * ( p_nh(jg)%metrics%z_ifc(jc,nlev,jb) +     &
                          p_nh(jg)%metrics%z_ifc(jc,nlev-1,jb) ) - &
                          p_nh(jg)%metrics%z_ifc(jc,nlevp1,jb)
          z3 = 0.5_wp * ( p_nh(jg)%metrics%z_ifc(jc,nlev-1,jb) +   &
                          p_nh(jg)%metrics%z_ifc(jc,nlev-2,jb) ) - &
                          p_nh(jg)%metrics%z_ifc(jc,nlevp1,jb)
          p_nh(jg)%metrics%wgtfacq_c(jc,3,jb) = z1*z2/(z2-z3)/(z1-z3)
          p_nh(jg)%metrics%wgtfacq_c(jc,2,jb) = &
            (z1-p_nh(jg)%metrics%wgtfacq_c(jc,3,jb)*(z1-z3))/(z1-z2)
          p_nh(jg)%metrics%wgtfacq_c(jc,1,jb) = 1._wp - &
           (p_nh(jg)%metrics%wgtfacq_c(jc,2,jb) +       &
            p_nh(jg)%metrics%wgtfacq_c(jc,3,jb))
        ENDDO
        DO jc = 1, nlen
          z1 = 0.5_wp * ( p_nh(jg)%metrics%z_ifc(jc,2,jb) -   &
                          p_nh(jg)%metrics%z_ifc(jc,1,jb) )
          z2 = 0.5_wp * ( p_nh(jg)%metrics%z_ifc(jc,2,jb) +   &
                          p_nh(jg)%metrics%z_ifc(jc,3,jb) ) - &
                          p_nh(jg)%metrics%z_ifc(jc,1,jb)
          z3 = 0.5_wp * ( p_nh(jg)%metrics%z_ifc(jc,3,jb) +   &
                          p_nh(jg)%metrics%z_ifc(jc,4,jb) ) - &
                          p_nh(jg)%metrics%z_ifc(jc,1,jb)
          p_nh(jg)%metrics%wgtfacq1_c(jc,3,jb) = z1*z2/(z2-z3)/(z1-z3)
          p_nh(jg)%metrics%wgtfacq1_c(jc,2,jb) = &
            (z1-p_nh(jg)%metrics%wgtfacq1_c(jc,3,jb)*(z1-z3))/(z1-z2)
          p_nh(jg)%metrics%wgtfacq1_c(jc,1,jb) = 1._wp - &
           (p_nh(jg)%metrics%wgtfacq1_c(jc,2,jb) +       &
            p_nh(jg)%metrics%wgtfacq1_c(jc,3,jb))
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (p_patch(jg)%cell_type == 3 .AND. msg_level >= 10) THEN
        z_offctr = MAXVAL(p_nh(jg)%metrics%vwind_impl_wgt,MASK=p_patch(jg)%cells%owner_mask(:,:))
        z_offctr = global_max(z_offctr) - 0.5_wp
        DO jk = 1, nlev
          WHERE(.NOT.p_patch(jg)%cells%owner_mask(:,:))
            z_maxslp(:,jk,:)  = -HUGE(0._wp)
            z_maxhgtd(:,jk,:) = -HUGE(0._wp)
          END WHERE
        ENDDO
        z_maxslope = MAXVAL(z_maxslp)
        z_maxslope = global_max(z_maxslope)
        z_maxhdiff = MAXVAL(z_maxhgtd)
        z_maxhdiff = global_max(z_maxhdiff)
        WRITE(message_text,'(a,f8.4)') 'Maximum vertical wind offcentering: ', z_offctr
        CALL message(TRIM(routine),message_text)
        WRITE(message_text,'(a,f8.4)') 'Maximum slope: ', z_maxslope
        CALL message(TRIM(routine),message_text)
        WRITE(message_text,'(a,f8.1)') 'Maximum height difference between adjacent points: ', &
          z_maxhdiff
        CALL message(TRIM(routine),message_text)
      ENDIF

      DEALLOCATE (z_maxslp,z_maxhgtd)

      IF (p_test_run) THEN
        p_nh(jg)%metrics%wgtfac_e = 0._wp
        p_nh(jg)%metrics%wgtfacq_e = 0._wp
        p_nh(jg)%metrics%wgtfacq1_e = 0._wp
      ENDIF

      ! Interpolate weighting coefficients to edges
      CALL cells2edges_scalar(p_nh(jg)%metrics%wgtfac_c, p_patch(jg), &
             p_int(jg)%c_lin_e, p_nh(jg)%metrics%wgtfac_e, 1, nlevp1)

      CALL cells2edges_scalar(p_nh(jg)%metrics%wgtfacq_c, p_patch(jg), &
             p_int(jg)%c_lin_e, p_nh(jg)%metrics%wgtfacq_e, 1, 3)

      CALL cells2edges_scalar(p_nh(jg)%metrics%wgtfacq1_c, p_patch(jg), &
             p_int(jg)%c_lin_e, p_nh(jg)%metrics%wgtfacq1_e, 1, 3)

      CALL sync_patch_array_mult(SYNC_E,p_patch(jg),3,p_nh(jg)%metrics%wgtfac_e,       &
                                 p_nh(jg)%metrics%wgtfacq_e,p_nh(jg)%metrics%wgtfacq1_e)

      IF (p_patch(jg)%cell_type == 3) THEN

      ! Reference atmosphere fields for triangular code
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, z_help, z_temp, z_aux1, z_aux2) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1,nblks_c
          IF (jb /= nblks_c) THEN
             nlen = nproma
          ELSE
             nlen = npromz_c
          ENDIF

          ! Reference surface temperature
          p_nh(jg)%metrics%tsfc_ref(1:nlen,jb) = (t0sl_bg-del_t_bg)+del_t_bg*   &
            & EXP(-p_nh(jg)%metrics%z_ifc(1:nlen,nlevp1,jb)/h_scal_bg)

          DO jk = 1, nlev
            ! Reference pressure, full level mass points
            z_aux1(1:nlen) = p0sl_bg*EXP(-grav/rd*h_scal_bg/(t0sl_bg-del_t_bg) &
              &  *LOG((EXP(p_nh(jg)%metrics%z_mc(1:nlen,jk,jb)/h_scal_bg)      &
              &  *(t0sl_bg-del_t_bg) +del_t_bg)/t0sl_bg))

            ! Reference Exner pressure, full level mass points
            p_nh(jg)%metrics%exner_ref_mc(1:nlen,jk,jb) = (z_aux1(1:nlen)/p0ref)**rd_o_cpd

            ! Reference temperature, full level mass points
            z_temp(1:nlen) = (t0sl_bg-del_t_bg)+del_t_bg*        &
              & EXP(-p_nh(jg)%metrics%z_mc(1:nlen,jk,jb)/h_scal_bg)

            ! Reference density, full level mass points
            p_nh(jg)%metrics%rho_ref_mc(1:nlen,jk,jb) = z_aux1(1:nlen)/(rd*z_temp(1:nlen))

            ! Reference Potential temperature, full level mass points
            p_nh(jg)%metrics%theta_ref_mc(1:nlen,jk,jb) = z_temp(1:nlen) &
              & /p_nh(jg)%metrics%exner_ref_mc(1:nlen,jk,jb)
          ENDDO

          IF (igradp_method <= 3) THEN
            DO jk = 1, nlev
              ! First vertical derivative of reference Exner pressure, full level mass points,
              ! divided by theta_ref
              ! Note: for computational efficiency, this field is in addition divided by
              ! the vertical layer thickness
              p_nh(jg)%metrics%d2dexdz2_fac1_mc(1:nlen,jk,jb)   =             &
                & -grav/(cpd*p_nh(jg)%metrics%theta_ref_mc(1:nlen,jk,jb)**2)* &
                & p_nh(jg)%metrics%inv_ddqz_z_full(1:nlen,jk,jb)

              ! Vertical derivative of d_exner_dz/theta_ref, full level mass points
              p_nh(jg)%metrics%d2dexdz2_fac2_mc(1:nlen,jk,jb)   =                            &
                &  2._wp*grav/(cpd*p_nh(jg)%metrics%theta_ref_mc(1:nlen,jk,jb)**3)*(grav/cpd &
                & -del_t_bg/h_scal_bg*EXP(-p_nh(jg)%metrics%z_mc(1:nlen,jk,jb)/h_scal_bg))   &
                & /p_nh(jg)%metrics%exner_ref_mc(1:nlen,jk,jb)
            ENDDO
          ENDIF

          DO jk = 1, nlevp1
            ! Reference pressure, half level mass points
            z_aux1(1:nlen) = p0sl_bg*EXP(-grav/rd*h_scal_bg/(t0sl_bg-del_t_bg) &
              &  *LOG((EXP(p_nh(jg)%metrics%z_ifc(1:nlen,jk,jb)/h_scal_bg)      &
              &  *(t0sl_bg-del_t_bg) +del_t_bg)/t0sl_bg))

            ! Reference Exner pressure, half level mass points
            z_help(1:nlen) = (z_aux1(1:nlen)/p0ref)**rd_o_cpd

            ! Reference temperature, half level mass points
            z_temp(1:nlen) = (t0sl_bg-del_t_bg)+del_t_bg*          &
              & EXP(-p_nh(jg)%metrics%z_ifc(1:nlen,jk,jb)/h_scal_bg)

            ! Reference density, half level mass points
            z_aux2(1:nlen) = z_aux1(1:nlen)/(rd*z_temp(1:nlen))

            ! Reference Potential temperature, half level mass points
            p_nh(jg)%metrics%theta_ref_ic(1:nlen,jk,jb) = z_temp(1:nlen)/z_help(1:nlen)

            ! First vertical derivative of reference Exner pressure, half level mass points
            p_nh(jg)%metrics%d_exner_dz_ref_ic(1:nlen,jk,jb)   =       &
              & -grav/cpd/p_nh(jg)%metrics%theta_ref_ic(1:nlen,jk,jb)

          ENDDO
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        ALLOCATE(z_me(nproma,nlev,p_patch(jg)%nblks_e))
        IF (p_test_run) THEN
          z_me = 0._wp
        ENDIF

        ! Compute geometric height at edge points
        CALL cells2edges_scalar(p_nh(jg)%metrics%z_mc, p_patch(jg), &
               p_int(jg)%c_lin_e, z_me)

        CALL sync_patch_array(SYNC_E,p_patch(jg),z_me)

        i_startblk = p_patch(jg)%edges%start_blk(2,1)

        ! Reference fields on edges
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, je, z_aux1, z_temp) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,nblks_e

          CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                             i_startidx, i_endidx, 2)

          DO jk = 1, nlev
            DO je = i_startidx, i_endidx

              ! Reference pressure, full level edge points
              z_aux1(je) = p0sl_bg*EXP(-grav/rd*h_scal_bg/(t0sl_bg-del_t_bg) &
                &  *LOG((EXP(z_me(je,jk,jb)/h_scal_bg)      &
                &  *(t0sl_bg-del_t_bg) +del_t_bg)/t0sl_bg))

              ! Reference temperature, full level edge points
              z_temp(je) = (t0sl_bg-del_t_bg)+del_t_bg*        &
                & EXP(-z_me(je,jk,jb)/h_scal_bg)

              ! Reference density, full level edge points
              p_nh(jg)%metrics%rho_ref_me(je,jk,jb) = z_aux1(je)/(rd*z_temp(je))

              ! Reference Potential temperature, full level edge points
              p_nh(jg)%metrics%theta_ref_me(je,jk,jb) = z_temp(je)/((z_aux1(je)/p0ref)**rd_o_cpd)
            ENDDO
          ENDDO

        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ENDIF

      ! Compute information needed for Taylor-expansion-based pressure gradient calculation
      IF (igradp_method == 1) THEN
        nflat_gradp(jg) = nlev
        DEALLOCATE(z_me)
        CYCLE
      ENDIF

      IF (nflatlev(jg) <= 2 .AND. igradp_method >= 4) THEN
        CALL finish (TRIM(routine),'flat_height must be at least 2 levels below top&
          & for igradp_method>3')
      ENDIF

      ALLOCATE(flat_idx(nproma,p_patch(jg)%nblks_e))
      flat_idx = nlev

      iidx => p_patch(jg)%edges%cell_idx
      iblk => p_patch(jg)%edges%cell_blk


!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,nblks_e

        CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, 2)

        IF (igradp_method <= 3) THEN
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
              p_nh(jg)%metrics%vertidx_gradp(1:2,je,jk,jb) = jk
              p_nh(jg)%metrics%zdiff_gradp(1,je,jk,jb)     =  &
                z_me(je,jk,jb) - p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk,iblk(je,jb,1))
              p_nh(jg)%metrics%zdiff_gradp(2,je,jk,jb)     =  &
                z_me(je,jk,jb) - p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk,iblk(je,jb,2))
            ENDDO
          ENDDO
        ENDIF

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

            ! Compute the highest vertical index (counted from top to bottom) for which
            ! the edge point lies inside the cell box of the adjacent grid points
            IF (z_me(je,jk,jb) <= p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),jk,iblk(je,jb,1))   .AND. &
                z_me(je,jk,jb) >= p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),jk+1,iblk(je,jb,1)) .AND. &
                z_me(je,jk,jb) <= p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),jk,iblk(je,jb,2))   .AND. &
                z_me(je,jk,jb) >= p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),jk+1,iblk(je,jb,2))) THEN
              flat_idx(je,jb) = jk
            ENDIF

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! Compute global minimum of flat_idx

      ! Attention: Edges in the halo of flat_idx are set to "random" values.
      ! So we have to consider only inner edges (or to do a boundary exchange)
      ! when searching the minimum.
      ! Please note also that a patch may be completely empty on some processors,
      ! in this case MINVAL returns HUGE() so there is no special care needed.

      nflat_gradp(jg) = MINVAL(flat_idx(:,:), MASK=p_patch(jg)%edges%owner_mask(:,:))
      nflat_gradp(jg) = NINT(global_min(REAL(nflat_gradp(jg),wp)))

      ! Compute level indices of neighbor cells within which the local edge is located
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jk1, jk_start, je, l_found, lfound_all, &
!$OMP z0, z1, z2, z3, coef1, coef2, coef3, dn1, dn2, dn3, dn4, dn5, &
!$OMP dn6) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,nblks_e

        CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, 2)

        IF (igradp_method <= 3) THEN 

          DO je = i_startidx, i_endidx

            jk_start = flat_idx(je,jb)
            DO jk = flat_idx(je,jb)+1,nlev
              DO jk1 = jk_start, nlev
                IF ( (jk1 == nlev) .OR. z_me(je,jk,jb) <=                         &
                    p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),jk1,iblk(je,jb,1)) .AND. &
                    z_me(je,jk,jb) >=                                             &
                    p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),jk1+1,iblk(je,jb,1))) THEN

                  p_nh(jg)%metrics%vertidx_gradp(1,je,jk,jb) = jk1
                  p_nh(jg)%metrics%zdiff_gradp(1,je,jk,jb)   = z_me(je,jk,jb) -   &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1,iblk(je,jb,1))
                  jk_start = jk1
                  EXIT
                ENDIF
              ENDDO
            ENDDO
            jk_start = flat_idx(je,jb)
            DO jk = flat_idx(je,jb)+1,nlev
              DO jk1 = jk_start, nlev
                IF ( (jk1 == nlev) .OR. z_me(je,jk,jb) <=                         &
                    p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),jk1,iblk(je,jb,2)) .AND. &
                    z_me(je,jk,jb) >=                                             &
                    p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),jk1+1,iblk(je,jb,2))) THEN

                  p_nh(jg)%metrics%vertidx_gradp(2,je,jk,jb) = jk1
                  p_nh(jg)%metrics%zdiff_gradp(2,je,jk,jb)   = z_me(je,jk,jb) -   &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1,iblk(je,jb,2))
                  jk_start = jk1
                  EXIT
                ENDIF
              ENDDO
            ENDDO

          ENDDO

        ELSE IF (igradp_method >= 4) THEN ! Coefficients for polynomial interpolation

          jk_start = nflatlev(jg) - 1
          DO jk = nflatlev(jg),nlev
            l_found(:) = .FALSE.
            lfound_all = .FALSE.
            DO jk1 = jk_start, nlev - 2
              DO je = i_startidx, i_endidx
                IF (z_me(je,jk,jb) <=                                            &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1,iblk(je,jb,1)) .AND. &
                    z_me(je,jk,jb) >                                             &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1+1,iblk(je,jb,1))) THEN

                  ! cubic interpolation
                  p_nh(jg)%metrics%vertidx_gradp(1,je,jk,jb) = jk1

                  z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1-1,iblk(je,jb,1))
                  z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1  ,iblk(je,jb,1))
                  z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1+1,iblk(je,jb,1))
                  z3 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1+2,iblk(je,jb,1))

                  coef1 =  z_me(je,jk,jb)-z0
                  coef2 = (z_me(je,jk,jb)-z0)*(z_me(je,jk,jb)-z1)
                  coef3 = (z_me(je,jk,jb)-z0)*(z_me(je,jk,jb)-z1)*(z_me(je,jk,jb)-z2)
                  
                  dn1 = 1._wp/(z0-z1)
                  dn2 = 1._wp/(z1-z2)
                  dn3 = 1._wp/(z2-z3)
                  dn4 = 1._wp/(z0-z2)
                  dn5 = 1._wp/(z0-z3)
                  dn6 = 1._wp/(z1-z3)
                  
                  p_nh(jg)%metrics%coeff_gradp(1,je,jk,jb) =            &
                    1._wp + coef1*dn1 + coef2*dn1*dn4 + coef3*dn1*dn4*dn5
                  p_nh(jg)%metrics%coeff_gradp(2,je,jk,jb) =                               &
                    -(coef1*dn1 + coef2*dn4*(dn1+dn2) + coef3*dn5*(dn1*dn4+dn2*dn4+dn2*dn6))
                  p_nh(jg)%metrics%coeff_gradp(3,je,jk,jb) =          &
                    coef2*dn2*dn4 + coef3*dn5*(dn2*dn4+dn2*dn6+dn3*dn6)
                  p_nh(jg)%metrics%coeff_gradp(4,je,jk,jb) = -coef3*dn3*dn5*dn6

                  l_found(je) = .TRUE.
                ELSE IF (z_me(je,jk,jb) <=                                        &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev-1,iblk(je,jb,1))) THEN

                  ! quadratic interpolation/extrapolation
                  p_nh(jg)%metrics%vertidx_gradp(1,je,jk,jb) = nlev - 1

                  z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev-2,iblk(je,jb,1))
                  z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev-1,iblk(je,jb,1))
                  z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev  ,iblk(je,jb,1))

                  coef1 =  z_me(je,jk,jb)-z0
                  coef2 = (z_me(je,jk,jb)-z0)*(z_me(je,jk,jb)-z1)
                  
                  dn1 = 1._wp/(z0-z1)
                  dn2 = 1._wp/(z1-z2)
                  dn4 = 1._wp/(z0-z2)
                  
                  p_nh(jg)%metrics%coeff_gradp(1,je,jk,jb) = 1._wp + coef1*dn1 + coef2*dn1*dn4
                  p_nh(jg)%metrics%coeff_gradp(2,je,jk,jb) = -(coef1*dn1 + coef2*dn4*(dn1+dn2))
                  p_nh(jg)%metrics%coeff_gradp(3,je,jk,jb) = coef2*dn2*dn4
                  p_nh(jg)%metrics%coeff_gradp(4,je,jk,jb) = 0._wp

                  l_found(je) = .TRUE.
                ENDIF
              ENDDO
              IF (ALL(l_found(i_startidx:i_endidx))) THEN
                lfound_all = .TRUE.
                EXIT
              ENDIF
            ENDDO
            IF (lfound_all) THEN
              jk_start = MIN(nlev-2, &
                MINVAL(p_nh(jg)%metrics%vertidx_gradp(1,i_startidx:i_endidx,jk,jb)))
            ENDIF
          ENDDO

          jk_start = nflatlev(jg) - 1
          DO jk = nflatlev(jg),nlev
            l_found(:) = .FALSE.
            lfound_all = .FALSE.
            DO jk1 = jk_start, nlev - 2
              DO je = i_startidx, i_endidx
                IF (z_me(je,jk,jb) <=                                            &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1,iblk(je,jb,2)) .AND. &
                    z_me(je,jk,jb) >                                             &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1+1,iblk(je,jb,2))) THEN

                  ! cubic interpolation
                  p_nh(jg)%metrics%vertidx_gradp(2,je,jk,jb) = jk1

                  z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1-1,iblk(je,jb,2))
                  z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1  ,iblk(je,jb,2))
                  z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1+1,iblk(je,jb,2))
                  z3 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1+2,iblk(je,jb,2))

                  coef1 =  z_me(je,jk,jb)-z0
                  coef2 = (z_me(je,jk,jb)-z0)*(z_me(je,jk,jb)-z1)
                  coef3 = (z_me(je,jk,jb)-z0)*(z_me(je,jk,jb)-z1)*(z_me(je,jk,jb)-z2)
                  
                  dn1 = 1._wp/(z0-z1)
                  dn2 = 1._wp/(z1-z2)
                  dn3 = 1._wp/(z2-z3)
                  dn4 = 1._wp/(z0-z2)
                  dn5 = 1._wp/(z0-z3)
                  dn6 = 1._wp/(z1-z3)
                  
                  p_nh(jg)%metrics%coeff_gradp(5,je,jk,jb) =            &
                    1._wp + coef1*dn1 + coef2*dn1*dn4 + coef3*dn1*dn4*dn5
                  p_nh(jg)%metrics%coeff_gradp(6,je,jk,jb) =                               &
                    -(coef1*dn1 + coef2*dn4*(dn1+dn2) + coef3*dn5*(dn1*dn4+dn2*dn4+dn2*dn6))
                  p_nh(jg)%metrics%coeff_gradp(7,je,jk,jb) =          &
                    coef2*dn2*dn4 + coef3*dn5*(dn2*dn4+dn2*dn6+dn3*dn6)
                  p_nh(jg)%metrics%coeff_gradp(8,je,jk,jb) = -coef3*dn3*dn5*dn6

                  l_found(je) = .TRUE.
                ELSE IF (z_me(je,jk,jb) <=                                        &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev-1,iblk(je,jb,2))) THEN

                  ! quadratic interpolation/extrapolation
                  p_nh(jg)%metrics%vertidx_gradp(2,je,jk,jb) = nlev - 1

                  z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev-2,iblk(je,jb,2))
                  z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev-1,iblk(je,jb,2))
                  z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev  ,iblk(je,jb,2))

                  coef1 =  z_me(je,jk,jb)-z0
                  coef2 = (z_me(je,jk,jb)-z0)*(z_me(je,jk,jb)-z1)
                  
                  dn1 = 1._wp/(z0-z1)
                  dn2 = 1._wp/(z1-z2)
                  dn4 = 1._wp/(z0-z2)
                  
                  p_nh(jg)%metrics%coeff_gradp(5,je,jk,jb) = 1._wp + coef1*dn1 + coef2*dn1*dn4
                  p_nh(jg)%metrics%coeff_gradp(6,je,jk,jb) = -(coef1*dn1 + coef2*dn4*(dn1+dn2))
                  p_nh(jg)%metrics%coeff_gradp(7,je,jk,jb) = coef2*dn2*dn4
                  p_nh(jg)%metrics%coeff_gradp(8,je,jk,jb) = 0._wp

                  l_found(je) = .TRUE.
                ENDIF
              ENDDO
              IF (ALL(l_found(i_startidx:i_endidx))) THEN
                lfound_all = .TRUE.
                EXIT
              ENDIF
            ENDDO
            IF (lfound_all) THEN
              jk_start = MIN(nlev-2, &
                MINVAL(p_nh(jg)%metrics%vertidx_gradp(2,i_startidx:i_endidx,jk,jb)))
            ENDIF
          ENDDO

        ENDIF

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (igradp_method == 3 .OR. igradp_method == 5) THEN

        i_startblk = p_patch(jg)%edges%start_blk(grf_bdywidth_e+1,1)

        ALLOCATE(imask(nproma,nlev,p_patch(jg)%nblks_e),icount(p_patch(jg)%nblks_e), &
                 z_shift(nproma,nlev,p_patch(jg)%nblks_e) )

        extrapol_dist = 5._wp ! maximum allowed extrapolation distance; may become a namelist variable later on

        ! Recompute indices and height differences if truly horizontal pressure gradient 
        ! computation would intersect the ground
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jk1, jk_start, je, z_aux1, &
!$OMP z_aux2) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,nblks_e

          CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                             i_startidx, i_endidx, grf_bdywidth_e+1)

          imask(:,:,jb) = 0
          icount(jb)    = 0

          DO je = i_startidx, i_endidx

            z_aux1(je) = MAX(p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),nlevp1,iblk(je,jb,1)), &
                             p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),nlevp1,iblk(je,jb,2)))

            z_aux2(je) = z_aux1(je) - extrapol_dist ! allow for some limited downward extrapolation

            DO jk = flat_idx(je,jb)+1,nlev
              IF ( z_me(je,jk,jb) < z_aux2(je)) THEN

                ! Save information needed for index list setup
                IF (p_patch(jg)%edges%owner_mask(je,jb)) THEN
                  imask(je,jk,jb) = 1
                  icount(jb)      = icount(jb) + 1
                  z_shift(je,jk,jb) = z_me(je,jk,jb) - z_aux2(je)
                ENDIF
              ENDIF
            ENDDO

            IF (igradp_method == 3) THEN
              jk_start = flat_idx(je,jb)
              DO jk = flat_idx(je,jb)+1,nlev
                IF ( z_me(je,jk,jb) < z_aux2(je)) THEN
                  DO jk1 = jk_start, nlev
                    IF ( jk1 == nlev .OR. z_aux2(je) <=                             &
                      p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),jk1,iblk(je,jb,1)) .AND. &
                      z_aux2(je) >=                                                 &
                      p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),jk1+1,iblk(je,jb,1))) THEN

                      p_nh(jg)%metrics%vertidx_gradp(1,je,jk,jb) = jk1
                      p_nh(jg)%metrics%zdiff_gradp(1,je,jk,jb)   = z_aux2(je) -     &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1,iblk(je,jb,1))
                      jk_start = jk1
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
              jk_start = flat_idx(je,jb)
              DO jk = flat_idx(je,jb)+1,nlev
                IF ( z_me(je,jk,jb) < z_aux2(je)) THEN
                  DO jk1 = jk_start, nlev
                    IF ( jk1 == nlev .OR. z_aux2(je) <=                             &
                      p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),jk1,iblk(je,jb,2)) .AND. &
                      z_aux2(je) >=                                                 &
                      p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),jk1+1,iblk(je,jb,2))) THEN

                      p_nh(jg)%metrics%vertidx_gradp(2,je,jk,jb) = jk1
                      p_nh(jg)%metrics%zdiff_gradp(2,je,jk,jb)   = z_aux2(je) -     &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1,iblk(je,jb,2))
                      jk_start = jk1
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ELSE IF (igradp_method == 5) THEN
              jk_start = flat_idx(je,jb)
              DO jk = flat_idx(je,jb)+1,nlev
                IF ( z_me(je,jk,jb) < z_aux2(je)) THEN
                  DO jk1 = jk_start, nlev - 2
                    IF (z_aux2(je) <=                                                &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1,iblk(je,jb,1)) .AND. &
                        z_aux2(je) >                                                 &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1+1,iblk(je,jb,1))) THEN

                      ! cubic interpolation
                      p_nh(jg)%metrics%vertidx_gradp(1,je,jk,jb) = jk1

                      z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1-1,iblk(je,jb,1))
                      z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1  ,iblk(je,jb,1))
                      z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1+1,iblk(je,jb,1))
                      z3 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1+2,iblk(je,jb,1))

                      coef1 =  z_aux2(je)-z0
                      coef2 = (z_aux2(je)-z0)*(z_aux2(je)-z1)
                      coef3 = (z_aux2(je)-z0)*(z_aux2(je)-z1)*(z_aux2(je)-z2)
                  
                      dn1 = 1._wp/(z0-z1)
                      dn2 = 1._wp/(z1-z2)
                      dn3 = 1._wp/(z2-z3)
                      dn4 = 1._wp/(z0-z2)
                      dn5 = 1._wp/(z0-z3)
                      dn6 = 1._wp/(z1-z3)
                  
                      p_nh(jg)%metrics%coeff_gradp(1,je,jk,jb) =            &
                        1._wp + coef1*dn1 + coef2*dn1*dn4 + coef3*dn1*dn4*dn5
                      p_nh(jg)%metrics%coeff_gradp(2,je,jk,jb) =                               &
                        -(coef1*dn1 + coef2*dn4*(dn1+dn2) + coef3*dn5*(dn1*dn4+dn2*dn4+dn2*dn6))
                      p_nh(jg)%metrics%coeff_gradp(3,je,jk,jb) =          &
                        coef2*dn2*dn4 + coef3*dn5*(dn2*dn4+dn2*dn6+dn3*dn6)
                      p_nh(jg)%metrics%coeff_gradp(4,je,jk,jb) = -coef3*dn3*dn5*dn6

                      jk_start = jk1
                      EXIT

                    ELSE IF (z_aux2(je) <=                                        &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev-1,iblk(je,jb,1))) THEN

                      ! quadratic interpolation/extrapolation
                      p_nh(jg)%metrics%vertidx_gradp(1,je,jk,jb) = nlev - 1

                      z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev-2,iblk(je,jb,1))
                      z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev-1,iblk(je,jb,1))
                      z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev  ,iblk(je,jb,1))

                      coef1 =  z_aux2(je)-z0
                      coef2 = (z_aux2(je)-z0)*(z_aux2(je)-z1)
                  
                      dn1 = 1._wp/(z0-z1)
                      dn2 = 1._wp/(z1-z2)
                      dn4 = 1._wp/(z0-z2)
                  
                      p_nh(jg)%metrics%coeff_gradp(1,je,jk,jb) = &
                        1._wp + coef1*dn1 + coef2*dn1*dn4
                      p_nh(jg)%metrics%coeff_gradp(2,je,jk,jb) = &
                        -(coef1*dn1 + coef2*dn4*(dn1+dn2))
                      p_nh(jg)%metrics%coeff_gradp(3,je,jk,jb) = coef2*dn2*dn4
                      p_nh(jg)%metrics%coeff_gradp(4,je,jk,jb) = 0._wp

                      jk_start = nlev - 2
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO

              jk_start = flat_idx(je,jb)
              DO jk = flat_idx(je,jb)+1,nlev
                IF ( z_me(je,jk,jb) < z_aux2(je)) THEN
                  DO jk1 = jk_start, nlev - 2
                    IF (z_aux2(je) <=                                                &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1,iblk(je,jb,2)) .AND. &
                        z_aux2(je) >                                                 &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1+1,iblk(je,jb,2))) THEN

                      ! cubic interpolation
                      p_nh(jg)%metrics%vertidx_gradp(2,je,jk,jb) = jk1

                      z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1-1,iblk(je,jb,2))
                      z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1  ,iblk(je,jb,2))
                      z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1+1,iblk(je,jb,2))
                      z3 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1+2,iblk(je,jb,2))

                      coef1 =  z_aux2(je)-z0
                      coef2 = (z_aux2(je)-z0)*(z_aux2(je)-z1)
                      coef3 = (z_aux2(je)-z0)*(z_aux2(je)-z1)*(z_aux2(je)-z2)
                  
                      dn1 = 1._wp/(z0-z1)
                      dn2 = 1._wp/(z1-z2)
                      dn3 = 1._wp/(z2-z3)
                      dn4 = 1._wp/(z0-z2)
                      dn5 = 1._wp/(z0-z3)
                      dn6 = 1._wp/(z1-z3)
                  
                      p_nh(jg)%metrics%coeff_gradp(5,je,jk,jb) =            &
                        1._wp + coef1*dn1 + coef2*dn1*dn4 + coef3*dn1*dn4*dn5
                      p_nh(jg)%metrics%coeff_gradp(6,je,jk,jb) =                               &
                        -(coef1*dn1 + coef2*dn4*(dn1+dn2) + coef3*dn5*(dn1*dn4+dn2*dn4+dn2*dn6))
                      p_nh(jg)%metrics%coeff_gradp(7,je,jk,jb) =          &
                        coef2*dn2*dn4 + coef3*dn5*(dn2*dn4+dn2*dn6+dn3*dn6)
                      p_nh(jg)%metrics%coeff_gradp(8,je,jk,jb) = -coef3*dn3*dn5*dn6

                      jk_start = jk1
                      EXIT

                    ELSE IF (z_aux2(je) <=                                        &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev-1,iblk(je,jb,2))) THEN

                      ! quadratic interpolation/extrapolation
                      p_nh(jg)%metrics%vertidx_gradp(2,je,jk,jb) = nlev - 1

                      z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev-2,iblk(je,jb,2))
                      z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev-1,iblk(je,jb,2))
                      z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev  ,iblk(je,jb,2))

                      coef1 =  z_aux2(je)-z0
                      coef2 = (z_aux2(je)-z0)*(z_aux2(je)-z1)
                  
                      dn1 = 1._wp/(z0-z1)
                      dn2 = 1._wp/(z1-z2)
                      dn4 = 1._wp/(z0-z2)
                  
                      p_nh(jg)%metrics%coeff_gradp(5,je,jk,jb) = &
                        1._wp + coef1*dn1 + coef2*dn1*dn4
                      p_nh(jg)%metrics%coeff_gradp(6,je,jk,jb) = &
                        -(coef1*dn1 + coef2*dn4*(dn1+dn2))
                      p_nh(jg)%metrics%coeff_gradp(7,je,jk,jb) = coef2*dn2*dn4
                      p_nh(jg)%metrics%coeff_gradp(8,je,jk,jb) = 0._wp

                      jk_start = nlev - 2
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO

            ENDIF
          ENDDO
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        ! Generate index list for grid points requiring downward extrapolation of the pressure gradient
        icount_total = SUM(icount(i_startblk:nblks_e))
        p_nh(jg)%metrics%pg_listdim = icount_total
        ic = 0

        ALLOCATE (p_nh(jg)%metrics%pg_edgeidx(icount_total),&
                  p_nh(jg)%metrics%pg_edgeblk(icount_total),&
                  p_nh(jg)%metrics%pg_vertidx(icount_total),&
                  p_nh(jg)%metrics%pg_exdist (icount_total))

        DO jb = i_startblk,nblks_e

          CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                             i_startidx, i_endidx, grf_bdywidth_e+1)

          DO je = i_startidx, i_endidx
            DO jk = flat_idx(je,jb)+1,nlev
              IF (imask(je,jk,jb) == 1) THEN
                ic = ic + 1
                p_nh(jg)%metrics%pg_edgeidx(ic) = je
                p_nh(jg)%metrics%pg_edgeblk(ic) = jb
                p_nh(jg)%metrics%pg_vertidx(ic) = jk
                p_nh(jg)%metrics%pg_exdist (ic) = z_shift(je,jk,jb)
              ENDIF
            ENDDO
          ENDDO
        ENDDO

        DEALLOCATE(imask,icount,z_shift)

      ENDIF

      IF (igradp_method <= 3) THEN
        CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%zdiff_gradp(1,:,:,:))
        CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%zdiff_gradp(2,:,:,:))
      ELSE
        DO ic = 1, 8
          CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%coeff_gradp(ic,:,:,:))
        ENDDO
      ENDIF
      DEALLOCATE(z_me,flat_idx)

    ENDDO

  END SUBROUTINE set_nh_metrics
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !>
  !! Computation of coefficients and index lists needed for truly horizontal
  !! temperature diffusion.
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2010-10-19)
  !!
  SUBROUTINE prepare_zdiffu(p_patch, p_nh, p_int, maxslp, maxhgtd)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_vertical_grid:prepare_zdiffu'

    TYPE(t_patch), TARGET, INTENT(INOUT) :: p_patch
    TYPE(t_nh_state), INTENT(INOUT)      :: p_nh
    TYPE(t_int_state), TARGET,INTENT(IN) :: p_int
    REAL(wp),        INTENT(IN)        :: maxslp(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp),        INTENT(IN)        :: maxhgtd(nproma,p_patch%nlev,p_patch%nblks_c)

    INTEGER :: jk, jb, jc, nblks_c, jk1, jk_start, ji, ji1
    INTEGER :: nlev                  !< number of full levels
    INTEGER :: i_startidx, i_endidx, i_startblk

    INTEGER,  DIMENSION(nproma,p_patch%nblks_c) :: i_masklist, k_start, &
                                                   k_end

    INTEGER, DIMENSION(nproma*p_patch%nblks_c) :: i_indlist, i_blklist
    INTEGER :: i_listdim, numpoints, i_listreduce(p_patch%nblks_c)

    REAL(wp) :: max_nbhgt, z_vintcoeff(nproma,p_patch%nlev,p_patch%nblks_c,3), &
                z_maxslp_avg(nproma,p_patch%nlev,p_patch%nblks_c), z_maxhgtd_avg(nproma,p_patch%nlev,p_patch%nblks_c)

    INTEGER  :: nbidx(nproma,p_patch%nlev,p_patch%nblks_c,3)

    INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk

    nblks_c    = p_patch%nblks_c

    ! number of vertical levels
    nlev   = p_patch%nlev

    ji = 0
    i_masklist(:,:) = 0

    iidx => p_patch%cells%neighbor_idx
    iblk => p_patch%cells%neighbor_blk

    ! Apply cell averaging to slope and height difference fields
    i_startblk = p_patch%cells%start_blk(2,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, i_startidx, i_endidx, 2)

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          z_maxslp_avg(jc,jk,jb) = maxslp(jc,jk,jb)                      *p_int%c_bln_avg(jc,1,jb) + &
                                   maxslp(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_int%c_bln_avg(jc,2,jb) + &
                                   maxslp(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_int%c_bln_avg(jc,3,jb) + &
                                   maxslp(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_int%c_bln_avg(jc,4,jb)

          z_maxhgtd_avg(jc,jk,jb) = maxhgtd(jc,jk,jb)                      *p_int%c_bln_avg(jc,1,jb) + &
                                    maxhgtd(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_int%c_bln_avg(jc,2,jb) + &
                                    maxhgtd(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_int%c_bln_avg(jc,3,jb) + &
                                    maxhgtd(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_int%c_bln_avg(jc,4,jb)

        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! Horizontal index list for truly horizontal diffusion
    i_startblk = p_patch%cells%start_blk(grf_bdywidth_c+1,1)

    ! Attention: this loop is not suitable for OpenMP parallelization
    DO jb = i_startblk,nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, grf_bdywidth_c+1)

      DO jc = i_startidx, i_endidx
        IF ((z_maxslp_avg(jc,nlev,jb)>=thslp_zdiffu .OR. z_maxhgtd_avg(jc,nlev,jb)>=thhgtd_zdiffu) &
            .AND. p_patch%cells%owner_mask(jc,jb)) THEN
          ji = ji+1
          i_blklist(ji) = jb
          i_indlist(ji) = jc
          i_masklist(jc,jb) = 1
        ENDIF
      ENDDO
    ENDDO

    i_listdim = ji
    i_listreduce(:) = 0

!$OMP PARALLEL
    ! Vertical start and end indices for each grid point
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jc, ji, max_nbhgt) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, grf_bdywidth_c+1)

      DO jc = i_startidx, i_endidx
        IF (i_masklist(jc,jb) == 1) THEN
          ! Criterion for setting the end index: Stencil points must not
          ! intersect the ground
          max_nbhgt = MAX(p_nh%metrics%z_mc(iidx(jc,jb,1),nlev,iblk(jc,jb,1)), &
                          p_nh%metrics%z_mc(iidx(jc,jb,2),nlev,iblk(jc,jb,2)), &
                          p_nh%metrics%z_mc(iidx(jc,jb,3),nlev,iblk(jc,jb,3)))
          DO jk = nlev, 1, -1
            IF( p_nh%metrics%z_mc(jc,jk,jb) >= max_nbhgt) THEN
              k_end(jc,jb) = jk
              EXIT
            ENDIF
          ENDDO
          DO jk = 1, nlev
            IF(z_maxslp_avg(jc,jk,jb) >= thslp_zdiffu .OR. z_maxhgtd_avg(jc,jk,jb)>=thhgtd_zdiffu) THEN
              k_start(jc,jb) = jk
              EXIT
            ENDIF
          ENDDO
          ! Reset mask list if vertical index range turns out to be empty
          IF (k_start(jc,jb) > k_end(jc,jb)) THEN
            i_masklist(jc,jb) = 0
            i_listreduce(jb) = i_listreduce(jb) + 1
            k_start(jc,jb) = nlev
          ENDIF
        ENDIF
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    i_listdim = i_listdim - SUM(i_listreduce)

    ! Recompute index list after removal of empty points
    ji = 0
    DO jb = i_startblk,nblks_c
        
      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, grf_bdywidth_c+1)

      DO jc = i_startidx, i_endidx
        IF (i_masklist(jc,jb) == 1) THEN
          ji = ji+1
          i_blklist(ji) = jb
          i_indlist(ji) = jc
        ENDIF
      ENDDO
    ENDDO

!$OMP PARALLEL

    ! Vertical indices for neighbor cells and related vertical interpolation coefficients
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jc, jk1, jk_start) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, grf_bdywidth_c+1)

      DO jc = i_startidx, i_endidx
        IF (i_masklist(jc,jb) == 1) THEN
          jk_start = nlev - 1
          DO jk = k_end(jc,jb), k_start(jc,jb), -1
            DO jk1 = jk_start, 1, -1
              IF ( p_nh%metrics%z_mc(jc,jk,jb) <=                        &
                p_nh%metrics%z_mc(iidx(jc,jb,1),jk1,iblk(jc,jb,1)) .AND. &
                p_nh%metrics%z_mc(jc,jk,jb) >=                           &
                p_nh%metrics%z_mc(iidx(jc,jb,1),jk1+1,iblk(jc,jb,1))) THEN

                ! upper point for linear vertical interpolation
                nbidx(jc,jk,jb,1) = jk1
                ! interpolation weight for upper point
                z_vintcoeff(jc,jk,jb,1) = (p_nh%metrics%z_mc(jc,jk,jb) -  &
                  p_nh%metrics%z_mc(iidx(jc,jb,1),jk1+1,iblk(jc,jb,1))) / &
                  (p_nh%metrics%z_mc(iidx(jc,jb,1),jk1,iblk(jc,jb,1)) -   &
                  p_nh%metrics%z_mc(iidx(jc,jb,1),jk1+1,iblk(jc,jb,1)) )
                jk_start = jk1
                EXIT
              ENDIF
            ENDDO
          ENDDO
          jk_start = nlev - 1
          DO jk = k_end(jc,jb), k_start(jc,jb), -1
            DO jk1 = jk_start, 1, -1
              IF ( p_nh%metrics%z_mc(jc,jk,jb) <=                        &
                p_nh%metrics%z_mc(iidx(jc,jb,2),jk1,iblk(jc,jb,2)) .AND. &
                p_nh%metrics%z_mc(jc,jk,jb) >=                           &
                p_nh%metrics%z_mc(iidx(jc,jb,2),jk1+1,iblk(jc,jb,2))) THEN

                ! upper point for linear vertical interpolation
                nbidx(jc,jk,jb,2) = jk1
                ! interpolation weight for upper point
                z_vintcoeff(jc,jk,jb,2) = (p_nh%metrics%z_mc(jc,jk,jb) -  &
                  p_nh%metrics%z_mc(iidx(jc,jb,2),jk1+1,iblk(jc,jb,2))) / &
                  (p_nh%metrics%z_mc(iidx(jc,jb,2),jk1,iblk(jc,jb,2)) -   &
                  p_nh%metrics%z_mc(iidx(jc,jb,2),jk1+1,iblk(jc,jb,2)) )
                jk_start = jk1
                EXIT
              ENDIF
            ENDDO
          ENDDO
          jk_start = nlev - 1
          DO jk = k_end(jc,jb), k_start(jc,jb), -1
            DO jk1 = jk_start, 1, -1
              IF ( p_nh%metrics%z_mc(jc,jk,jb) <=                        &
                p_nh%metrics%z_mc(iidx(jc,jb,3),jk1,iblk(jc,jb,3)) .AND. &
                p_nh%metrics%z_mc(jc,jk,jb) >=                           &
                p_nh%metrics%z_mc(iidx(jc,jb,3),jk1+1,iblk(jc,jb,3))) THEN

                ! upper point for linear vertical interpolation
                nbidx(jc,jk,jb,3) = jk1
                ! interpolation weight for upper point
                z_vintcoeff(jc,jk,jb,3) = (p_nh%metrics%z_mc(jc,jk,jb) -  &
                  p_nh%metrics%z_mc(iidx(jc,jb,3),jk1+1,iblk(jc,jb,3))) / &
                  (p_nh%metrics%z_mc(iidx(jc,jb,3),jk1,iblk(jc,jb,3)) -   &
                  p_nh%metrics%z_mc(iidx(jc,jb,3),jk1+1,iblk(jc,jb,3)) )
                jk_start = jk1
                EXIT
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! Compute total number of grid points contained in the 3D index
    ! lists to know how much memory to allocate
    numpoints = 0
    DO ji = 1, i_listdim
      jc = i_indlist(ji)
      jb = i_blklist(ji)
      numpoints = numpoints+(k_end(jc,jb)-k_start(jc,jb)+1)
    ENDDO

    ! Now allocate fields
    ALLOCATE(p_nh%metrics%zd_indlist(4,numpoints), &
             p_nh%metrics%zd_blklist(4,numpoints), &
             p_nh%metrics%zd_vertidx(4,numpoints), &
             p_nh%metrics%zd_edgeidx(3,numpoints), &
             p_nh%metrics%zd_edgeblk(3,numpoints), &
             p_nh%metrics%zd_geofac (4,numpoints), &
             p_nh%metrics%zd_e2cell (3,numpoints), &
             p_nh%metrics%zd_intcoef(3,numpoints), &
             p_nh%metrics%zd_diffcoef(numpoints)    )

    p_nh%metrics%zd_listdim = numpoints

    ! Fill index lists
    ji1 = 0
    DO ji = 1, i_listdim
      jc = i_indlist(ji)
      jb = i_blklist(ji)
      DO jk = k_start(jc,jb),k_end(jc,jb)
        ji1 = ji1 + 1
        p_nh%metrics%zd_indlist(1,ji1) = jc
        p_nh%metrics%zd_indlist(2,ji1) = p_patch%cells%neighbor_idx(jc,jb,1)
        p_nh%metrics%zd_indlist(3,ji1) = p_patch%cells%neighbor_idx(jc,jb,2)
        p_nh%metrics%zd_indlist(4,ji1) = p_patch%cells%neighbor_idx(jc,jb,3)
        p_nh%metrics%zd_blklist(1,ji1) = jb
        p_nh%metrics%zd_blklist(2,ji1) = p_patch%cells%neighbor_blk(jc,jb,1)
        p_nh%metrics%zd_blklist(3,ji1) = p_patch%cells%neighbor_blk(jc,jb,2)
        p_nh%metrics%zd_blklist(4,ji1) = p_patch%cells%neighbor_blk(jc,jb,3)
        p_nh%metrics%zd_edgeidx(1,ji1) = p_patch%cells%edge_idx(jc,jb,1)
        p_nh%metrics%zd_edgeidx(2,ji1) = p_patch%cells%edge_idx(jc,jb,2)
        p_nh%metrics%zd_edgeidx(3,ji1) = p_patch%cells%edge_idx(jc,jb,3)
        p_nh%metrics%zd_edgeblk(1,ji1) = p_patch%cells%edge_blk(jc,jb,1)
        p_nh%metrics%zd_edgeblk(2,ji1) = p_patch%cells%edge_blk(jc,jb,2)
        p_nh%metrics%zd_edgeblk(3,ji1) = p_patch%cells%edge_blk(jc,jb,3)
        p_nh%metrics%zd_e2cell(1,ji1)  = p_int%e_bln_c_s(jc,1,jb)
        p_nh%metrics%zd_e2cell(2,ji1)  = p_int%e_bln_c_s(jc,2,jb)
        p_nh%metrics%zd_e2cell(3,ji1)  = p_int%e_bln_c_s(jc,3,jb)
        p_nh%metrics%zd_geofac (1,ji1) = p_int%geofac_n2s(jc,1,jb)
        p_nh%metrics%zd_geofac (2,ji1) = p_int%geofac_n2s(jc,2,jb)
        p_nh%metrics%zd_geofac (3,ji1) = p_int%geofac_n2s(jc,3,jb)
        p_nh%metrics%zd_geofac (4,ji1) = p_int%geofac_n2s(jc,4,jb)
        p_nh%metrics%zd_vertidx(1,ji1) = jk
        p_nh%metrics%zd_vertidx(2,ji1) = nbidx(jc,jk,jb,1)
        p_nh%metrics%zd_vertidx(3,ji1) = nbidx(jc,jk,jb,2)
        p_nh%metrics%zd_vertidx(4,ji1) = nbidx(jc,jk,jb,3)
        p_nh%metrics%zd_intcoef(1,ji1) = z_vintcoeff(jc,jk,jb,1)
        p_nh%metrics%zd_intcoef(2,ji1) = z_vintcoeff(jc,jk,jb,2)
        p_nh%metrics%zd_intcoef(3,ji1) = z_vintcoeff(jc,jk,jb,3)

        p_nh%metrics%zd_diffcoef(ji1)   =                                   &
          MAX(0._wp,SQRT(MAX(0._wp,z_maxslp_avg(jc,jk,jb)-thslp_zdiffu))/250._wp, &
                    2.e-4_wp*SQRT(MAX(0._wp,z_maxhgtd_avg(jc,jk,jb)-thhgtd_zdiffu)) )
        p_nh%metrics%zd_diffcoef(ji1)   = MIN(1._wp/500._wp,p_nh%metrics%zd_diffcoef(ji1))

      ENDDO
    ENDDO

    numpoints = global_sum_array(p_nh%metrics%zd_listdim)

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a,i10)') 'Number of z-diffusion points: ', numpoints
      CALL message(TRIM(routine),message_text)
    ENDIF

  END SUBROUTINE prepare_zdiffu
  !----------------------------------------------------------------------------

END MODULE mo_vertical_grid

