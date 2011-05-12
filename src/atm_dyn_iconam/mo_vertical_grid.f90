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

MODULE mo_vertical_grid

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data,            ONLY: ext_data
  USE mo_model_domain_import, ONLY: n_dom
  USE mo_nonhydrostatic_nml,  ONLY: rayleigh_coeff,damp_height, igradp_method, ivctype,  & 
    &                               vwind_offctr, exner_expol, l_zdiffu_t, thslp_zdiffu, &
    &                               thhgtd_zdiffu, htop_moist_proc, htop_qvadv,          &
    &                               kstart_moist, kstart_qv
  USE mo_sleve_nml,           ONLY: sleve_nml_setup, min_lay_thckn, top_height, decay_scale_1,   &
    &                               decay_scale_2, decay_exp, flat_height, stretch_fac
  USE mo_run_nml,             ONLY: ltestcase, nproma, i_cell_type, msg_level
  USE mo_vertical_coord_table,ONLY: vct_a, vct_b, vct, read_vct
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH, max_dom, &
    &                               min_rlcell_int, min_rlcell
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c, grf_bdywidth_e, grf_fbk_start_c
  USE mo_nh_testcases,        ONLY: layer_thickness, n_flat_level
  USE mo_physical_constants,  ONLY: grav, p0ref, grav, rd, rd_o_cpd, cpd
  USE mo_math_operators,      ONLY: grad_fd_norm, grad_fd_tang, nabla2_scalar
  USE mo_timer,               ONLY: init_timer, cleanup_timer
  USE mo_interpolation,       ONLY: t_int_state, cells2edges_scalar, &
                                    cells2verts_scalar, verts2edges_scalar
  USE mo_math_constants,      ONLY: pi_2, pi
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c
  USE mo_nonhydro_state,      ONLY: t_nh_state
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array, &
                                    sync_patch_array_mult, global_min, global_max
  USE mo_parallel_nml,        ONLY: p_test_run

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  ! Constants used for the computation of the background reference atmosphere of the nh-model
  !
  REAL(wp), PARAMETER :: h_scal_bg = 10000._wp    ! [m]      scale height
  REAL(wp), PARAMETER :: p0sl_bg   = 101325._wp   ! [Pa]     sea level pressure
  REAL(wp), PARAMETER :: t0sl_bg   = 288.15_wp    ! [K]      sea level temperature
  REAL(wp), PARAMETER :: del_t_bg  = 75._wp       ! [K]      difference between sea level
  !                                                          temperature and asymptotic
  !                                                          stratospheric temperature
  REAL(wp), PARAMETER :: grav_o_cpd      = grav/cpd

  INTEGER:: nflat, nflatlev(max_dom), nrdmax(max_dom), nflat_gradp(max_dom)

  PUBLIC :: set_nh_metrics, init_hybrid_coord, init_sleve_coord, &
            nflat, nflatlev, nrdmax, nflat_gradp

  CONTAINS

  !---------------------------------------------------------------------------
  !>
  !! Initialize hybrid coords by reading the 'a' and 'b'. They are assumed to
  !! be HEIGHT BASED in contrast to the hydrostatic model version. The file
  !! name which contains those data has the same name and structure as its
  !! hydrostatic counterpart, namely 'HYB_PARAMS_XX'.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M, (2009-04-14)
  !!
  SUBROUTINE init_hybrid_coord(nlev)

    INTEGER, INTENT(IN) :: nlev  !< number of full levels
    REAL(wp) :: z_height, z_flat
    INTEGER  :: jk, ist
    INTEGER  :: nlevp1            !< number of half levels

    !-------------------------------------------------------------------------

    ! number of vertical half levels
    nlevp1 = nlev+1


    ! If the run is a testrun, the layer thickness might have been given.
    ! If the run is not a testcase run, or the layer thickness is not explicitly
    ! given, the vertical level coordinates are read from HYB_PARAMS_XX.

    IF ( .NOT. ltestcase ) layer_thickness = -999.0_wp

    ! read hybrid parameters as in the hydrostatic model
    IF ( layer_thickness < 0.0_wp) THEN

      CALL read_vct (nlev)

      DO jk = 1, nlevp1
        IF (vct_b(jk) /= 0.0_wp) THEN
          nflat = jk-1
          nflatlev(1) = nflat
          EXIT
        ENDIF
      ENDDO

    ELSE

      ALLOCATE(vct_a(nlevp1), STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish ('mo_vertical_grid:init_hybrid_coord', &
                     'allocation of vct_a failed')
      ENDIF

      ALLOCATE(vct_b(nlevp1), STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish ('mo_vertical_grid:init_hybrid_coord', &
                     'allocation of vct_b failed')
      ENDIF

      ALLOCATE(vct(nlevp1*2), STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish ('mo_vertical_grid:init_hybrid_coord', &
                     'allocation of vct failed')
      ENDIF

      nflat  = -1
      z_flat = REAL(nlevp1-n_flat_level,wp) * layer_thickness
      DO jk = 1, nlevp1
        z_height  = layer_thickness*REAL(nlevp1-jk,wp)
        vct_a(jk) = z_height
        IF ( z_height >= z_flat) THEN
          vct_b(jk) = 0.0_wp
        ELSE
          vct_b(jk) = (z_flat - z_height)/z_flat
          IF (nflat == -1) THEN
            nflat = jk-1
          ENDIF
        ENDIF
      ENDDO
      nflatlev(1) = nflat

      vct(       1:       nlevp1) = vct_a(:)
      vct(nlevp1+1:nlevp1+nlevp1) = vct_b(:)

    ENDIF

    CALL message('mo_vertical_grid: init_hybrid_coord', '')

  END SUBROUTINE init_hybrid_coord

  !---------------------------------------------------------------------------
  !>
  !! Initialize SLEVE coordinate for nonhydrostatic model.
  !! In this initial version, the layer distribution is generated by an analytic
  !! formula based on a couple of namelist variables, but an option to read
  !! in a table may be added.
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2010-07-21)
  !!
  SUBROUTINE init_sleve_coord(nlev)

    INTEGER, INTENT(IN) :: nlev  !< number of full levels
    REAL(wp) :: z_exp
    INTEGER  :: jk, ist
    INTEGER  :: nlevp1        !< number of full and half levels

    !-------------------------------------------------------------------------

    ! number of vertical levels
    nlevp1 = nlev+1


    ! Unlike the hybrid Gal-Chen coordinate, the SLEVE coordinate uses only one
    ! field with coordinate parameters
    ALLOCATE(vct_a(nlevp1), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish ('mo_vertical_grid:init_sleve_coord', &
                   'allocation of vct_a failed')
    ENDIF

    ! Read namelist for SLEVE coordinate
    CALL sleve_nml_setup

    ! However, vct_b also needs to be defined because it is used elsewhere
    ALLOCATE(vct_b(nlevp1), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish ('mo_vertical_grid:init_sleve_coord', &
                   'allocation of vct_b failed')
    ENDIF

    ! And vct is allocated because mo_io_vlist accesses it
    ALLOCATE(vct(nlevp1*2), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish ('mo_vertical_grid:init_hybrid_coord', &
                   'allocation of vct failed')
    ENDIF

    DO jk = 1, nlevp1
      vct_b(jk) = (REAL(jk-1,wp)/REAL(nlev,wp))**1.25_wp
    ENDDO

    z_exp = LOG(min_lay_thckn/top_height)/LOG(2._wp/pi*ACOS(REAL(nlev-1,wp)**stretch_fac/&
      &     REAL(nlev,wp)**stretch_fac))

    ! Set up distribution of coordinate surfaces according to the analytical formula
    ! vct = h_top*(2/pi*arccos(jk-1/nlev))**z_exp (taken from the COSMO model, src_artifdata)
    ! z_exp has been calculated above in order to return min_lay_thckn as thickness
    ! of the lowest model layer
    DO jk = 1, nlevp1
      vct_a(jk)      = top_height*(2._wp/pi*ACOS(REAL(jk-1,wp)**stretch_fac/ &
        &              REAL(nlev,wp)**stretch_fac))**z_exp
      vct(jk)        = vct_a(jk)
      vct(jk+nlevp1) = vct_b(jk)
    ENDDO

    ! Determine nflat (first model level for which coordinate surfaces have a
    ! terrain-following component)
    DO jk = 1, nlevp1
      IF (vct_a(jk) < flat_height) THEN
        nflat = jk-1
        EXIT
      ENDIF
    ENDDO

    ! nflat must not be zero for global domain
    nflat = MAX(1,nflat)
    nflatlev(1) = nflat

    CALL message('mo_vertical_grid: init_sleve_coord', '')

    IF (msg_level >= 15) THEN
     WRITE(message_text,'(a)') 'Heights of coordinate half levels (m):'
        CALL message('', TRIM(message_text))

      DO jk = 1, nlevp1
       WRITE(message_text,'(a,i4,F12.3)') 'jk, vct_a: ',jk, vct_a(jk)
        CALL message('', TRIM(message_text))
      ENDDO
    ENDIF

  END SUBROUTINE init_sleve_coord

  !---------------------------------------------------------------------------
  !>
  !! Compuutes the smoothed topography needed for the SLEVE coordinate.
  !! May be bypassed once an option for reading the smooth topography from data
  !! is available
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2010-07-21)
  !!
  SUBROUTINE compute_smooth_topo(p_patch, p_int)

    TYPE(t_patch),TARGET,INTENT(IN) :: p_patch
    TYPE(t_int_state), INTENT(IN) :: p_int

    INTEGER  :: jg, jb, jc, iter, niter
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx
    REAL(wp) :: z_topo(nproma,1,p_patch%nblks_c),nabla2_topo(nproma,1,p_patch%nblks_c)
    REAL(wp) :: z_topo_v(nproma,1,p_patch%nblks_v)

    !-------------------------------------------------------------------------

    jg = p_patch%id
    niter = 25 ! 25 smoothing iterations (do we need this to be a namelist variable?)

    ! Initialize auxiliary fields for topography with data and nullify nabla2 field
    z_topo(:,1,:)      = ext_data(jg)%atm%topography_c(:,:)
    z_topo_v(:,1,:)    = ext_data(jg)%atm%topography_v(:,:)
    nabla2_topo(:,1,:) = 0._wp

    i_startblk = p_patch%cells%start_blk(2,1)
    nblks_c    = p_patch%nblks_c

    CALL sync_patch_array(SYNC_C,p_patch,z_topo)

    ! Apply nabla2-diffusion niter times to create smooth topography
    DO iter = 1, niter

      CALL nabla2_scalar(z_topo, p_patch, p_int, nabla2_topo, 1, 1)

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 2)

        DO jc = i_startidx, i_endidx
          z_topo(jc,1,jb) = z_topo(jc,1,jb) + 0.125_wp*nabla2_topo(jc,1,jb) &
            &                               * p_patch%cells%area(jc,jb)
        ENDDO
      ENDDO

      CALL sync_patch_array(SYNC_C,p_patch,z_topo)

    ENDDO

    ! Interpolate smooth topography from cells to vertices
    CALL cells2verts_scalar(z_topo,p_patch,p_int%cells_aw_verts,z_topo_v,1,1)

    CALL sync_patch_array(SYNC_V,p_patch,z_topo_v)

    ext_data(jg)%atm%topography_smt_c(:,:) = z_topo(:,1,:)
    ext_data(jg)%atm%topography_smt_v(:,:) = z_topo_v(:,1,:)

  END SUBROUTINE compute_smooth_topo

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
  SUBROUTINE set_nh_metrics(p_patch, p_nh, p_int)

    TYPE(t_patch), TARGET, INTENT(INOUT) :: p_patch(n_dom)  !< patch
    TYPE(t_nh_state), INTENT(INOUT)      :: p_nh(n_dom)
    TYPE(t_int_state),  INTENT(IN)       :: p_int(n_dom)

    INTEGER :: jg, jgp, jk, jk1, jk_start, jb, jc, je, jv, jn, jgc, nlen, &
               nblks_c, npromz_c, nblks_e, npromz_e, nblks_v, npromz_v, ic
    INTEGER :: nlev, nlevp1              !< number of full levels
    INTEGER :: nshift_total(n_dom)       !< Total shift of model top w.r.t. global domain
    ! Note: the present way of setting up coordinate surfaces will not work for vertical refinement
    INTEGER :: i_startidx, i_endidx, i_startblk, i_endblk, i_nchdom, icount_total
    INTEGER :: ica(max_dom)

    REAL(wp) :: z_diff, z1, z2, z3, z_help(nproma), z_maxhgt(p_patch(1)%nlevp1), &
      &         z_temp(nproma), z_rho(nproma,p_patch(1)%nlev), z_aux1(nproma),   &
      &         z_aux2(nproma)
    REAL(wp) :: z_fac1, z_fac2, z_topo_dev(nproma), z_maxslope, z_maxhdiff, z_offctr
    REAL(wp), ALLOCATABLE :: z_ifv(:,:,:), z_mfv(:,:,:)
    REAL(wp), ALLOCATABLE :: z_me(:,:,:),z_maxslp(:,:,:),z_maxhgtd(:,:,:),z_shift(:,:,:)
    REAL(wp) :: extrapol_dist
    INTEGER,  ALLOCATABLE :: flat_idx(:,:), imask(:,:,:),icount(:)
    INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
    LOGICAL :: l_half_lev_centr

    !------------------------------------------------------------------------

    SELECT CASE (i_cell_type)
    CASE (6)
      l_half_lev_centr = .TRUE.
      ! The HALF LEVEL where the model layer are flat, moves one layer upward.
      ! there could also be a zero there
      nflat = nflat-1
    CASE (3)
      l_half_lev_centr = .FALSE.
    END SELECT

    ! model timer is not yet initialized here and grad_fd_norm asks for a
    ! running timer
    CALL init_timer

    nshift_total(1) = 0

    DO jg = 1,n_dom

      nblks_c   = p_patch(jg)%nblks_int_c
      npromz_c  = p_patch(jg)%npromz_int_c

      i_nchdom   = MAX(1,p_patch(jg)%n_childdom)

      ! number of vertical levels
      nlev   = p_patch(jg)%nlev
      nlevp1 = p_patch(jg)%nlevp1

      ! total shift of model top with respect to global domain
      IF (jg > 1) THEN
        jgp = p_patch(jg)%parent_id
        nshift_total(jg) = nshift_total(jgp) + p_patch(jg)%nshift
        nflatlev(jg)     = nflatlev(1) - nshift_total(jg)
      ENDIF
      IF (i_cell_type == 6) nflatlev(jg) = nflat
      IF (jg > 1 .AND. nshift_total(jg) > 0 .AND. nflatlev(jg) < 1) THEN
        CALL finish ('mo_vertical_grid:set_nh_metrics', &
                     'nflat must be more than the top of the innermost nested domain')
      ENDIF

      ! Compute smooth topography when SLEVE coordinate is used
      IF (ivctype == 2) THEN
        CALL compute_smooth_topo(p_patch(jg), p_int(jg))
      ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, jk1, z_fac1, z_fac2, z_topo_dev)
      DO jb = 1,nblks_c

        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
           ! To avoid uninitialized field elements going in the output routine
           p_nh(jg)%metrics%geopot(nlen+1:nproma,:,jb) = 0._wp
           p_nh(jg)%metrics%z_ifc (nlen+1:nproma,:,jb) = 0._wp
        ENDIF

        p_nh(jg)%metrics%z_ifc(1:nlen,nlevp1,jb) = ext_data(jg)%atm%topography_c(1:nlen,jb)

        ! vertical interface height
        IF (ivctype == 1) THEN ! hybrid Gal-Chen coordinate
          DO jk = 1, nlev
            jk1 = jk + nshift_total(jg)
            p_nh(jg)%metrics%z_ifc(1:nlen,jk,jb) = &
            & vct_a(jk1) + vct_b(jk1)*p_nh(jg)%metrics%z_ifc(1:nlen,nlevp1,jb)
          ENDDO
        ELSE IF (ivctype == 2) THEN ! SLEVE coordinate (Leuenberger et al. MWR 2010)
          DO jk = 1, nflatlev(jg)
            jk1 = jk + nshift_total(jg)
            p_nh(jg)%metrics%z_ifc(1:nlen,jk,jb) = vct_a(jk1)
          ENDDO
          DO jk = nflatlev(jg) + 1, nlev
            jk1 = jk + nshift_total(jg)
            ! Scaling factors for large-scale and small-scale topography
            z_fac1 = SINH((top_height/decay_scale_1)**decay_exp - &
              &           (vct_a(jk1)/decay_scale_1)**decay_exp)/ &
              &      SINH((top_height/decay_scale_1)**decay_exp)
            z_fac2 = SINH((top_height/decay_scale_2)**decay_exp - &
              &           (vct_a(jk1)/decay_scale_2)**decay_exp)/ &
              &      SINH((top_height/decay_scale_2)**decay_exp)

            ! Small-scale topography (i.e. full topo - smooth topo)
            z_topo_dev(1:nlen) = ext_data(jg)%atm%topography_c(1:nlen,jb) - &
              &                  ext_data(jg)%atm%topography_smt_c(1:nlen,jb)

            p_nh(jg)%metrics%z_ifc(1:nlen,jk,jb) = vct_a(jk1) +       &
              & ext_data(jg)%atm%topography_smt_c(1:nlen,jb)*z_fac1 + &
              & z_topo_dev(1:nlen)*z_fac2
          ENDDO
          ! Ensure that layers do not intersect; except for the surface layer, the interface level
          ! distance is limited to min_lay_thckn
          DO jk = nlev-1, 1, -1
            DO jc = 1, nlen
              p_nh(jg)%metrics%z_ifc(jc,jk,jb) = MAX(p_nh(jg)%metrics%z_ifc(jc,jk,jb), &
                p_nh(jg)%metrics%z_ifc(jc,jk+1,jb)+min_lay_thckn)
            ENDDO
          ENDDO
        ENDIF

        DO jk = 1, nlevp1
          ! geopotential (at interfaces)
          p_nh(jg)%metrics%geopot_ifc(1:nlen,jk,jb) = grav * &
          &                p_nh(jg)%metrics%z_ifc(1:nlen,jk,jb)
        ENDDO

        DO jk = 1, nlev
          ! geometric height of full levels
          p_nh(jg)%metrics%z_mc(1:nlen,jk,jb) = 0.5_wp    &
            &  * (p_nh(jg)%metrics%z_ifc(1:nlen,jk  ,jb)  &
            &  +  p_nh(jg)%metrics%z_ifc(1:nlen,jk+1,jb))
          ! geopotential on full levels
          p_nh(jg)%metrics%geopot(1:nlen,jk,jb) = grav &
          & * p_nh(jg)%metrics%z_mc(1:nlen,jk,jb)
    
        ENDDO

        !-------------------------------------------------------------
        !-------------------------------------------------------------
        IF (l_half_lev_centr) THEN
          ! redefine half levels to be in the center between full levels
          DO jk = 2, nlev
            p_nh(jg)%metrics%z_ifc(1:nlen,jk,jb) = 0.5_wp * &
            & (p_nh(jg)%metrics%z_mc(1:nlen,jk-1,jb) &
            & +p_nh(jg)%metrics%z_mc(1:nlen,jk  ,jb))
          ENDDO
        ENDIF
        !-------------------------------------------------------------
        !-------------------------------------------------------------

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
          p_nh(jg)%metrics%diff_1_o_dz(1:nlen,jk,jb) = &
          &  1.0_wp/p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk-1,jb) & 
          & -1.0_wp/p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk  ,jb)
          p_nh(jg)%metrics%mult_1_o_dz(1:nlen,jk,jb) = &
          &  1.0_wp/(p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk-1,jb) & 
          &         *p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk  ,jb))
        ENDDO
        p_nh(jg)%metrics%diff_1_o_dz(1:nlen,1,jb) = 0.0_wp
        p_nh(jg)%metrics%mult_1_o_dz(1:nlen,1,jb) = 0.0_wp
        p_nh(jg)%metrics%ddqz_z_half(1:nlen,1,jb) =   &
        & 2.0_wp*(p_nh(jg)%metrics%z_ifc(1:nlen,1,jb) &
        &       - p_nh(jg)%metrics%z_mc (1:nlen,1,jb))
        p_nh(jg)%metrics%diff_1_o_dz(1:nlen,nlevp1,jb) = 0.0_wp
        p_nh(jg)%metrics%mult_1_o_dz(1:nlen,nlevp1,jb) = 0.0_wp
        p_nh(jg)%metrics%ddqz_z_half(1:nlen,nlevp1,jb) =    &
        & 2.0_wp*(p_nh(jg)%metrics%z_mc (1:nlen,nlev  ,jb)  &
        &       - p_nh(jg)%metrics%z_ifc(1:nlen,nlevp1,jb))
        IF (i_cell_type==3) THEN
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
            p_nh(jg)%metrics%geopot_agl(1:nlen,jk-1,jb) =  0.5 *              &
           &              ( p_nh(jg)%metrics%geopot_agl_ifc(1:nlen,jk-1  ,jb )&
           &              + p_nh(jg)%metrics%geopot_agl_ifc(1:nlen,jk    ,jb))
          ENDDO
      ENDDO
!$OMP END DO
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
        CALL message('', TRIM(message_text))
      ENDDO
    ENDIF

      ! For the tangential slope we need temporarily also the height
      ! at the vertices
      ALLOCATE(z_ifv(nproma,nlevp1,p_patch(jg)%nblks_v))
      ALLOCATE(z_mfv(nproma,nlev  ,p_patch(jg)%nblks_v))
      nblks_v   = p_patch(jg)%nblks_int_v
      npromz_v  = p_patch(jg)%npromz_int_v

      ! Start index for slope computations
      i_startblk = p_patch(jg)%edges%start_blk(2,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, jk1, z_fac1, z_fac2, z_topo_dev)
      DO jb = 1,nblks_v
        IF (jb /= nblks_v) THEN
           nlen = nproma
        ELSE
           nlen = npromz_v
        ENDIF
        z_ifv(1:nlen,nlevp1,jb) = ext_data(jg)%atm%topography_v(1:nlen,jb)
        ! vertical interface height at vertices
        IF (ivctype == 1) THEN ! hybrid Gal-Chen coordinate
          DO jk = 1, nlev
            jk1 = jk + nshift_total(jg)
            z_ifv(1:nlen,jk,jb) = vct_a(jk1) + vct_b(jk1)*z_ifv(1:nlen,nlevp1,jb)
          ENDDO
        ELSE IF (ivctype == 2) THEN ! SLEVE coordinate (Leuenberger et al. MWR 2010)

          DO jk = 1, nflatlev(jg)
            jk1 = jk + nshift_total(jg)
            z_ifv(1:nlen,jk,jb) = vct_a(jk1)
          ENDDO
          DO jk = nflatlev(jg) + 1, nlev
            jk1 = jk + nshift_total(jg)
            ! Scaling factors for large-scale and small-scale topography
            z_fac1 = SINH((top_height/decay_scale_1)**decay_exp - &
              &           (vct_a(jk1)/decay_scale_1)**decay_exp)/ &
              &      SINH((top_height/decay_scale_1)**decay_exp)
            z_fac2 = SINH((top_height/decay_scale_2)**decay_exp - &
              &           (vct_a(jk1)/decay_scale_2)**decay_exp)/ &
              &      SINH((top_height/decay_scale_2)**decay_exp)

            ! Small-scale topography (i.e. full topo - smooth topo)
            z_topo_dev(1:nlen) = ext_data(jg)%atm%topography_v(1:nlen,jb) - &
              &                  ext_data(jg)%atm%topography_smt_v(1:nlen,jb)

            z_ifv(1:nlen,jk,jb) = vct_a(jk1) +                        &
              & ext_data(jg)%atm%topography_smt_v(1:nlen,jb)*z_fac1 + &
              & z_topo_dev(1:nlen)*z_fac2
          ENDDO
        ENDIF

        IF (l_half_lev_centr) THEN
          DO jk = 1, nlev
            ! full level height
            z_mfv(1:nlen,jk,jb)=0.5_wp*(z_ifv(1:nlen,jk,jb)+z_ifv(1:nlen,jk+1,jb))
          ENDDO
          ! redefine half levels to be in the center between full levels
          DO jk = 2, nlev
            z_ifv(1:nlen,jk,jb) = 0.5_wp*(z_mfv(1:nlen,jk-1,jb)+z_mfv(1:nlen,jk,jb))
          ENDDO
        ENDIF
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      p_nh(jg)%metrics%ddxt_z_half = 0._wp
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

      IF (p_test_run) p_nh(jg)%metrics%z_mc_e = 0._wp

      CALL cells2edges_scalar(p_nh(jg)%metrics%z_mc, &
                              p_patch(jg),p_int(jg)%c_lin_e, &
                              p_nh(jg)%metrics%z_mc_e )

      CALL sync_patch_array_mult(SYNC_E,p_patch(jg),4,p_nh(jg)%metrics%ddxt_z_half,        &
                                 p_nh(jg)%metrics%ddxn_z_half,p_nh(jg)%metrics%ddxn_z_full,&
                                 p_nh(jg)%metrics%z_mc_e)

      ! vertically averaged metrics
      nblks_e   = p_patch(jg)%nblks_int_e
      npromz_e  = p_patch(jg)%npromz_int_e
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, je)
      DO jb = i_startblk,nblks_e

        CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, 2)

        IF (.NOT.l_half_lev_centr) THEN
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
              p_nh(jg)%metrics%ddxn_z_full(je,jk,jb) = 0.5_wp * &
              & (p_nh(jg)%metrics%ddxn_z_half(je,jk  ,jb) + &
              &  p_nh(jg)%metrics%ddxn_z_half(je,jk+1,jb))
            ENDDO
          ENDDO
        ENDIF

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      IF (p_test_run) p_nh(jg)%metrics%ddqz_z_full_e = 0._wp

      ! functional determinant at full level edges
      CALL cells2edges_scalar(p_nh(jg)%metrics%ddqz_z_full, &
           &                  p_patch(jg), p_int(jg)%c_lin_e, &
           &                  p_nh(jg)%metrics%ddqz_z_full_e )

      CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%ddqz_z_full_e)

      IF(i_cell_type==6)THEN

        IF (p_test_run) THEN
          p_nh(jg)%metrics%ddqz_z_half_e = 0._wp
          p_nh(jg)%metrics%ddqz_z_half_v = 0._wp
          p_nh(jg)%metrics%ddqz_z_full_v = 0._wp
          p_nh(jg)%metrics%ddqz_z_half_r = 0._wp
        ENDIF
        ! functional determinant at half level the edges
        CALL cells2edges_scalar(p_nh(jg)%metrics%ddqz_z_half, &
                                p_patch(jg),p_int(jg)%c_lin_e, &
                                p_nh(jg)%metrics%ddqz_z_half_e, 1, nlevp1 )

        CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%ddqz_z_half_e)

        ! functional determinant at half level dual centers
        CALL cells2verts_scalar(p_nh(jg)%metrics%ddqz_z_half, &
                                p_patch(jg),p_int(jg)%cells_aw_verts, &
                                p_nh(jg)%metrics%ddqz_z_half_v, 1, nlevp1)

        ! functional determinant at full level dual centers
        CALL cells2verts_scalar(p_nh(jg)%metrics%ddqz_z_full, &
                                p_patch(jg),p_int(jg)%cells_aw_verts, &
                                p_nh(jg)%metrics%ddqz_z_full_v, 1, nlev)

        CALL sync_patch_array_mult(SYNC_V,p_patch(jg),2,p_nh(jg)%metrics%ddqz_z_half_v, &
                                   p_nh(jg)%metrics%ddqz_z_full_v)

        CALL verts2edges_scalar(p_nh(jg)%metrics%ddqz_z_half_v, &
                                p_patch(jg), p_int(jg)%tria_aw_rhom, &
                                p_nh(jg)%metrics%ddqz_z_half_r, 1, nlevp1)

        CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%ddqz_z_half_r)

        iidx => p_patch(jg)%verts%edge_idx
        iblk => p_patch(jg)%verts%edge_blk
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, jv)
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
!$OMP END DO
!$OMP END PARALLEL
      ENDIF

      ! Determine start level for moist physics processes (specified by htop_moist_proc)
      DO jk = 1, nlev
        jk1 = jk + nshift_total(jg)
        IF (0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < htop_moist_proc) THEN
          kstart_moist(jg) = jk
          EXIT
        ENDIF
      ENDDO

      IF (kstart_moist(jg) > 1 .AND. msg_level >= 10) THEN
        WRITE(message_text,'(2(a,i4))') 'Domain', jg, &
          '; computation of moist physics processes starts in layer ', kstart_moist(jg)
        CALL message('mo_vertical_grid',message_text)
      ENDIF

      ! Determine start level for QV advection (specified by htop_qvadv)
      DO jk = 1, nlev
        jk1 = jk + nshift_total(jg)
        IF (0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < htop_qvadv) THEN
          kstart_qv(jg) = jk
          EXIT
        ENDIF
      ENDDO

      IF (kstart_qv(jg) > 1 .AND. msg_level >= 10) THEN
        WRITE(message_text,'(2(a,i4))') 'Domain', jg, &
          '; computation of QV advection starts in layer ', kstart_qv(jg)
        CALL message('mo_vertical_grid',message_text)
      ENDIF

      ! Rayleigh damping properties
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, jc, z_diff)
      DO jb = 1,nblks_c
        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF
        DO jk = 1, nlevp1
          DO jc = 1, nlen
            z_diff = MAX(0.0_wp,p_nh(jg)%metrics%z_ifc(jc,jk,jb)-damp_height(jg))
            p_nh(jg)%metrics%rayleigh_w(jc,jk,jb)= &
            &   rayleigh_coeff(jg) * (SIN ( pi_2 * z_diff &
            &   /MAX(1.e-3_wp,p_nh(jg)%metrics%z_ifc(jc,1,jb)-damp_height(jg))))**2
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO

      ! offcentering in vertical mass flux 
!$OMP WORKSHARE
      p_nh(jg)%metrics%vwind_expl_wgt(:,:) = 0.5_wp - vwind_offctr
      p_nh(jg)%metrics%vwind_impl_wgt(:,:) = 0.5_wp + vwind_offctr
!$OMP END WORKSHARE

!$OMP END PARALLEL

      ! Determine depth of damping layer
      DO jk = 1, nlevp1
        z_maxhgt(jk) =MAXVAL(p_nh(jg)%metrics%z_ifc(:,jk,:),MASK=p_patch(jg)%cells%owner_mask(:,:))
      ENDDO

      ! Get maximum height of all PEs
      z_maxhgt = global_max(z_maxhgt)

      nrdmax(jg) = 0
      DO jk = 2, nlevp1
        IF (z_maxhgt(jk) >= damp_height(jg)) nrdmax(jg) = jk-1
      ENDDO

      IF (msg_level >= 10) THEN
        WRITE(message_text,'(2(a,i4))') 'Domain', jg, &
          '; end index of Rayleigh damping layer: ', nrdmax(jg)
        CALL message('mo_vertical_grid',message_text)
      ENDIF

      ! Compute variable Exner extrapolation factors and offcentering coefficients for the
      ! vertical wind solver to optimize numerical stability over steep orography
      iidx => p_patch(jg)%cells%edge_idx
      iblk => p_patch(jg)%cells%edge_blk

      i_startblk = p_patch(jg)%cells%start_blk(2,1)

      IF (i_cell_type == 3) THEN
        ALLOCATE (z_maxslp(nproma,nlev,p_patch(jg)%nblks_c), &
                  z_maxhgtd(nproma,nlev,p_patch(jg)%nblks_c) )

!$OMP PARALLEL
!$OMP WORKSHARE
        ! Initialization to ensure that values are properly set at lateral boundaries
        p_nh(jg)%metrics%exner_exfac(:,:,:) = exner_expol

        z_maxslp(:,:,:)  = 0._wp
        z_maxhgtd(:,:,:) = 0._wp
!$OMP END WORKSHARE

!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jc, z_maxslope, z_offctr)
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
              IF (z_maxslp(jc,jk,jb) > 1._wp) &
                p_nh(jg)%metrics%exner_exfac(jc,jk,jb) = &
                  MAX(-0.25_wp,0.50_wp*(1._wp-z_maxslp(jc,jk,jb)))
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
            z_offctr =   MIN(0.5_wp,z_offctr)

            p_nh(jg)%metrics%vwind_expl_wgt(jc,jb) = 0.5_wp - z_offctr
            p_nh(jg)%metrics%vwind_impl_wgt(jc,jb) = 0.5_wp + z_offctr

          ENDDO
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ENDIF

      ! The remaining computations are needed for the triangular grid only
      ! once the initialization in mo_nh_testcases is properly rewritten for hexagons
      IF (i_cell_type == 6) CYCLE

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
!$OMP DO PRIVATE(jb, nlen, jk, jc, z1, z2, z3)
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
!$OMP END DO
!$OMP END PARALLEL

      IF (i_cell_type == 3 .AND. msg_level >= 10) THEN
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
        CALL message('mo_vertical_grid',message_text)
        WRITE(message_text,'(a,f8.4)') 'Maximum slope: ', z_maxslope
        CALL message('mo_vertical_grid',message_text)
        WRITE(message_text,'(a,f8.1)') 'Maximum height difference between adjacent points: ', &
          z_maxhdiff
        CALL message('mo_vertical_grid',message_text)
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

      IF (i_cell_type == 3) THEN

      ! Reference atmosphere fields for triangular code
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, z_help, z_temp, z_aux1, z_aux2, z_rho)
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

            ! Reference temperature, full level mass points
            z_temp(1:nlen) = (t0sl_bg-del_t_bg)+del_t_bg*        &
              & EXP(-p_nh(jg)%metrics%z_mc(1:nlen,jk,jb)/h_scal_bg)

            ! Reference density, full level mass points
            z_rho(1:nlen,jk) = z_aux1(1:nlen)/(rd*z_temp(1:nlen))

            ! Reference Potential temperature, full level mass points
            p_nh(jg)%metrics%theta_ref_mc(1:nlen,jk,jb) = z_temp(1:nlen) &
              & /p_nh(jg)%metrics%exner_ref_mc(1:nlen,jk,jb)

            ! First vertical derivative of reference Exner pressure, full level mass points
            ! Note: for computational efficiency, this field is in addition divided by
            ! the vertical layer thickness
            p_nh(jg)%metrics%d_exner_dz_ref_mc(1:nlen,jk,jb)   =       &
              & -grav/cpd/p_nh(jg)%metrics%theta_ref_mc(1:nlen,jk,jb)* &
              & p_nh(jg)%metrics%inv_ddqz_z_full(1:nlen,jk,jb)

            ! Second vertical derivative of reference Exner pressure, full level mass points
            p_nh(jg)%metrics%d2_exner_dz2_ref_mc(1:nlen,jk,jb)   =                                &
              &  grav/cpd/p_nh(jg)%metrics%theta_ref_mc(1:nlen,jk,jb)**2                          &
              & *(grav/cpd-del_t_bg/h_scal_bg*EXP(-p_nh(jg)%metrics%z_mc(1:nlen,jk,jb)/h_scal_bg))&
              & /p_nh(jg)%metrics%exner_ref_mc(1:nlen,jk,jb)
          ENDDO

          ! rho_refcorr_ic is not needed at surface level => set to zero
          p_nh(jg)%metrics%rho_refcorr_ic(:,nlevp1,jb)    = 0._wp

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

            ! Correction term for vertical interpolation of density, half level mass points
            ! This (additive) term summarizes the following steps:
            ! - subtract reference density from full field at full levels
            ! - perform bilinear vertical interpolation to half levels
            ! - add reference density at half levels
            IF (jk == 1) THEN
              p_nh(jg)%metrics%rho_refcorr_ic(1:nlen,jk,jb) = z_aux2(1:nlen) - &
                & (p_nh(jg)%metrics%wgtfacq1_c(1:nlen,1,jb)*z_rho(1:nlen,1)    &
                & +p_nh(jg)%metrics%wgtfacq1_c(1:nlen,2,jb)*z_rho(1:nlen,2)    &
                & +p_nh(jg)%metrics%wgtfacq1_c(1:nlen,3,jb)*z_rho(1:nlen,3) )
            ELSE IF (jk <= nlev) THEN
              p_nh(jg)%metrics%rho_refcorr_ic(1:nlen,jk,jb)   =                               &
                & z_aux2(1:nlen) - (p_nh(jg)%metrics%wgtfac_c(1:nlen,jk,jb)*z_rho(1:nlen,jk)  &
                & +(1._wp-p_nh(jg)%metrics%wgtfac_c(1:nlen,jk,jb))*z_rho(1:nlen,jk-1))
            ENDIF

          ENDDO
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ENDIF

      ! Compute information needed for Taylor-expansion-based pressure gradient calculation
      IF (igradp_method == 1) THEN
        nflat_gradp(jg) = nlev
        CYCLE
      ENDIF
      ALLOCATE(z_me(nproma,nlev,p_patch(jg)%nblks_e), flat_idx(nproma,p_patch(jg)%nblks_e))
      IF (p_test_run) THEN
        z_me = 0._wp
        p_nh(jg)%metrics%zdiff_gradp = 0._wp
      ENDIF
      flat_idx = nlev

      ! Compute geometric height at edge points
      CALL cells2edges_scalar(p_nh(jg)%metrics%z_mc, p_patch(jg), &
             p_int(jg)%c_lin_e, z_me)

      CALL sync_patch_array(SYNC_E,p_patch(jg),z_me)

      iidx => p_patch(jg)%edges%cell_idx
      iblk => p_patch(jg)%edges%cell_blk

      i_startblk = p_patch(jg)%edges%start_blk(2,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, je)
      DO jb = i_startblk,nblks_e

        CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

            p_nh(jg)%metrics%vertidx_gradp(1:2,je,jk,jb) = jk
            p_nh(jg)%metrics%zdiff_gradp(1,je,jk,jb)     =  &
              z_me(je,jk,jb) - p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk,iblk(je,jb,1))
            p_nh(jg)%metrics%zdiff_gradp(2,je,jk,jb)     =  &
              z_me(je,jk,jb) - p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk,iblk(je,jb,2))

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
!$OMP END DO
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
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jk1, jk_start, je)
      DO jb = i_startblk,nblks_e

        CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, 2)

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
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      IF (igradp_method == 3) THEN

        i_startblk = p_patch(jg)%edges%start_blk(grf_bdywidth_e+1,1)

        ALLOCATE(imask(nproma,nlev,p_patch(jg)%nblks_e),icount(p_patch(jg)%nblks_e), &
                 z_shift(nproma,nlev,p_patch(jg)%nblks_e) )

        extrapol_dist = 20._wp ! maximum allowed extrapolation distance; may become a namelist variable later on

        ! Recompute indices and height differences if truly horizontal pressure gradient 
        ! computation would intersect the ground
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jk1, jk_start, je, z_aux1, z_aux2)
        DO jb = i_startblk,nblks_e

          CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                             i_startidx, i_endidx, grf_bdywidth_e+1)

          imask(:,:,jb) = 0
          icount(jb)    = 0

          DO je = i_startidx, i_endidx

            z_aux1(je) = MAX(p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),nlevp1,iblk(je,jb,1)), &
                             p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),nlevp1,iblk(je,jb,2)))

            z_aux2(je) = z_aux1(je) - extrapol_dist ! allow for some limited downward extrapolation

            jk_start = flat_idx(je,jb)
            DO jk = flat_idx(je,jb)+1,nlev
              IF ( z_me(je,jk,jb) < z_aux2(je)) THEN

                ! Save information needed for index list setup
                IF (p_patch(jg)%edges%owner_mask(je,jb)) THEN
                  imask(je,jk,jb) = 1
                  icount(jb)      = icount(jb) + 1
                  z_shift(je,jk,jb) = z_me(je,jk,jb) - z_aux2(je)
                ENDIF

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

          ENDDO
        ENDDO
!$OMP END DO
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

      CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%zdiff_gradp(1,:,:,:))
      CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%zdiff_gradp(2,:,:,:))

      DEALLOCATE(z_me,flat_idx)

    ENDDO

    CALL cleanup_timer

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

    REAL(wp) :: max_nbhgt, z_vintcoeff(nproma,p_patch%nlev,p_patch%nblks_c,3)
    INTEGER  :: nbidx(nproma,p_patch%nlev,p_patch%nblks_c,3)

    INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk

    nblks_c    = p_patch%nblks_c

    ! number of vertical levels
    nlev   = p_patch%nlev

    ji = 0
    i_masklist(:,:) = 0

    iidx => p_patch%cells%neighbor_idx
    iblk => p_patch%cells%neighbor_blk

    ! Horizontal index list for truly horizontal diffusion
    i_startblk = p_patch%cells%start_blk(grf_bdywidth_c+1,1)

    ! Attention: this loop is not suitable for OpenMP parallelization
    DO jb = i_startblk,nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, grf_bdywidth_c+1)

      DO jc = i_startidx, i_endidx
        IF ((maxslp(jc,nlev,jb)>=thslp_zdiffu .OR. maxhgtd(jc,nlev,jb)>=thhgtd_zdiffu) &
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
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jc, ji, max_nbhgt)
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
            IF(maxslp(jc,jk,jb) >= thslp_zdiffu .OR. maxhgtd(jc,jk,jb)>=thhgtd_zdiffu) THEN
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
!$OMP END DO
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
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jc, jk1, jk_start)
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
!$OMP END DO
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
          MAX(0._wp,SQRT(MAX(0._wp,maxslp(jc,jk,jb)-thslp_zdiffu))/250._wp, &
                    2.e-4_wp*SQRT(MAX(0._wp,maxhgtd(jc,jk,jb)-thhgtd_zdiffu)) )
        p_nh%metrics%zd_diffcoef(ji1)   = MIN(1._wp/200._wp,p_nh%metrics%zd_diffcoef(ji1))

      ENDDO
    ENDDO

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a,i10)') 'Number of z-diffusion points: ', numpoints
      CALL message('mo_vertical_grid:prepare_zdiffu',message_text)
    ENDIF

  END SUBROUTINE prepare_zdiffu
  !----------------------------------------------------------------------------

END MODULE mo_vertical_grid
