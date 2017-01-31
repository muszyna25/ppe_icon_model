!>
!!  This module contains parameters, initialization subroutines and.
!!
!!  This module contains parameters, initialization subroutines and
!!  functions to be used in the Pure 3D-Advection test case of the
!!  hydrostatic dynamical core.
!!
!! @par Revision History
!!  Developed by Jochen Foerstner, DWD (2008-05-26)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_pa_test
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2008
!
!-------------------------------------------------------------------------
!
!
!

USE mo_kind,                ONLY: wp
USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
USE mo_physical_constants,  ONLY: rgrav, rd
USE mo_math_constants,      ONLY: pi_2, pi
USE mo_vertical_coord_table,ONLY: vct_a, vct_b, ceta, cetah
USE mo_eta_coord_diag,      ONLY: half_level_pressure, full_level_pressure
USE mo_model_domain,        ONLY: t_patch
USE mo_ext_data_types,      ONLY: t_external_data
USE mo_icoham_dyn_types,    ONLY: t_hydro_atm_prog, t_hydro_atm_diag
USE mo_intp_data_strc,      ONLY: t_int_state
USE mo_intp,                ONLY: cells2edges_scalar
USE mo_parallel_config,     ONLY: nproma
USE mo_run_config,          ONLY: nsteps, ntracer
USE mo_ncar_testcases,      ONLY: init_pure_adv_wind, init_pure_adv_tracers, &
  & init_ncar_testcases_domain

USE mo_io_units,            ONLY: find_next_free_unit
USE mo_exception,           ONLY: finish
USE mo_mpi,                 ONLY: my_process_is_stdio
USE mo_datetime,            ONLY: rdaylen
USE mo_grid_config,         ONLY: grid_sphere_radius, grid_angular_velocity

IMPLICIT NONE

PRIVATE

REAL(wp), PARAMETER :: zp0     = 100000._wp
REAL(wp), PARAMETER :: zt0     = 300._wp
REAL(wp), PARAMETER :: ztau    = 4._wp * rdaylen     ! 4 days period
REAL(wp), PARAMETER :: zomega0 = pi*40000._wp / ztau

PUBLIC :: init_hydro_state_prog_patest,  &
  &       set_vertical_velocity

!--------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!
  !>
  !!               Initialization of prognostic state vector.
  !!
  !!
  !! @par Revision History
  !!  Original version by Jochen Foerstner, DWD (2008-05)
  !!  Code restructuring by Almut Gassmann, MPI-M (2008-10)
  !!
  SUBROUTINE init_hydro_state_prog_patest( ptr_patch, ptr_prog, ptr_diag,      &
    &                                      ptr_int, ptr_ext_data,              &
    &                                      p_rotate_axis_deg, linit_tracer_fv, &
    &                                      tracer_inidist_list)

    TYPE(t_patch),TARGET,INTENT(INOUT):: ptr_patch
    TYPE(t_int_state), INTENT(INOUT)  :: ptr_int
    TYPE(t_hydro_atm_prog), INTENT(INOUT) :: ptr_prog
    TYPE(t_hydro_atm_diag), INTENT(INOUT) :: ptr_diag
    TYPE(t_external_data), INTENT(INOUT) :: ptr_ext_data !< external data
    REAL(wp), INTENT(IN)            :: p_rotate_axis_deg
    LOGICAL, INTENT(IN) :: linit_tracer_fv  !< tracer finite volume initialization
    INTEGER, INTENT(IN) :: tracer_inidist_list(:) !< selected initial tracer distributions

    INTEGER  :: ikp1, nblks_e, nblks_c, nblks_v, npromz_e, npromz_c, npromz_v, &
                nlen, je, jk, jc, jt, jb, jv, ist, it4, it5, it6, it7, it8
    INTEGER  :: nlev, nlevp1              !< number of full and half levels

    REAL(wp) :: zlon, zlat, zpk, zpkp1, zpres, zheight ! location
    REAL(wp) :: zu, zv, zq4, zq5, zq6, zq7, zq8    ! initialized variables
    REAL(wp) :: z_aleph

    REAL(wp), ALLOCATABLE :: zhelp_c (:,:,:)

!--------------------------------------------------------------------
!
    CALL init_ncar_testcases_domain()

    z_aleph = p_rotate_axis_deg * pi/180.0_wp

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    ! allocate memory for temporary field
    nblks_c = ptr_patch%nblks_c
    ALLOCATE( zhelp_c(nproma,nlev,nblks_c), STAT=ist)
    IF(ist/=SUCCESS)THEN
       CALL finish('mo_pa_test:init_hydro_state_prog_patest',  &
         &         'allocation of ptr_c failed')
    ENDIF

    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c
    nblks_e  = ptr_patch%nblks_e
    npromz_e = ptr_patch%npromz_e
    nblks_v  = ptr_patch%nblks_v
    npromz_v = ptr_patch%npromz_v

    ! topography
    ptr_ext_data%atm%topography_c(:,:) = 0.0_wp

    ! init surface pressure
    ptr_prog%pres_sfc(:,:) = zp0

    ! init temperature (isothermal atmosphere)
    ptr_prog%temp(:,:,:) = zt0

    ! init surface geopotential
    ptr_diag%geo_ic(:,nlevp1,:) = 0._wp

    it4 = 0
    it5 = 0
    it6 = 0
    it7 = 0
    it8 = 0
    DO jt = 1, ntracer
      SELECT CASE(tracer_inidist_list(jt))
      CASE(4)
        it4 = jt
      CASE(5)
        it5 = jt
      CASE(6)
        it6 = jt
      CASE(7)
        it7 = jt
      CASE(8)
        it8 = jt
      END SELECT
    ENDDO
    IF ( it4==0 .AND. it5==0 .AND. it6==0 .AND. it7==0 .AND. it8==0) THEN
      CALL finish('mo_pa_test:init_hydro_state_prog_patest',  &
         &         'incorrect tracer list')
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen,ikp1)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      ! compute half-level pressure
      CALL half_level_pressure( ptr_prog%pres_sfc(:,jb), nproma, nlen, &! input
                                ptr_diag%pres_ic(:,:,jb)              ) ! output
      ! compute pressure value at full levels
      CALL full_level_pressure( ptr_diag%pres_ic(:,:,jb),nproma,nlen, &! in
                                ptr_diag%pres_mc(:,:,jb)             ) ! out

      DO jk = 1, nlev
        ikp1 = jk+1
        ! compute layer thickness at cell centers
        zhelp_c(1:nlen,jk,jb) =  ptr_diag%pres_ic(1:nlen,ikp1,jb)&
                               - ptr_diag%pres_ic(1:nlen,jk,jb)
        ptr_diag%delp_c(1:nlen,jk,jb)  = zhelp_c(1:nlen,jk,jb)
        ptr_diag%rdelp_c(1:nlen,jk,jb) = 1.0_wp/zhelp_c(1:nlen,jk,jb)
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! compute layer thickness at edge centers
    CALL cells2edges_scalar( zhelp_c, ptr_patch, ptr_int%c_lin_e, &
                             ptr_diag%delp_e )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,nlen,zlon,zlat,zu,zv)

    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
         nlen = nproma
      ELSE
         nlen = npromz_e
      ENDIF
      DO jk = 1, nlev

        DO je = 1, nlen

          ! location of edge midpoint
          zlon = ptr_patch%edges%center(je,jb)%lon
          zlat = ptr_patch%edges%center(je,jb)%lat

          CALL init_pure_adv_wind ( zlon, zlat, p_rotate_axis_deg, zu, zv )

          ! calculate normal wind component
          ptr_prog%vn(je,jk,jb) = &
                      zu * ptr_patch%edges%primal_normal(je,jb)%v1   &
                    + zv * ptr_patch%edges%primal_normal(je,jb)%v2

          ! Coriolis parameter
          ptr_patch%edges%f_e(je,jb) = 2.0_wp*grid_angular_velocity * &
             & (SIN(zlat)*COS(z_aleph)-COS(zlon)*COS(zlat)*SIN(z_aleph))

        ENDDO ! edge loop

      ENDDO

    ENDDO
!$OMP END DO


!$OMP DO PRIVATE(jb,jk,jc,jt,nlen,ikp1,zpk,zpkp1,zpres,zheight, &
!$OMP            zlon,zlat,zq4,zq5,zq6,zq7,zq8)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jk = 1, nlev

        ! calculate pressure at main levels and height
        ikp1 = jk+1
        zpk     = vct_a(jk)   + vct_b(jk)   * zp0
        zpkp1   = vct_a(ikp1) + vct_b(ikp1) * zp0
        IF(zpk/=0.0_wp) THEN
          zpres   = EXP( (zpkp1*LOG(zpkp1)-zpk*LOG(zpk))/(zpkp1-zpk) - 1._wp)
        ELSE
          zpres   = zpkp1*0.5_wp
        ENDIF
        zheight = rd * zt0 * rgrav * LOG( zp0 / zpres )

        DO jc = 1, nlen

          ! location of quadrature points
          IF (linit_tracer_fv) THEN
            zlon = ptr_int%gquad%qpts_tri_l(jc,jb)%lon
            zlat = ptr_int%gquad%qpts_tri_l(jc,jb)%lat
          ELSE
            zlon = ptr_patch%cells%center(jc,jb)%lon
            zlat = ptr_patch%cells%center(jc,jb)%lat
          ENDIF


          CALL init_pure_adv_tracers ( tracer_inidist_list(1:5), zlon, zlat, zheight, &
            &                          p_rotate_axis_deg, zq4, zq5, zq6,              &
            &                          zq7, zq8 )

          IF (it4 /= 0) ptr_prog%tracer(jc,jk,jb,it4) = zq4
          IF (it5 /= 0) ptr_prog%tracer(jc,jk,jb,it5) = zq5
          IF (it6 /= 0) ptr_prog%tracer(jc,jk,jb,it6) = zq6
          IF (it7 /= 0) ptr_prog%tracer(jc,jk,jb,it7) = zq7
          IF (it8 /= 0) ptr_prog%tracer(jc,jk,jb,it8) = zq8

          ! location of cell circumcenter
          zlon = ptr_patch%cells%center(jc,jb)%lon
          zlat = ptr_patch%cells%center(jc,jb)%lat

          ! Coriolis parameter
          ptr_patch%cells%f_c(jc,jb) = 2.0_wp*grid_angular_velocity * &
            & (SIN(zlat)*COS(z_aleph)-COS(zlon)*COS(zlat)*SIN(z_aleph))
        ENDDO ! cell loop

      ENDDO ! vertical level loop

    ENDDO ! block loop
!$OMP END DO

!$OMP DO PRIVATE(jb,jv,nlen,zlat,zlon)
    DO jb = 1, nblks_v
      IF (jb /= nblks_v) THEN
         nlen = nproma
      ELSE
         nlen = npromz_v
      ENDIF
      DO jv = 1, nlen
        zlat   = ptr_patch%verts%vertex(jv,jb)%lat
        zlon   = ptr_patch%verts%vertex(jv,jb)%lon
        ! Coriolis parameter
        ptr_patch%verts%f_v(jv,jb) = 2.0_wp*grid_angular_velocity * &
          & (SIN(zlat)*COS(z_aleph)-COS(zlon)*COS(zlat)*SIN(z_aleph))
      ENDDO

    ENDDO
!$OMP END PARALLEL


    ! deallocate memory for temporary field
    DEALLOCATE( zhelp_c, STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish('mo_pa_test:init_hydro_state_prog_patest',  &
        &         'deallocation of zhelp_c failed')
    ENDIF

  END SUBROUTINE init_hydro_state_prog_patest

!-------------------------------------------------------------------------
!
  !>
  !!               Set the time-variant vertical velocity for.
  !!
  !!               Set the time-variant vertical velocity for
  !!               pure 3d advection test.
  !!
  !! @par Revision History
  !!  Original version by Jochen Foerstner, DWD (2008-05)
  !!  Code restructuring by Almut Gassmann, MPI-M (2008-10)
  !!
  SUBROUTINE set_vertical_velocity( ptr_patch, ptr_diag, k_step, p_sim_time )


    TYPE(t_patch), INTENT(IN)         :: ptr_patch
    TYPE(t_hydro_atm_diag), INTENT(INOUT) :: ptr_diag

    INTEGER,  INTENT(IN) :: k_step      ! actual time step
    REAL(wp), INTENT(IN) :: p_sim_time  ! simulation time in seconds since start

    INTEGER  :: jb, jk, jc          ! loop indices
    INTEGER  :: ikm1, nlen, nblks_c, npromz_c, ist, ikcenter
    INTEGER  :: nlev, nlevp1            !< number of full and half levels

    REAL(wp) :: zetah                   ! location
    REAL(wp) :: zetah_top               ! top level coordinate
    REAL(wp) :: zshape                  ! shape function
    REAL(wp) :: zdpdeta                 ! dp / deta
    REAL(wp) :: zwpres, zweta, zweta0   ! vertical velocity
    REAL(wp) :: zeta, zetam1            ! full level vertical coordinates
    REAL(wp) :: zp, zpm1                ! full level pressures

    ! for writing data
    INTEGER, SAVE :: n_file_ti          ! file identifier
    CHARACTER (LEN=MAX_CHAR_LENGTH) :: file_ti
    REAL(wp) :: zsim_days, zeta_dot, zeta_dot_star

!--------------------------------------------------------------------
!
    zdpdeta = 0.0_wp
    zwpres  = 0.0_wp

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    ikcenter = nlev/2

    ! Open the datafile in first time step
    IF ( my_process_is_stdio() .AND. k_step == 1 ) THEN
      file_ti   = 'vertical_velocity.dat'
      n_file_ti = find_next_free_unit(10,20)
      OPEN( UNIT=n_file_ti, FILE=TRIM(file_ti), FORM='FORMATTED', IOSTAT=ist )
      IF (ist/=SUCCESS) THEN
         CALL finish('mo_pa_test:set_vertical_velocity',&
                     'could not open datafile')
      ENDIF
      WRITE (n_file_ti,'(4A22,A8,I2,A4,F9.5)') &
        &  ' TIMESTEP            ,',      &
        &  ' ELAPSED TIME  (days),',      &
        &  ' ETA_DOT        (1/s),',      &
        &  ' grid_angular_velocity         (Pa/s) ',      &
        &  ' at ETA(', ikcenter, ') = ',  &
        &   cetah(ikcenter)
    ENDIF

    ! compute vertical velocity time dependent part
    zweta0 = zomega0 / zp0 * COS( 2._wp*pi/ztau * p_sim_time )


    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c

    ! top vertical half level coordinate
    zetah_top = cetah(1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen,ikm1,zetah,zetam1,zeta,zshape, &
!$OMP            zpm1,zp,zwpres,zdpdeta,zweta,zeta_dot, &
!$OMP            zeta_dot_star)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      ! top and bottom boundary values
      ptr_diag%weta(:,1,jb)        = 0._wp
      ptr_diag%weta(:,nlevp1,jb)   = 0._wp
      ptr_diag%wpres_mc(:,1,jb)    = 0._wp
      ptr_diag%wpres_mc(:,nlev,jb) = 0._wp

      DO jk = 2, nlev

        ikm1 = jk-1

        ! vertical coordinate at half level
        zetah  = cetah(jk)
        ! vertical coordinate at top and bottom main level
        zetam1 = ceta(ikm1)
        zeta   = ceta(jk)

        ! shape function
        zshape = MIN( 1._wp,   &
                 2._wp*SQRT( SIN( (zetah-zetah_top)/(1._wp-zetah_top)*pi) ) )


        DO jc=1, nlen

          ! pressure at top and bottom main level
          zpm1 = ptr_diag%pres_mc(jc,ikm1,jb)
          zp   = ptr_diag%pres_mc(jc,jk  ,jb)

          ! compute vertical velocity
          zwpres  = zweta0 * SIN( zshape * pi_2 )
          zdpdeta = ( zp - zpm1 ) / ( zeta - zetam1 )
          zweta   = zwpres * zdpdeta

          ! set vertical velocity value
          ptr_diag%wpres_mc(jc,jk,jb) = zwpres
          ptr_diag%weta(jc,jk,jb)     = zweta

        ENDDO ! cell loop

        IF ( jk == ikcenter ) THEN
          zeta_dot      = zwpres
          zeta_dot_star = zwpres * zdpdeta
        ENDIF

      ENDDO ! vertical level loop

    ENDDO ! block loop
!$OMP END DO
!$OMP END PARALLEL

    IF(my_process_is_stdio()) THEN
      zsim_days = p_sim_time / rdaylen
      ! write vertical velocity to data file
      WRITE(n_file_ti,'(i21,f22.8,e22.8,e22.8)') &
            k_step, zsim_days, zeta_dot, zeta_dot_star

      IF ( k_step == nsteps ) THEN
        CLOSE(n_file_ti)
      ENDIF
    ENDIF

  END SUBROUTINE set_vertical_velocity

!--------------------------------------------------------------------
  END MODULE mo_pa_test
