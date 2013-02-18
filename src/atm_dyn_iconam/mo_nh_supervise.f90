!>
!! Initializes and controls the time stepping in the nonhydrostatic model.
!!
!!
!! @par Revision History
!! Initial release by Almut Gassmann, MPI-M (2009-02-06)
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
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_supervise

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog, t_nh_diag
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: n_dom
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime, nsteps, msg_level, &
    &                               ltransport, ntracer, lforcing, iforcing
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH, inwp, min_rlcell_int, &
                                    min_rledge_int
  USE mo_physical_constants,  ONLY: cvd
  USE mo_mpi,                 ONLY: my_process_is_stdio, get_my_mpi_all_id, &
    &                               process_mpi_stdio_id
  USE mo_io_units,            ONLY: find_next_free_unit
  USE mo_sync,                ONLY: global_sum_array, global_max
 
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: supervise_total_integrals_nh, print_maxwinds
  
  ! Needed by supervise_total_integrals_nh to keep data between steps
  REAL(wp), ALLOCATABLE, SAVE :: z_total_tracer_old(:)
  REAL(wp), ALLOCATABLE, SAVE :: z_total_tracer_0(:)

  CONTAINS

  

  !-----------------------------------------------------------------------------
  !>
  !! supervise_total_integrals_nh
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-06-03)
  !! Modification by Daniel Reinert, DWD (2010-04-30):
  !! - computation of tracer mass error
  !!
  SUBROUTINE supervise_total_integrals_nh( k_step, p_patch, p_nh_state, ntimlev, ntimlev_rcf)

  INTEGER,                INTENT(in) :: k_step            ! actual time step
  TYPE(t_patch),            INTENT(IN) :: p_patch(n_dom)    ! Patch
  TYPE(t_nh_state), TARGET, INTENT(in) :: p_nh_state(n_dom) ! State
  INTEGER,                INTENT(in) :: ntimlev(n_dom)    ! time level
  INTEGER,                INTENT(in) :: ntimlev_rcf(n_dom)! rcf time level

  REAL(wp),SAVE :: z_total_mass_0, z_total_energy_0
  REAL(wp) :: z_total_mass, z_help, z_mean_surfp,   &
  & z_total_energy, z_kin_energy, z_pot_energy, z_int_energy, &
  & z_kin_energy_re, z_int_energy_re, z_pot_energy_re, z_total_mass_re, z_total_energy_re
#ifndef NOMPI
  REAL(wp) :: z0(nproma),&
  &           z1(nproma,p_patch(1)%nblks_c), &
  &           z2(nproma,p_patch(1)%nblks_c), &
  &           z3(nproma,p_patch(1)%nblks_c), &
  &           z4(nproma,p_patch(1)%nblks_c)
#endif

  REAL (wp):: z_total_tracer(ntracer)       ! total tracer mass
  REAL (wp):: z_aux_tracer(nproma,p_patch(1)%nblks_c,ntracer)
  REAL (wp):: z_rel_err_tracer_s1(ntracer) ! relative error of total tracer
  REAL (wp):: z_rel_err_tracer(ntracer)    ! relative error of total tracer

  INTEGER :: jg, jb, jk, jc, jt             ! loop indices
  INTEGER :: nlen, npromz_c, nblks_c, istat
  INTEGER :: nlev                           ! number of full levels
  INTEGER :: ist                            ! status variable
  INTEGER, SAVE :: n_file_ti = -1, n_file_tti = -1        ! file identifier
  CHARACTER (len=MAX_CHAR_LENGTH) :: file_ti, file_tti    ! file name
  TYPE(t_nh_prog), POINTER :: p_prog       ! prog state
  TYPE(t_nh_prog), POINTER :: p_prog_rcf   ! prog_rcf state
  TYPE(t_nh_diag), POINTER :: p_diag       ! diag state
  REAL (wp):: z_total_moist, z_elapsed_time
  REAL (wp), TARGET :: z_ekin(nproma,p_patch(1)%nlev,p_patch(1)%nblks_c)
  REAL(wp), DIMENSION(:,:,:),   POINTER :: ptr_ekin
  !-----------------------------------------------------------------------------

  ! Hack [ha]:
  IF (.NOT. ALLOCATED (z_total_tracer_old)) THEN
     ALLOCATE (z_total_tracer_old(ntracer), STAT=ist)
     IF(ist/=SUCCESS)THEN
        CALL finish ('mo_nh_stepping:supervise_total_integrals_nh', &
                     'allocation of z_total_tracer_old failed')
     ENDIF
     ALLOCATE (z_total_tracer_0(ntracer), STAT=ist)
     IF(ist/=SUCCESS)THEN
        CALL finish ('mo_nh_stepping:supervise_total_integrals_nh', &
                     'allocation of z_total_tracer_0 failed')
     ENDIF
     z_total_tracer_old = 0.0_wp
     z_total_tracer_0   = 0.0_wp
  END IF

  ! Open the datafile
  IF (k_step == 1 .AND. my_process_is_stdio()) THEN
    file_ti   = 'total_integrals.dat'
    n_file_ti = find_next_free_unit(10,20)
    OPEN(UNIT=n_file_ti,FILE=TRIM(file_ti),FORM='FORMATTED',IOSTAT=istat)
    IF (istat/=SUCCESS) THEN
      CALL finish('supervise_total_integrals_nh','could not open datafile')
    ENDIF
    WRITE (n_file_ti,'(A8,6A20)')'TIMESTEP',&
             '            m/m0 -1,',&
             '            e/e0 -1,',&
             '             % kine,',&
             '             % inne,',&
             '             % pote,',&
             '    mean surf press.'

    ! Open the datafile for tracer diagnostic
    IF (ltransport .OR. ( iforcing == inwp ) ) THEN

!!$       ALLOCATE(z_total_tracer_old(ntracer), STAT=ist)
!!$       IF(ist/=SUCCESS)THEN
!!$          CALL finish ('mo_nh_stepping:supervise_total_integrals_nh', &
!!$               'allocation of z_total_tracer_old failed')
!!$       ENDIF
!!$       ALLOCATE(z_total_tracer_0(ntracer), STAT=ist)
!!$       IF(ist/=SUCCESS)THEN
!!$          CALL finish ('mo_nh_stepping:supervise_total_integrals_nh', &
!!$               'allocation of z_total_tracer_0 failed')
!!$       ENDIF


       file_tti   = 'tracer_total_integrals.dat'
       n_file_tti = find_next_free_unit(10,20)
       OPEN(UNIT=n_file_tti,FILE=TRIM(file_tti),FORM='FORMATTED',IOSTAT=istat)
       IF (istat/=SUCCESS) THEN
         CALL finish('supervise_total_integrals_nh','could not open datafile')
       ENDIF
       WRITE (n_file_tti,'(A22,A22,A22,3A40)') &
            ' TIMESTEP            ,',&
            ' ELAPSED TIME    (hr),',&
            ' TRACER NR        (#),',&
            ' TOTAL TRACER   (kg),',&
            ' RELATIVE ERROR to step N-1(TRACER)',&
            ' RELATIVE ERROR to step 1 (TRACER)'
      ENDIF

  ENDIF

  z_elapsed_time = dtime*REAL(k_step,wp)/3600.0_wp

  jg = 1 ! It does not make sense to double-count nested domains!

  p_prog     => p_nh_state(jg)%prog(ntimlev(jg))
  p_prog_rcf => p_nh_state(jg)%prog(ntimlev_rcf(jg))
  p_diag     => p_nh_state(jg)%diag

  IF (p_patch(jg)%cell_type == 3) THEN
    ptr_ekin => z_ekin
  ELSE
    ptr_ekin => p_diag%e_kin
  ENDIF

  nblks_c   = p_patch(jg)%nblks_int_c
  npromz_c  = p_patch(jg)%npromz_int_c

  ! number of vertical levels
  nlev = p_patch(jg)%nlev

  DO jb = 1, nblks_c
    IF (jb /= nblks_c) THEN
      nlen = nproma
    ELSE
      nlen = npromz_c
    ENDIF
    IF (p_patch(jg)%cell_type == 3) THEN
      DO jk = 1, nlev
        DO jc = 1, nlen
          ptr_ekin(jc,jk,jb) = p_diag%e_kinh(jc,jk,jb) + 0.25_wp* &
            (p_prog%w(jc,jk,jb)**2 + p_prog%w(jc,jk+1,jb)**2)
        ENDDO
      ENDDO
    ELSE
      DO jk = 1, nlev
        ptr_ekin(1:nlen,jk,jb) = p_diag%e_kinh(1:nlen,jk,jb) +0.25_wp &
          &*(p_prog%w(1:nlen,jk  ,jb)**2*p_nh_state(jg)%metrics%ddqz_z_half(1:nlen,jk  ,jb) &
          & +p_prog%w(1:nlen,jk+1,jb)**2*p_nh_state(jg)%metrics%ddqz_z_half(1:nlen,jk+1,jb))&
          & /p_nh_state(jg)%metrics%ddqz_z_full(1:nlen,jk,jb)
      ENDDO
    ENDIF
  ENDDO

#ifdef NOMPI

  z_total_mass = 0.0_wp
  z_kin_energy = 0.0_wp
  z_pot_energy = 0.0_wp
  z_int_energy = 0.0_wp
  z_mean_surfp = 0.0_wp
  DO jb = 1, nblks_c
    IF (jb /= nblks_c) THEN
      nlen = nproma
    ELSE
      nlen = npromz_c
    ENDIF
    DO jk = 1, nlev
      DO jc = 1, nlen
        z_help = p_patch(jg)%cells%area(jc,jb)*p_nh_state(jg)%metrics%ddqz_z_full(jc,jk,jb) &
        &       /p_patch(jg)%n_patch_cells_g
        z_total_mass = z_total_mass + &
        &              p_prog%rho(jc,jk,jb)*z_help
        z_kin_energy = z_kin_energy + &
        &              p_prog%rho(jc,jk,jb)*ptr_ekin(jc,jk,jb)*z_help
        z_int_energy = z_int_energy + &
        &              cvd*p_prog%exner(jc,jk,jb)*p_prog%rhotheta_v(jc,jk,jb)*z_help
        z_pot_energy = z_pot_energy + &
        &              p_prog%rho(jc,jk,jb)*p_nh_state(jg)%metrics%geopot(jc,jk,jb)*z_help
      ENDDO
    ENDDO
    DO jc = 1, nlen
      z_mean_surfp = z_mean_surfp + p_diag%pres_sfc(jc,jb)/      &
                     &  REAL(p_patch(jg)%n_patch_cells_g,wp)
    ENDDO
  ENDDO
  z_total_energy = z_int_energy+z_kin_energy+z_pot_energy

#else

  DO jb = 1, nblks_c
    IF (jb /= nblks_c) THEN
      nlen = nproma
    ELSE
      nlen = npromz_c
    ENDIF
    z1(1:nlen,jb) = 0.0_wp
    z2(1:nlen,jb) = 0.0_wp
    z3(1:nlen,jb) = 0.0_wp
    z4(1:nlen,jb) = 0.0_wp
    DO jk = 1,nlev
      z0(1:nlen) = p_patch(jg)%cells%area(1:nlen,jb)      &
      & *p_nh_state(jg)%metrics%ddqz_z_full(1:nlen,jk,jb) &
      & /REAL(p_patch(jg)%n_patch_cells_g,wp)
      z1(1:nlen,jb) = z1(1:nlen,jb)&
      & +p_prog%rho(1:nlen,jk,jb)*z0(1:nlen)
      z2(1:nlen,jb) = z2(1:nlen,jb)&
      & +p_prog%rho(1:nlen,jk,jb)*ptr_ekin(1:nlen,jk,jb)*z0(1:nlen)
      z3(1:nlen,jb) = z3(1:nlen,jb)&
      & +cvd*p_prog%exner(1:nlen,jk,jb)*p_prog%rhotheta_v(1:nlen,jk,jb)*z0(1:nlen)
      z4(1:nlen,jb) = z4(1:nlen,jb)&
      & +p_prog%rho(1:nlen,jk,jb)*p_nh_state(jg)%metrics%geopot(1:nlen,jk,jb)*z0(1:nlen)
    ENDDO
    WHERE(.NOT.p_patch(jg)%cells%owner_mask(:,jb)) z1(:,jb) = 0._wp
    WHERE(.NOT.p_patch(jg)%cells%owner_mask(:,jb)) z2(:,jb) = 0._wp
    WHERE(.NOT.p_patch(jg)%cells%owner_mask(:,jb)) z3(:,jb) = 0._wp
    WHERE(.NOT.p_patch(jg)%cells%owner_mask(:,jb)) z4(:,jb) = 0._wp
    WHERE(.NOT.p_patch(jg)%cells%owner_mask(:,jb)) p_diag%pres_sfc(:,jb) = 0._wp
  ENDDO
  z_mean_surfp = global_sum_array( p_diag%pres_sfc )/      &
                 &     REAL(p_patch(jg)%n_patch_cells_g,wp)
  z_total_mass = global_sum_array( z1 )
  z_kin_energy = global_sum_array( z2 )
  z_int_energy = global_sum_array( z3 )
  z_pot_energy = global_sum_array( z4 )
  z_total_energy = z_int_energy+z_kin_energy+z_pot_energy

#endif

  IF (k_step == 1) THEN
    z_total_mass_0   = z_total_mass
    z_total_energy_0 = z_total_energy
  ENDIF

  ! percentage of single energies out of the total energy
  z_kin_energy_re = z_kin_energy/z_total_energy*100.0_wp
  z_int_energy_re = z_int_energy/z_total_energy*100.0_wp
  z_pot_energy_re = z_pot_energy/z_total_energy*100.0_wp
  ! changes compared to first step
  z_total_mass_re   = z_total_mass  /z_total_mass_0  -1.0_wp
  z_total_energy_re = z_total_energy/z_total_energy_0-1.0_wp

  IF (my_process_is_stdio()) THEN
    if (n_file_ti >= 0) &
    WRITE(n_file_ti,'(i8,6e20.12)') &
    &   k_step, z_total_mass_re, z_total_energy_re, z_kin_energy_re, z_int_energy_re, &
    &   z_pot_energy_re, z_mean_surfp
    IF (k_step == nsteps) THEN
      if (n_file_ti >= 0) &
      CLOSE(n_file_ti)
    ENDIF
  ENDIF

  IF (ltransport  .OR. ( iforcing == inwp ) ) THEN

    z_total_tracer(:)   = 0.0_wp
    z_aux_tracer(:,:,:) = 0.0_wp
    z_total_moist       = 0.0_wp

    DO jt=1, ntracer

      DO jb = 1, nblks_c

        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF

        ! compute tracer mass in each vertical column
        DO jk = 1, nlev
          DO jc = 1, nlen
            z_help = p_patch(jg)%cells%area(jc,jb)             &
              &    * p_nh_state(jg)%metrics%ddqz_z_full(jc,jk,jb) &
              &    * p_prog%rho(jc,jk,jb)

            z_aux_tracer(jc,jb,jt) = z_aux_tracer(jc,jb,jt)    &
              &    + p_prog_rcf%tracer(jc,jk,jb,jt) * z_help
          ENDDO
        ENDDO

      ENDDO

      WHERE(.NOT.p_patch(jg)%cells%owner_mask(:,:)) z_aux_tracer(:,:,jt) = 0._wp
      z_total_tracer(jt) = global_sum_array(z_aux_tracer(:,:,jt))

    ENDDO  ! ntracer


    IF (lforcing) THEN
       ! when run in physics mode, compute total mass of all tracers
       ! (water substances)
       z_total_moist = (SUM(z_total_tracer(:)))
    ENDIF


    ! Save total tracer mass at first time step
    IF (k_step == 1) THEN
      z_total_tracer_0(:)    = z_total_tracer(:)
      z_rel_err_tracer_s1(:) = 0._wp
      z_rel_err_tracer(:)    = 0._wp

    ELSE

      ! compute relative mass error
      ! a) relative to time step 1 (z_rel_err_tracer_s1)
      ! b) relative to the previous time step n-1 (z_rel_err_tracer)

      DO jt=1, ntracer

        IF (z_total_tracer_old(jt) == 0._wp) THEN
          z_rel_err_tracer(jt) = 0._wp
        ELSE
          z_rel_err_tracer(jt) = (z_total_tracer(jt)/z_total_tracer_old(jt)) - 1._wp
        ENDIF
        IF (z_total_tracer_0(jt) == 0._wp) THEN
          z_rel_err_tracer_s1(jt) = 0._wp
        ELSE
          IF (lforcing) THEN
            z_rel_err_tracer_s1(jt) = (z_total_moist/z_total_tracer_0(jt)) - 1._wp
          ELSE
            z_rel_err_tracer_s1(jt) = (z_total_tracer(jt)/z_total_tracer_0(jt)) - 1._wp
          ENDIF
        ENDIF

      ENDDO

    ENDIF


    ! save total tracer mass for the next step
    z_total_tracer_old(:) = z_total_tracer(:)


    DO jt=1, ntracer
      if (n_file_tti > 0) &
      WRITE(n_file_tti,'(i21,f22.8,i21,3e40.16)')             &
              k_step, z_elapsed_time, jt, z_total_tracer(jt), &
              z_rel_err_tracer(jt), z_rel_err_tracer_s1(jt)
    ENDDO

    IF (k_step == nsteps) THEN
      if (n_file_ti >= 0) &
      CLOSE(n_file_ti)
      if (n_file_tti > 0) &
      CLOSE(n_file_tti)

       DEALLOCATE(z_total_tracer_old, z_total_tracer_0, STAT=ist)
       IF(ist/=SUCCESS)THEN
          CALL finish ('mo_nh_stepping:supervise_total_integrals_nh', &
               'deallocation of z_total_tracer_old, z_total_tracer_0 failed')
       ENDIF
    ENDIF

  ENDIF    ! ltransport


  END SUBROUTINE supervise_total_integrals_nh


  !-------------------------------------------------------------------------
  !>
  !! Computation of maximum horizontal and vertical wind speed for runtime diagnostics
  !! Was included in mo_nh_stepping before
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2013-01-07)
  !!
  SUBROUTINE print_maxwinds(p_patch, vn, w)

  TYPE(t_patch), INTENT(IN) :: p_patch    ! Patch
  REAL(wp),      INTENT(IN) :: vn(:,:,:), w(:,:,:) ! horizontal and vertical wind speed

  ! local variables
  REAL(wp) :: vn_aux(p_patch%edges%end_blk(min_rledge_int,MAX(1,p_patch%n_childdom)),p_patch%nlev)
  REAL(wp) :: w_aux (p_patch%cells%end_blk(min_rlcell_int,MAX(1,p_patch%n_childdom)),p_patch%nlevp1)
  REAL(wp) :: vn_aux_lev(p_patch%nlev), w_aux_lev(p_patch%nlevp1), vmax(2)

  INTEGER  :: i_nchdom
  INTEGER  :: jb, jk, nlen, iendblk_c, iendblk_e, jg
  INTEGER  :: proc_id(2), keyval(2)


!-----------------------------------------------------------------------

    i_nchdom = MAX(1,p_patch%n_childdom)
    iendblk_c = p_patch%cells%end_blk(min_rlcell_int,i_nchdom)
    iendblk_e = p_patch%edges%end_blk(min_rledge_int,i_nchdom)
    jg        = p_patch%id

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, iendblk_e
      IF (jb /= iendblk_e) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%edges%end_idx(min_rledge_int,i_nchdom)
      ENDIF
      DO jk = 1, p_patch%nlev
        vn_aux(jb,jk) = MAXVAL(ABS(vn(1:nlen,jk,jb)))
      ENDDO
    END DO
!$OMP END DO

!$OMP DO PRIVATE(jb, jk, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, iendblk_c
      IF (jb /=  iendblk_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%cells%end_idx(min_rlcell_int,i_nchdom)
      ENDIF
      DO jk = 1, p_patch%nlevp1
        w_aux(jb,jk) = MAXVAL(ABS(w(1:nlen,jk,jb)))
      ENDDO
    END DO
!$OMP END DO

!$OMP DO PRIVATE(jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jk = 1, p_patch%nlev
      vn_aux_lev(jk) = MAXVAL(vn_aux(:,jk))
      w_aux_lev (jk) = MAXVAL(w_aux(:,jk))
    END DO
!$OMP END DO

!$OMP END PARALLEL

    ! Add surface level for w
    jk = p_patch%nlevp1
    w_aux_lev (jk) = MAXVAL(w_aux(:,jk))

    !--- Get max over all PEs
    IF (msg_level >= 8) THEN

      ! print a detailed information on global maxima:
      ! containing the process ID and the level where the
      ! maximum occurred.
      vmax(1)   = MAXVAL(vn_aux_lev(:))
      keyval(1) = MAXLOC(vn_aux_lev(:),1)
      vmax(2)   = MAXVAL(w_aux_lev(:))
      keyval(2) = MAXLOC(w_aux_lev(:),1)

      IF (msg_level >= 13) THEN
        proc_id(:) = get_my_mpi_all_id()
        vmax       = global_max(vmax, proc_id=proc_id, keyval=keyval)
        WRITE(message_text,'(a,i3,a,2(e18.10,a,i5,a,i3,a))') 'MAXABS VN, W in domain', jg, ':', &
          & vmax(1), " (on proc #", proc_id(1), ", level ", keyval(1), "), ", &
          & vmax(2), " (on proc #", proc_id(2), ", level ", keyval(2), "), "
      ELSE
        vmax = global_max(vmax, keyval=keyval)
        WRITE(message_text,'(a,i3,a,2(e18.10,a,i3,a))') 'MAXABS VN, W in domain', jg, ':', &
          & vmax(1), " at level ", keyval(1), ", ", &
          & vmax(2), " at level ", keyval(2), ", "
      END IF

    ELSE

      ! on PE0 print a short information on global maxima:
      vmax(1) = MAXVAL(vn_aux_lev)
      vmax(2) = MAXVAL(w_aux_lev)
      vmax = global_max(vmax, iroot=process_mpi_stdio_id)
      WRITE(message_text,'(a,2e18.10)') 'MAXABS VN, W ', vmax(1), vmax(2)

    END IF

    CALL message('',message_text)

  END SUBROUTINE print_maxwinds


END MODULE mo_nh_supervise


