!>
!! Diagnostic subroutines for the hydrostatic core.
!!
!! Currently, only mass
!! and total energy are supervised.
!!
!! @par Revision History
!! Initial release by Almut Gassmann (2008-04-25)
!!
!! @par Copyright
!! 2002-2008 by DWD and MPI-M
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
MODULE mo_ha_diagnostics

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: SUCCESS, MAX_CHAR_LENGTH,inwp,iecham,ildf_echam
  USE mo_model_domain_import,ONLY: n_dom
  USE mo_model_domain,       ONLY: t_patch
  USE mo_ext_data,           ONLY: ext_data
  USE mo_parallel_config,    ONLY: nproma
  USE mo_run_config,         ONLY: dtime, nsteps, nlev, ntracer,iforcing
  USE mo_io_config,          ONLY: no_output
  USE mo_dynamics_config,    ONLY: lshallow_water
  USE mo_physical_constants, ONLY: rgrav, cpd, grav
  USE mo_vertical_coord_table, ONLY: dela, delb
  USE mo_icoham_dyn_types,   ONLY: t_hydro_atm, t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_io_units,           ONLY: find_next_free_unit
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_sync,               ONLY: global_sum_array

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  ! ASCII files that contain the global integrals

  CHARACTER (len=MAX_CHAR_LENGTH) :: file_ti, file_tti  ! file names
  INTEGER :: n_file_ti, n_file_tti                      ! I/O units

  ! Total mass, energy and tracer at time step n-1. Used in the calculation of
  ! the relative errors respect to the previous time step.

  REAL(wp) :: total_mass_old, total_energy_old
  REAL(wp), ALLOCATABLE :: total_tracer_old(:)

  ! Total mass, energy and tracer at time step 1. Used in the calculation of
  ! the relative errors respect to time step 1.

  REAL(wp) :: total_mass_ini, total_energy_ini
  REAL(wp), ALLOCATABLE :: total_tracer_ini(:)

  PUBLIC :: init_total_integrals
  PUBLIC :: supervise_total_integrals

  CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Initialization subroutine for global integrals:
  !! opens ascii output file; allocates memory.
  !!
  !! @Revision History
  !! Separated from subroutine supervise_total_integrals by Hui Wan (MPI-M, 2011-05-24)
  !! when implementing the restart functionality.
  !!
  SUBROUTINE init_total_integrals

    INTEGER :: ist
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: thisroutine = "init_total_integrals"
    
    IF ( no_output ) RETURN

    IF (.NOT.lshallow_water) THEN
      !---------------------------------------------------------
      ! Open file and allocate memory for total mass and energy
      
      file_ti   = 'total_integrals.dat'
      n_file_ti = find_next_free_unit(10,20)
      
      OPEN(UNIT=n_file_ti,FILE=TRIM(file_ti),FORM='FORMATTED',IOSTAT=ist)
      IF (ist/=SUCCESS) THEN
        CALL finish( TRIM(thisroutine),'could not create data file' )
      ENDIF
      
      WRITE (n_file_ti,'(A22,A22,6A40)') &
           ' TIMESTEP            ,',&
           ' ELAPSED TIME    (hr),',&
           ' TOTAL MASS      (kg),',&
           ' RELATIVE ERROR to step N-1 (MASS)',&
           ' RELATIVE ERROR to step 1 (MASS)',&
           ' TOTAL ENERGY     (J),',&
           ' RELATIVE ERROR to step N-1 (ENERGY)',&
           ' RELATIVE ERROR to step 1 (ENERGY)'

      !---------------------------------------------------------
      ! Open file and allocate memory for tracer diagnostics
      
      IF (ntracer > 0) THEN

         ALLOCATE(total_tracer_old(ntracer), STAT=ist)
         IF(ist/=SUCCESS)THEN
           CALL finish ( TRIM(thisroutine), 'allocation of total_tracer_old failed')
         ENDIF
         ALLOCATE(total_tracer_ini(ntracer), STAT=ist)
         IF(ist/=SUCCESS)THEN
           CALL finish ( TRIM(thisroutine), 'allocation of total_tracer_ini failed')
         ENDIF

         file_tti   = 'tracer_total_integrals.dat'
         n_file_tti = find_next_free_unit(10,20)

         OPEN(UNIT=n_file_tti,FILE=TRIM(file_tti),FORM='FORMATTED',IOSTAT=ist)
         IF (ist/=SUCCESS) THEN
           CALL finish( TRIM(thisroutine),'could not open data file for tracers')
         ENDIF
         WRITE (n_file_tti,'(A22,A22,A22,3A40)') &
              ' TIMESTEP            ,',&
              ' ELAPSED TIME    (hr),',&
              ' TRACER NR        (#),',&
              ' TOTAL TRACER   (kg),',&
              ' RELATIVE ERROR to step N-1(TRACER)',&
              ' RELATIVE ERROR to step 1 (TRACER)'

      ENDIF ! ntracer > 0

   !=====================
   ! Shallow water model
   !=====================   
   ELSE

      file_ti   = 'total_integrals.dat'
      n_file_ti = find_next_free_unit(10,20)
      OPEN(UNIT=n_file_ti,FILE=TRIM(file_ti),FORM='FORMATTED',IOSTAT=ist)
      IF (ist/=SUCCESS) THEN
        CALL finish( TRIM(thisroutine),'could not open datafile')
      ENDIF
      WRITE (n_file_ti,'(A24,A24,A24,A24,A24,A24)') &
           ' TIMESTEP              ,',&
           ' ELAPSED TIME      (hr),',&
           ' TOTAL MASS       (m^3),',&
           ' TOTAL ENERGY (m^5/s^2),',&
           ' TOT. CIRCULAT. (m^2/s),',&
           ' TOT. ENSTROPHY (m/s^2),'

    ENDIF ! .NOT.lshallow_water
        
  END SUBROUTINE init_total_integrals
  !-------------
  !>
  !! This routine computes and prints the total mass and energy.
  !!
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2008-04-25)
  !! Modification by Hui Wan, MPI-M (2009-03-23):
  !! - calculation of the relative error of total mass and energy
  !! Modification by Daniel Reinert, DWD (2009-06-16):
  !! - calculation of the relative error of tracer concentration
  !! Modification by Pilar Ripodas, DWD (2009-07-1):
  !! - added calculation of the relative error of total mass and energy
  !!   respect to the first time step
  !! Modification by Pilar Ripodas, DWD (2010-14-1):
  !! - added calculation of the relative error of total tracer
  !!   respect to the first time step. Also a rgrav factor was missing
  !!   was missing in the computation of total tracer
  !! - substructed pressure at the top level in the calculation
  !!   of total mass (it is not always zero)
  !! Modification by Almut Gassmann, MPI-M (2010-06-01)
  !! - Compute enstrophy and circulation based on the dual grid philosophy.
  !!   Dual grid means rhombi for hexagon model and hexagons for triangular model.
  !! Modification by Hui Wan, MPI-M (2011-05-24)
  !! - Moved the initialization part into a separate subroutine.
  !!
  SUBROUTINE supervise_total_integrals (k_step,p_patch,p_hydro_state,ntimlev)

  INTEGER,                   INTENT(in)  :: k_step          ! actual time step
  TYPE(t_patch), TARGET,     INTENT(IN)  :: p_patch(n_dom)  ! Patch
  TYPE(t_hydro_atm), TARGET, INTENT(in)  :: p_hydro_state(n_dom) ! State
  INTEGER,                   INTENT(in)  :: ntimlev(n_dom)  ! time level

  TYPE(t_hydro_atm_prog), POINTER :: p_prog   ! prog state
  TYPE(t_hydro_atm_diag), POINTER :: p_diag   ! diag state

  INTEGER  :: jg, jb, jc, jk, jt, jx, & ! loop counters
       &      nlen,                   & ! array shapes
       &      npromz_c, nblks_c,      & ! array shapes
       &      npromz_v, nblks_v         ! array shapes
  
  REAL (wp):: z_total_mass, z_total_energy, z_elapsed_time, z_help, &
              z_total_circulation, z_total_enstrophy
  REAL (wp):: z_total_tracer(ntracer)
  REAL (wp):: z_mass(nproma,p_patch(1)%nblks_c)
  REAL (wp):: z_energy(nproma,p_patch(1)%nblks_c)
  REAL (wp):: z_aux_tracer(nproma,p_patch(1)%nblks_c,ntracer)
  
  REAL (wp), ALLOCATABLE :: z_circulation(:,:), z_enstrophy(:,:)
  REAL (wp), POINTER :: ptr_cori(:,:), ptr_area(:,:)

  ! Relative errors respect to step n-1

  REAL (wp) :: z_rel_err_mass             ! relative error of total mass
  REAL (wp) :: z_rel_err_energy           ! relative error of total energy
  REAL (wp) :: z_rel_err_tracer(ntracer)  ! relative error of total tracer

  ! Relative errors respect to step 1

  REAL (wp) :: z_rel_err_mass_s1             ! relative error of total mass
  REAL (wp) :: z_rel_err_energy_s1           ! relative error of total energy
  REAL (wp) :: z_rel_err_tracer_s1(ntracer)  ! relative error of total tracer
  REAL (wp) :: z_total_moist

  IF ( no_output ) RETURN
  !================================================
  ! Hydrostatic model
  !================================================
  IF (.NOT. lshallow_water) THEN

      z_elapsed_time = dtime*REAL(k_step,wp)/3600.0_wp
      
      !---------------------------
      ! Total air mass and energy
      !---------------------------
      
      z_total_mass   = 0.0_wp
      z_total_energy = 0.0_wp
      z_mass(:,:)    = 0.0_wp
      z_energy(:,:)  = 0.0_wp
      IF (ntracer > 0) THEN
        z_total_tracer(:)   = 0.0_wp
        z_aux_tracer(:,:,:) = 0.0_wp
      ENDIF

      jg = 1 ! It does not make sense to double-count nested domains!

      p_prog => p_hydro_state(jg)%prog(ntimlev(jg))
      p_diag => p_hydro_state(jg)%diag

      nblks_c   = p_patch(jg)%nblks_int_c
      npromz_c  = p_patch(jg)%npromz_int_c

#ifndef __SUNPRO_F95
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jc,z_help)
#endif
      DO jb = 1, nblks_c

        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF

        DO jc = 1, nlen
          z_help          = p_patch(jg)%cells%area(jc,jb)* &
                            (p_prog%pres_sfc(jc,jb)-p_diag%pres_ic(jc,1,jb))
          z_mass(jc,jb)   = z_help*rgrav
          z_energy(jc,jb) = z_help*ext_data(jg)%atm%topography_c(jc,jb)
        ENDDO
      ENDDO
#ifndef __SUNPRO_F95
!$OMP END DO
!$OMP DO PRIVATE(jb,nlen,jk,jc)
#endif
      DO jb = 1, nblks_c

        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF

        DO jk = 1, nlev
          DO jc = 1, nlen
            z_energy(jc,jb) = z_energy(jc,jb) + &
                 rgrav*p_patch(jg)%cells%area(jc,jb)*&
                 (dela(jk)+delb(jk)*p_prog%pres_sfc(jc,jb))*&
                 (p_prog%temp(jc,jk,jb)*cpd+&
                 p_diag%e_kin(jc,jk,jb))
          ENDDO
        ENDDO

      ENDDO
#ifndef __SUNPRO_F95
!$OMP END DO
#endif

      !---------------------------
      ! Total mass of tracers
      !---------------------------
      IF (ntracer >0) THEN

        DO jt=1, ntracer

#ifndef __SUNPRO_F95
!$OMP DO PRIVATE(jb,jk,jc,nlen)
#endif
          DO jb = 1, nblks_c

            IF (jb /= nblks_c) THEN
              nlen = nproma
            ELSE
              nlen = npromz_c
            ENDIF

            DO jk = 1, nlev
              DO jc = 1, nlen
                z_aux_tracer(jc,jb,jt) = z_aux_tracer(jc,jb,jt) + p_prog%tracer(jc,jk,jb,jt) * &
                  p_patch(jg)%cells%area(jc,jb)*rgrav * p_diag%delp_c(jc,jk,jb)
              ENDDO
            ENDDO

          ENDDO

#ifndef __SUNPRO_F95
!$OMP END DO
#endif

        ENDDO  ! ntracer
      ENDIF

#ifndef __SUNPRO_F95
!$OMP END PARALLEL
#endif

      WHERE(.NOT.p_patch(jg)%cells%owner_mask(:,:)) z_mass(:,:) = 0._wp
      z_total_mass = global_sum_array(z_mass)

      WHERE(.NOT.p_patch(jg)%cells%owner_mask(:,:)) z_energy(:,:) = 0._wp
      z_total_energy = global_sum_array(z_energy)

      IF (ntracer > 0) THEN
        DO jt=1, ntracer
          WHERE(.NOT.p_patch(jg)%cells%owner_mask(:,:)) z_aux_tracer(:,:,jt) = 0._wp
          z_total_tracer(jt) = global_sum_array(z_aux_tracer(:,:,jt))
        ENDDO
      ENDIF

      IF (iforcing==inwp.OR.iforcing==iecham.OR.iforcing==ildf_echam) THEN
          z_total_moist = (SUM(z_total_tracer(1:3)))
      ENDIF


     IF(my_process_is_stdio()) THEN

      ! calculate the relative error of total mass, energy and tracer
      IF (k_step == 1) THEN
        z_rel_err_mass     = 0._wp
        z_rel_err_energy   = 0._wp
        z_rel_err_mass_s1  = 0._wp
        z_rel_err_energy_s1= 0._wp
        total_mass_ini = z_total_mass
        total_energy_ini = z_total_energy
        IF (ntracer > 0) THEN
          z_rel_err_tracer(:)= 0._wp
          z_rel_err_tracer_s1(:)= 0._wp
          total_tracer_ini(:) = z_total_tracer(:)
        ENDIF

      ELSE
        z_rel_err_mass     = z_total_mass/total_mass_old - 1._wp
        z_rel_err_energy   = z_total_energy/total_energy_old - 1._wp
        z_rel_err_mass_s1  = z_total_mass/total_mass_ini - 1._wp
        z_rel_err_energy_s1= z_total_energy/total_energy_ini - 1._wp

        IF (ntracer > 0) THEN
        
          DO jt=1, ntracer
            IF(total_tracer_old(jt) == 0._wp) THEN
              z_rel_err_tracer(jt) = 0._wp
            ELSE
              z_rel_err_tracer(jt) = z_total_tracer(jt)/total_tracer_old(jt) -1._wp
            ENDIF
            IF(total_tracer_ini(jt) == 0._wp) THEN
              z_rel_err_tracer_s1(jt) = 0._wp
            ELSE
              IF (iforcing > 0) THEN
                z_rel_err_tracer_s1(jt) = z_total_moist/total_tracer_ini(jt) -1._wp
              ELSE
                z_rel_err_tracer_s1(jt) = z_total_tracer(jt)/total_tracer_ini(jt) -1._wp
              ENDIF
            ENDIF
          ENDDO

        ENDIF ! ntracer > 0
      ENDIF ! my_process_is_stdio()

      ! save the total integrals for the next step
      total_mass_old = z_total_mass
      total_energy_old = z_total_energy

      WRITE(n_file_ti,'(i21,f22.8,6e40.16)') &
           k_step,z_elapsed_time, z_total_mass,  z_rel_err_mass,  &
                                  z_rel_err_mass_s1,              &
                                  z_total_energy,z_rel_err_energy,&
                                  z_rel_err_energy_s1


      IF (ntracer > 0) THEN
      ! save the total integrals for the next step
      total_tracer_old(:) = z_total_tracer(:)

        DO jt=1, ntracer
          WRITE(n_file_tti,'(i21,f22.8,i21,3e40.16)') &
                  k_step, z_elapsed_time, jt, z_total_tracer(jt), &
                  z_rel_err_tracer(jt), z_rel_err_tracer_s1(jt)
        ENDDO
      ENDIF


      IF (k_step == nsteps) THEN

        CLOSE(n_file_ti)
        IF (ntracer > 0) THEN
          CLOSE(n_file_tti)
          DEALLOCATE(total_tracer_old)
        ENDIF
      ENDIF
 
     ENDIF ! my_process_is_stdio()

  !=====================================================
  ! Shallow water model model (lshallow_water = .TRUE.)
  !=====================================================
  ELSE

      z_elapsed_time = dtime*REAL(k_step,wp)/3600.0_wp
      z_total_mass   = 0.0_wp
      z_total_energy = 0.0_wp
      z_total_circulation = 0.0_wp
      z_total_enstrophy = 0.0_wp

      jg = 1 ! It does not make sense to double-count nested domains!
      p_prog => p_hydro_state(jg)%prog(ntimlev(jg))
      p_diag => p_hydro_state(jg)%diag

      nblks_c   = p_patch(jg)%nblks_int_c
      npromz_c  = p_patch(jg)%npromz_int_c
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        DO jc = 1, nlen
          z_help          = p_patch(jg)%cells%area(jc,jb)*&
                            p_prog%pres_sfc(jc,jb)
          z_mass(jc,jb)   = z_help
          z_energy(jc,jb) = z_help* &
                            (p_diag%geo_ic(jc,2,jb)+p_diag%e_kin(jc,1,jb)+&
                             0.5_wp*grav*p_prog%pres_sfc(jc,jb))
        ENDDO
      ENDDO

      WHERE(.NOT.p_patch(jg)%cells%owner_mask(:,:)) z_mass(:,:) = 0._wp
      z_total_mass = global_sum_array(z_mass)

      WHERE(.NOT.p_patch(jg)%cells%owner_mask(:,:)) z_energy(:,:) = 0._wp
      z_total_energy = global_sum_array(z_energy)

      nblks_v   = p_patch(jg)%nblks_int_v
      npromz_v  = p_patch(jg)%npromz_int_v
      ptr_area  => p_patch(jg)%verts%dual_area
      ptr_cori  => p_patch(jg)%verts%f_v
      ALLOCATE(z_circulation(nproma,nblks_v))
      ALLOCATE(z_enstrophy(nproma,nblks_v))
      DO jb = 1, nblks_v
        IF (jb /= nblks_v) THEN
          nlen = nproma
        ELSE
          nlen = npromz_v
        ENDIF
        DO jx = 1, nlen
          z_help              = ptr_area(jx,jb)*   &
                                (p_diag%rel_vort(jx,1,jb)+ptr_cori(jx,jb))
          z_circulation(jx,jb)= z_help
          z_enstrophy(jx,jb)  = 0.5_wp * z_help * &
                                (p_diag%rel_vort(jx,1,jb)+ptr_cori(jx,jb))     &
                                /p_diag%delp_v(jx,1,jb)
        ENDDO
      ENDDO

      WHERE(.NOT.p_patch(jg)%verts%owner_mask(:,:)) z_circulation(:,:) = 0._wp
      z_total_circulation = global_sum_array(z_circulation)
      WHERE(.NOT.p_patch(jg)%verts%owner_mask(:,:)) z_enstrophy(:,:) = 0._wp
      z_total_enstrophy = global_sum_array(z_enstrophy)

      DEALLOCATE(z_circulation, z_enstrophy)

      IF(my_process_is_stdio()) THEN
        WRITE(n_file_ti,'(i24,f24.10,e24.10,e24.10,e24.10,e24.10)') &
              k_step,z_elapsed_time,z_total_mass,z_total_energy,&
              z_total_circulation,z_total_enstrophy

        IF (k_step == nsteps) THEN
          CLOSE(UNIT=n_file_ti)
        ENDIF
      ENDIF ! my_process_is_stdio()

  ENDIF !(.NOT.lshallow_water)

  END SUBROUTINE supervise_total_integrals
  !-------------

!SUBROUTINE vertical_integral (kstep,p_patch,p_hydro_state,ntimlev)

!  INTEGER,                   INTENT(in) :: k_step     ! actual time step
!  TYPE(patch),               INTENT(IN) :: p_patch(n_dom) ! Patch
!  TYPE(t_hydro_atm), TARGET, INTENT(in) :: p_hydro_state(n_dom) ! State
!  INTEGER,                   INTENT(in) :: ntimlev(n_dom)  ! time level

!  TYPE(t_hydro_atm_prog), POINTER :: p_prog   ! prog state
!  TYPE(t_hydro_atm_diag), POINTER :: p_diag   ! diag state

!
!  INTEGER  :: jg, jb, jc, jk, jv, jt,    & ! loop counters
!              nlen,                      & ! array shapes
!              npromz_c, nblks_c            ! array shapes

!    p_diag%tracer_vint(:,:,:) = 0._wp

!      jg = 1 ! It does not make sense to double-count nested domains!

!      p_prog => p_hydro_state(jg)%prog(ntimlev(jg))
!      p_diag => p_hydro_state(jg)%diag

!      nblks_c   = p_patch(jg)%nblks_int_c
!      npromz_c  = p_patch(jg)%npromz_int_c

!       DO jt=1,ntracer
!        DO jb = 1, nblks_c
!          IF (jb /= nblks_c) THEN
!            nlen = nproma
!          ELSE
!            nlen = npromz_c
!          ENDIF
!        DO jk = nlev,1,-1
!          DO jc = 1, nlen
!         p_diag%tracer_vint(jc,jb,jt)= p_diag%tracer_vint(jc,jb,jt) &
!                                     +  p_prog%tracer(jc,jk,jb,jt)           &
!                                     * vct_a(jk) + vct_b(jk)p_prog%pres_sfc(jc,jb)
!          ENDDO
!        ENDDO

!        DO jc = 1, nlen
!         p_diag%tracer_vint(jc,jb,jt)= p_diag%tracer_vint(jc,jb,jt)*rgrav
   !     ENDDO
!       ENDDO
!      ENDDO
!END SUBROUTINE vertical_integral

!-----------------------------------------------------------------------
END MODULE mo_ha_diagnostics


