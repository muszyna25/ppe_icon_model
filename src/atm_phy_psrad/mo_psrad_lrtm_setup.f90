
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_lrtm_setup

  USE mo_psrad_general, ONLY : wp, ngpt_orig, nbndlw
  USE mo_psrad_reduce
  USE mo_psrad_lrtm_netcdf, ONLY : lrtm_read
  USE mo_psrad_lrtm_kgs


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setup_lrtm

  TYPE(ptr3), DIMENSION(2,nbndlw) :: kmajor3_in
  TYPE(ptr4), DIMENSION(2,nbndlw) :: kmajor4_in
  TYPE(ptr2), DIMENSION(2,nbndlw) :: h2oref_in
  TYPE(ptr3), DIMENSION(max_minor_species,2,nbndlw) :: kgas_in
  TYPE(ptr2), DIMENSION(2,nbndlw) :: planck_fraction2_in

  REAL(wp) :: pa1, pa2, pa3, pb1, pb2, pc1, pc2, pc3
  PUBLIC :: pa1, pa2, pa3, pb1, pb2, pc1, pc2, pc3

CONTAINS

  SUBROUTINE allocate_all
    INTEGER :: band, atm, which, i, j, gas, M, N
    DO band = 1,nbndlw
      DO atm = 1,2
        nullify(planck_fraction2(atm,band)%v, &
          planck_fraction2_in(atm,band)%v)
        nullify(kmajor(atm,band)%v, kmajor3_in(atm,band)%v, &
          kmajor4_in(atm,band)%v)
        DO i= 1, max_minor_species
          nullify(kgas(i,atm,band)%v)
          nullify(kgas_in(i,atm,band)%v)
        ENDDO
      ENDDO
    ENDDO
    DO band = 1,nbndlw
      DO which = 1,2
        M = nh2oref(which,band)
        ALLOCATE(h2oref(which,band)%v(M,ngpt(band)))
        ALLOCATE(h2oref_in(which,band)%v(M,ngpt_orig))
      END DO
      DO atm = 1,2
        n = nsp_fraction(atm,band)
        allocate(planck_fraction2(atm,band)%v(n,ngpt(band)), &
          planck_fraction2_in(atm,band)%v(n,ngpt_orig))
        n = nsp(atm,band)
        IF (n /= 0) THEN
          IF (n == 1) THEN
            allocate(kmajor3_in(atm,band)%v(5,npressure(atm),ngpt_orig))
          ELSE
            allocate(kmajor4_in(atm,band)%v(n,5,npressure(atm),ngpt_orig))
          ENDIF
          allocate(kmajor(atm,band)%v(5 * n * npressure(atm),ngpt(band)))
        ENDIF
        DO i = 1, max_minor_species
          gas = minor_species(i,atm,band)%gas
          IF (gas == 0) EXIT
          IF (gas <= ngas) THEN
            M = minor_species(i,atm,band)%M
            N = minor_species(i,atm,band)%N
            IF (N == 0) THEN
              allocate(kgas(i,atm,band)%v(M,ngpt(band)))
              allocate(kgas_in(i,atm,band)%v(M,ngpt_orig,1))
            ELSE
              allocate(kgas(i,atm,band)%v(M*N,ngpt(band)))
              allocate(kgas_in(i,atm,band)%v(M,N,ngpt_orig))
            ENDIF
          ELSE
            IF (atm == 1) THEN
              allocate(kgas(i,1,band)%v(1,ngpt(band)), &
                kgas_in(i,1,band)%v(ngpt_orig,1,1))
            ELSE
              DO j = 1,max_minor_species
                IF (minor_species(j,1,band)%gas == gas) THEN
                  kgas(i,2,band)%v => kgas(j,1,band)%v
                  EXIT
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE allocate_all

  SUBROUTINE deallocate_input
    INTEGER :: band, atm, which, i
    DO band = 1,nbndlw
      DO atm = 1,2
        IF (associated(kmajor3_in(atm,band)%v)) THEN
          deallocate(kmajor3_in(atm,band)%v)
        ENDIF
        IF (associated(kmajor4_in(atm,band)%v)) THEN
          deallocate(kmajor4_in(atm,band)%v)
        ENDIF
        IF (associated(planck_fraction2_in(atm,band)%v)) THEN
          deallocate(planck_fraction2_in(atm,band)%v)
        ENDIF
        DO i= 1, max_minor_species
          IF (associated(kgas_in(i,atm,band)%v)) THEN
            deallocate(kgas_in(i,atm,band)%v)
          ENDIF
        ENDDO
      ENDDO
      DO which = 1,2
        IF (associated(h2oref_in(which,band)%v)) THEN
          deallocate(h2oref_in(which,band)%v)
        ENDIF
      END DO
    ENDDO
  END SUBROUTINE deallocate_input

  SUBROUTINE reduce_all
    INTEGER :: band, i, gas, n, atm, which
    CALL init_rwgt(1)
    DO band = 1,nbndlw
      DO which = 1,2 ! self/foreign
        n = nh2oref(which,band)
        CALL reduce(band, h2oref_in(which,band)%v, &
          .true., h2oref(which,band)%v)
      ENDDO
      DO atm = 1,2 ! lower/upper
        n = nsp(atm,band)
        IF (n /= 0) THEN
          IF (n == 1) THEN
            CALL reduce_pack(band, &
              kmajor3_in(atm,band)%v, kmajor(atm,band)%v)
          ELSE
            CALL reduce_pack(band, &
              kmajor4_in(atm,band)%v, kmajor(atm,band)%v)
          ENDIF
        ENDIF
        n = nsp_fraction(atm,band)
        IF (n /= 0 .and. &
          .not. copy_planck_fraction_from_other_band(band) == atm) THEN
          CALL reduce(band, &
            planck_fraction2_in(atm,band)%v, &
            .false., planck_fraction2(atm,band)%v)
        ENDIF
        DO i= 1, max_minor_species
          gas = minor_species(i,atm,band)%gas
          IF (gas == 0) EXIT
          IF (gas <= ngas) THEN
            n = minor_species(i,atm,band)%N
            IF (n == 0) THEN
              CALL reduce(band, kgas_in(i,atm,band)%v(:,:,1), &
                .true., kgas(i,atm,band)%v)
            ELSE
              CALL reduce_pack(band, kgas_in(i,atm,band)%v, &
                kgas(i,atm,band)%v)
            ENDIF
          ELSEIF (atm == 1) THEN
            CALL reduce(band, kgas_in(i,1,band)%v(:,1,1), &
              .true., kgas(i,1,band)%v(1,:))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE reduce_all

  ! Calculate reference ratio to be used in calculation of Planck
  ! fraction in lower/upper atmosphere.
  SUBROUTINE fill_tables
    INTEGER :: i, j, k
    DO k = 1,59
    DO j = 1,ngas
    DO i = 1,ngas
    !from taumol03:
    !refrat_planck_a = chi_mls(ih2o,9)/chi_mls(ico2,9) !  P = 212.725 mb
    !refrat_planck_b = chi_mls(ih2o,13)/chi_mls(ico2,13) !  P = 95.58 mb
    !refrat_m_a = chi_mls(ih2o,3)/chi_mls(ico2,3) !  P = 706.270mb
    !refrat_m_b = chi_mls(ih2o,13)/chi_mls(ico2,13) !  P = 95.58 mb 
    !from taumol04:
    !refrat_planck_a = chi_mls(ih2o,11)/chi_mls(ico2,11) ! P =   142.5940 mb
    !refrat_planck_b = chi_mls(io3,13)/chi_mls(ico2,13) ! P = 95.58350 mb
    !from taumol05:
    !refrat_planck_a = chi_mls(ih2o,5)/chi_mls(ico2,5) ! P = 473.420 mb
    !refrat_planck_b = chi_mls(io3,43)/chi_mls(ico2,43) ! P = 0.2369 mb
    !refrat_m_a = chi_mls(ih2o,7)/chi_mls(ico2,7) ! P = 317.3480
    !from taumol07:
    !refrat_planck_a = chi_mls(ih2o,3)/chi_mls(io3,3) ! P = 706.2620 mb
    !refrat_m_a = chi_mls(ih2o,3)/chi_mls(io3,3) ! P = 706.2720 mb
    !from taumol09:
    !refrat_planck_a = chi_mls(ih2o,9)/chi_mls(ich4,9) ! P = 212 mb
    !refrat_m_a = chi_mls(ih2o,3)/chi_mls(ich4,3) ! P = 706.272 mb 
    !from taumol12:
    !refrat_planck_a = chi_mls(ih2o,10)/chi_mls(ico2,10) ! P =   174.164 mb 
    !from taumol13:
    !refrat_planck_a = chi_mls(ih2o,5)/chi_mls(in2o,5) ! P = 473.420 mb (Level 5)
    !refrat_m_a = chi_mls(ih2o,1)/chi_mls(in2o,1) ! P = 1053. (Level 1)
    !refrat_m_a3 = chi_mls(ih2o,3)/chi_mls(in2o,3) ! P = 706. (Level 3)
    !from taumol15:
    !refrat_planck_a = chi_mls(in2o,1)/chi_mls(ico2,1) ! P = 1053. mb (Level 1)
    !refrat_m_a = chi_mls(in2o,1)/chi_mls(ico2,1) ! P = 1053.
    !from taumol16:
    !refrat_planck_a = chi_mls(ih2o,6)/chi_mls(ich4,6) ! P = 387. mb (Level 6)

      IF (chi_mls(j,k) /= 0) THEN
        planck_ratio(i,j,k) = chi_mls(i,k) / chi_mls(j,k) 
      ELSE
        planck_ratio(i,j,k) = 0
      ENDIF
    END DO
    END DO
    END DO
  END SUBROUTINE fill_tables

  SUBROUTINE fill_in_missing
    INTEGER :: band, present, absent
    DO band = 1,nbndlw
      absent = copy_planck_fraction_from_other_band(band)
      IF (absent /= 0) THEN
        present = 3 - absent
        planck_fraction2(absent,band)%v = planck_fraction2(present,band)%v
      END IF
    END DO
  END SUBROUTINE fill_in_missing

  SUBROUTINE setup_tau_info
    INTEGER :: band, atm, i, t
    n_major_species = 0
    n_minor_species = 0
    stratosphere_fudge_flag = 0
    DO band = 1, nbndlw
      DO atm = 1, 2
        n_major_species(atm,band) = COUNT(major_species(1:2,atm,band) /= 0)
        t = 0
        DO i = 1,max_minor_species
          IF (minor_species(i,atm,band)%gas /= 0) t = i
        END DO
        n_minor_species(atm,band) = t
      END DO
      IF (ANY(stratosphere_fudge(:,band) /= 0)) THEN
        stratosphere_fudge_flag(band) = 1
      ENDIF
    END DO
  END SUBROUTINE setup_tau_info

  SUBROUTINE setup_lrtm(pressure_factor)
    REAL(wp), INTENT(IN) :: pressure_factor
    INTEGER :: band
    ! NOTE: hacked out of lrtm_solver...
    ! This secant and weight corresponds to the standard diffusivity 
    ! angle.  The angle is redefined for some bands.
    REAL(wp), PARAMETER :: wtdiff = 0.5_wp !, rec_6 = 0.166667_wp
    REAL(wp) :: fudge
    fudge = 2.0e+04_wp * pi

    CALL setup_tau_info
    CALL allocate_all
    CALL lrtm_read(kmajor3_in, kmajor4_in, h2oref_in, kgas_in, &
      planck_fraction2_in)
    CALL reduce_all
    CALL deallocate_input
    CALL fill_in_missing
    CALL fill_tables

   DO band = 1,nbndlw
      totplanck(:,band) = totplanck(:,band) * delwave(band) * wtdiff * fudge
   ENDDO
   pwvcm_factor = (10.0 * amw) / (amd * grav) / pressure_factor

   pa1 = 0.15_wp / pressure_factor
   pa2 = 250._wp * pressure_factor
   pa3 = 154.4_wp

   pb1 = 0.15_wp / pressure_factor
   pb2 = 95.6_wp

   pc1 = .05_wp / pressure_factor
   pc2 = 100._wp * pressure_factor
   pc3 = 900._wp

  END SUBROUTINE setup_lrtm

END MODULE mo_psrad_lrtm_setup
