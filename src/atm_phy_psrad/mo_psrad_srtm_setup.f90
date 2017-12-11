!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_srtm_setup

  USE mo_psrad_general, ONLY : wp, nbndsw, ngpt_orig, npressure, &
    ptr1, ptr2, ptr3, ptr4
  USE mo_psrad_reduce
  USE mo_psrad_srtm_netcdf,    ONLY : srtm_read
  USE mo_psrad_srtm_kgs, ONLY : ngpt, nh2oref, h2oref, &
    nsp, kmajor, nsfluxref, sfluxref, rayl1, rayl2, rayl_type, &
    minor_species, minor_species_missing_data, kgas

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setup_srtm

  TYPE(ptr3), DIMENSION(2,nbndsw) :: kmajor3_in
  TYPE(ptr4), DIMENSION(2,nbndsw) :: kmajor4_in
  TYPE(ptr2), DIMENSION(2,nbndsw) :: h2oref_in
  TYPE(ptr1), DIMENSION(2,nbndsw) :: kgas_in
  TYPE(ptr1) :: rayl1_in(nbndsw)
  TYPE(ptr2) :: rayl2_in(2,nbndsw)
  TYPE(ptr2), DIMENSION(nbndsw) :: sfluxref_in

CONTAINS

  SUBROUTINE setup_srtm

    CALL allocate_all
    CALL srtm_read(kmajor3_in, kmajor4_in, h2oref_in, kgas_in, &
      rayl1_in, rayl2_in, sfluxref_in)
    CALL init_rwgt(2)
    CALL reduce_all
    CALL deallocate_input

    !Data hacks from srtm_gas_optics::setup

    kgas(2,5)%v => kgas(1,5)%v

    kmajor(2,7)%v = kmajor(2,7)%v *1.6_wp
    kgas(1,7)%v = 4.35e-4_wp/(350.0_wp*2.0_wp)
    kgas(2,7)%v => kgas(1,7)%v

    kmajor(1,8)%v = kmajor(1,8)%v * 1.029_wp

    sfluxref(12)%v(1,:) = sfluxref(12)%v(1,:) * 50.15_wp/48.37_wp

  END SUBROUTINE setup_srtm

  SUBROUTINE reduce_all
    INTEGER :: band, gas, n, which
    DO band = 1,nbndsw
      IF (nsfluxref(band) /= 0) THEN
        CALL reduce(band, sfluxref_in(band)%v, .false., sfluxref(band)%v)
      ENDIF
      DO which = 1,2 ! self/foreign
        n = nh2oref(which,band);
        IF (n /= 0) THEN
          CALL reduce(band, &
            h2oref_in(which,band)%v(1:n,:), .true., &
              h2oref(which,band)%v(1:n,:))
          IF (which == 2 .and. n == 3) THEN
            h2oref(which,band)%v(4,:) = h2oref(which,band)%v(3,:)
          ENDIF
        ENDIF
      ENDDO
      DO which = 1,2 ! lower/upper
        n = nsp(which,band)
        IF (n /= 0) THEN
          IF (n == 1) THEN
            CALL reduce_pack(band, &
              kmajor3_in(which,band)%v, kmajor(which,band)%v)
          ELSE
            CALL reduce_pack(band, &
              kmajor4_in(which,band)%v, kmajor(which,band)%v)
          ENDIF
        ENDIF
        gas = minor_species(which,band)
        IF (gas /= 0 .and. &
          minor_species_missing_data(which,band) == 0 .and. &
          band /= 7) THEN
          CALL reduce(band, kgas_in(which,band)%v, .true., kgas(which,band)%v)
        ENDIF
      END DO
      DO which = 1,2 ! lower/upper
        n = rayl_type(which,band)
        IF (n /= 0) THEN
          IF (n == 1) THEN
            CALL reduce(band, rayl1_in(band)%v, .true., rayl1(band)%v)
            IF (which == 1) EXIT
          ELSE
            CALL reduce(band, &
              rayl2_in(which,band)%v, .true., rayl2(which,band)%v)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE reduce_all

  SUBROUTINE allocate_all
    INTEGER :: band, atm, which, m

    DO band = 1,nbndsw
      NULLIFY(rayl1(band)%v)
      NULLIFY(rayl1_in(band)%v)
      DO atm = 1,2
        NULLIFY(kmajor(atm,band)%v)
        NULLIFY(kmajor3_in(atm,band)%v)
        NULLIFY(kmajor4_in(atm,band)%v)
        NULLIFY(rayl2(atm,band)%v)
        NULLIFY(rayl2_in(atm,band)%v)
        NULLIFY(kgas(atm,band)%v)
        NULLIFY(kgas_in(atm,band)%v)
      END DO
      DO which = 1,2
        NULLIFY(h2oref(which,band)%v)
        NULLIFY(h2oref_in(which,band)%v)
      ENDDO
    END DO

    DO band = 1,nbndsw
      m = nsfluxref(band)
      ALLOCATE(sfluxref_in(band)%v(m,ngpt_orig))
      ALLOCATE(sfluxref(band)%v(m,ngpt(band)))

      DO which = 1,2
        m = nh2oref(which,band) ! for which == 2, may be 3 or 4
        IF (m /= 0) THEN
          ALLOCATE(h2oref_in(which,band)%v(m, ngpt_orig))
          !NOTE: Must always allocate 4 or get out-of-range.
          IF (m == 3) m = 4
          ALLOCATE(h2oref(which,band)%v(m, ngpt(band)))
        ENDIF
      ENDDO
      DO atm = 1,2
        m = rayl_type(atm,band)
        IF (m == 1 .and. .not. ASSOCIATED(rayl1(band)%v)) THEN
          ALLOCATE(rayl1(band)%v(ngpt(band)))
          ALLOCATE(rayl1_in(band)%v(ngpt_orig))
        ENDIF
        IF (m > 1) THEN
          ALLOCATE(rayl2(atm,band)%v(m,ngpt(band)))
          ALLOCATE(rayl2_in(atm,band)%v(m,ngpt_orig))
        ENDIF

        m = minor_species(atm,band)
        IF (m /= 0) THEN
          ALLOCATE(kgas(atm,band)%v(ngpt(band)))
          !FIXME: Temporary fix to compensate for hacks...
          IF (atm == 1 .or. band /= 7) THEN
            ALLOCATE(kgas_in(atm,band)%v(ngpt_orig))
          END IF
        END IF

        m = nsp(atm,band)
        IF (m /= 0) THEN
          IF (m == 1) THEN
            ALLOCATE(kmajor3_in(atm,band)%v(5,npressure(atm),ngpt_orig))
          ELSE
            ALLOCATE(kmajor4_in(atm,band)%v(m,5,npressure(atm),ngpt_orig))
          ENDIF
          ALLOCATE(kmajor(atm,band)%v(m * 5 * npressure(atm),ngpt(band)))
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE allocate_all

  SUBROUTINE deallocate_input
    INTEGER :: band, atm, which
    DO band = 1,nbndsw
      IF (ASSOCIATED(rayl1_in(band)%v)) THEN
        DEALLOCATE(rayl1_in(band)%v)
      ENDIF
      IF (ASSOCIATED(sfluxref_in(band)%v)) THEN
        DEALLOCATE(sfluxref_in(band)%v)
      ENDIF
      DO atm = 1,2
        IF (ASSOCIATED(kmajor3_in(atm,band)%v)) THEN
          DEALLOCATE(kmajor3_in(atm,band)%v)
        ENDIF
        IF (ASSOCIATED(kmajor4_in(atm,band)%v)) THEN
          DEALLOCATE(kmajor4_in(atm,band)%v)
        ENDIF
        IF (ASSOCIATED(rayl2_in(atm,band)%v)) THEN
          DEALLOCATE(rayl2_in(atm,band)%v)
        ENDIF
        IF (ASSOCIATED(kgas_in(atm,band)%v)) THEN
          DEALLOCATE(kgas_in(atm,band)%v)
        ENDIF
      ENDDO
      DO which = 1,2
        IF (associated(h2oref_in(which,band)%v)) THEN
          deallocate(h2oref_in(which,band)%v)
        ENDIF
      END DO
    ENDDO
  END SUBROUTINE deallocate_input

END MODULE mo_psrad_srtm_setup

