!>
!! Initialization routines for vertical coordinate
!!
!! Initialization routines for two height-based terrain following coordinates
!! - Hybrid Gal-Chen coordiate (Gal-Chen and Somerville, 1975; Simmons and Burridge, 1981)
!! - Smooth level vertical (SLEVE) coordinate (Schär et al, 2002; Leuenberger et al, 2010) 
!!
!! @author Guenther Zaengl, DWD
!!
!!
!! @par Revision History
!! First version by Guenther Zaengl, DWD (2011-06-29)
!! Modification by Daniel Reinert, DWD (2017-01-12)
!! - Encapsulated into separate module
!!
!!
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

MODULE mo_init_vgrid

  USE mo_kind,                  ONLY: wp
  USE mo_impl_constants,        ONLY: MAX_CHAR_LENGTH
  USE mo_math_constants,        ONLY: pi
  USE mo_exception,             ONLY: message, message_text, finish
  USE mo_run_config,            ONLY: msg_level
  USE mo_parallel_config,       ONLY: nproma
  USE mo_nonhydrostatic_config, ONLY: ivctype
  USE mo_sleve_config,          ONLY: itype_laydistr, min_lay_thckn, max_lay_thckn, htop_thcknlimit, top_height, &
                                      decay_scale_1, decay_scale_2, decay_exp, flat_height, stretch_fac
  USE mo_dynamics_config,       ONLY: iequations
  USE mo_vertical_coord_table,  ONLY: read_vct


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_hybrid_coord
  PUBLIC :: init_sleve_coord
  PUBLIC :: prepare_hybrid_coord
  PUBLIC :: prepare_sleve_coord
  PUBLIC :: init_vert_coord

  ! Dirty stuff
  ! So far used by mo_nh_testcase_nml, where these guys are set.
  REAL(wp) :: layer_thickness        ! constant layer thickness (A(k)-A(k+1)) for 
                                     ! Gal-Chen hybrid coordinate. (m)
                                     ! If layer_thickness<0,  A(k), B(k) are read 
                                     ! from file.                                    
  INTEGER  :: n_flat_level           ! Number of flat levels, i.e. where B=0.
  !
  PUBLIC :: layer_thickness, n_flat_level

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
  SUBROUTINE init_hybrid_coord(nlev, vct_a, vct_b)

!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
!!$      &  routine = 'mo_nh_init_utils:init_hybrid_coord'

    INTEGER,  INTENT(IN)    :: nlev  !< number of full levels
    REAL(wp), INTENT(INOUT) :: vct_a(:), vct_b(:)

    REAL(wp) :: z_height, z_flat
    INTEGER  :: jk
    INTEGER  :: nlevp1            !< number of half levels

    ! number of vertical half levels
    nlevp1 = nlev+1

    ! read hybrid parameters as in the hydrostatic model
    IF ( layer_thickness < 0.0_wp) THEN

      CALL read_vct (iequations,nlev)

    ELSE

      z_flat = REAL(nlevp1-n_flat_level,wp) * layer_thickness
      DO jk = 1, nlevp1
        z_height  = layer_thickness*REAL(nlevp1-jk,wp)
        vct_a(jk) = z_height
        IF ( z_height >= z_flat) THEN
          vct_b(jk) = 0.0_wp
        ELSE
          vct_b(jk) = (z_flat - z_height)/z_flat
        ENDIF
      ENDDO

    ENDIF

  END SUBROUTINE init_hybrid_coord



  !---------------------------------------------------------------------------
  !> Utility routine: computes some values based on vct_a, vct_b
  !
  SUBROUTINE prepare_hybrid_coord(nlev, vct_a, vct_b, vct, nflatlev)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_init_utils:init_hybrid_coord'

    INTEGER,  INTENT(IN)  :: nlev  !< number of full levels
    REAL(wp), INTENT(IN)  :: vct_a(:), vct_b(:)
    REAL(wp), INTENT(OUT) :: vct(:)
    INTEGER,  INTENT(OUT) :: nflatlev(:)

    REAL(wp) :: z_height, z_flat
    INTEGER  :: jk, nflat
    INTEGER  :: nlevp1            !< number of half levels

    ! number of vertical half levels
    nlevp1 = nlev+1

    IF ( layer_thickness < 0.0_wp) THEN

      CALL read_vct (iequations,nlev)

      DO jk = 1, nlevp1
        IF (vct_b(jk) /= 0.0_wp) THEN
          nflat = jk-1
          nflatlev(1) = nflat
          EXIT
        ENDIF
      ENDDO

    ELSE

      nflat  = -1
      z_flat = REAL(nlevp1-n_flat_level,wp) * layer_thickness
      DO jk = 1, nlevp1
        z_height  = layer_thickness*REAL(nlevp1-jk,wp)
        IF ( z_height < z_flat) THEN
          IF (nflat == -1) THEN
            nflat = jk-1
          ENDIF
        ENDIF
      ENDDO
      nflatlev(1) = nflat

      vct(       1:       nlevp1) = vct_a(:)
      vct(nlevp1+1:nlevp1+nlevp1) = vct_b(:)
    ENDIF

!!$    CALL message(TRIM(routine), ' coordinate setup finished')
    IF (msg_level >= 7) THEN
      CALL print_vcoord_info(nlev, vct_a, "init_hybrid_coord")
    ENDIF

  END SUBROUTINE prepare_hybrid_coord


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
  SUBROUTINE init_sleve_coord(nlev, vct_a, vct_b)

    INTEGER,  INTENT(IN)    :: nlev  !< number of full levels
    REAL(wp), INTENT(INOUT) :: vct_a(:), vct_b(:)

    REAL(wp) :: z_exp, dvct(nlev), zvcta(nlev+1), stretchfac, zdvct, x1, a, b, c, jkr
    INTEGER  :: jk, jk1, jks, jk2
    INTEGER  :: nlevp1        !< number of full and half levels

    ! number of vertical levels
    nlevp1 = nlev+1

    IF (min_lay_thckn > 0.01_wp) THEN
      IF (itype_laydistr == 1) THEN
        z_exp = LOG(min_lay_thckn/top_height)/LOG(2._wp/pi*ACOS(REAL(nlev-1,wp)**stretch_fac/&
          &     REAL(nlev,wp)**stretch_fac))

        ! Set up distribution of coordinate surfaces according to the analytical formula
        ! vct = h_top*(2/pi*arccos(jk-1/nlev))**z_exp (taken from the COSMO model, src_artifdata)
        ! z_exp has been calculated above in order to return min_lay_thckn as thickness
        ! of the lowest model layer
        DO jk = 1, nlevp1
          vct_a(jk)      = top_height*(2._wp/pi*ACOS(REAL(jk-1,wp)**stretch_fac/ &
            &              REAL(nlev,wp)**stretch_fac))**z_exp
          vct_b(jk)      = EXP(-vct_a(jk)/5000._wp)
        ENDDO
      ELSE
        ! use third-order polynomial
        x1 = (2._wp*stretch_fac-1._wp)*min_lay_thckn
        b = (top_height-x1/6._wp*nlev**3-(min_lay_thckn-x1/6._wp)*nlev)/&
            (nlev**2-1._wp/3._wp*nlev**3-2._wp/3._wp*nlev)
        a = (x1-2._wp*b)/6._wp
        c = min_lay_thckn-(a+b)
        DO jk = 1, nlevp1
          jkr = REAL(nlevp1-jk,wp)
          vct_a(jk) = a*jkr**3 + b*jkr**2 + c*jkr
          vct_b(jk) = EXP(-vct_a(jk)/5000._wp)
        ENDDO
      ENDIF
      ! Apply additional limitation on layer thickness in the middle and upper troposphere if the paramter
      ! max_lay_thckn is specified appropriately
      IF (max_lay_thckn > 2._wp*min_lay_thckn .AND. max_lay_thckn < 0.5_wp*htop_thcknlimit) THEN
        jk1 = 0
        DO jk = 1, nlev
          dvct(jk) = vct_a(jk) - vct_a(jk+1)
          IF (dvct(jk) > max_lay_thckn) jk1 = jk ! lowest layer in which the original layer thickness exceeds
                                                 ! the specified value
        ENDDO
        jks = 0
        jk2 = 0
        zvcta(nlevp1) = 0._wp
        DO jk = nlev, 1, -1
          IF (zvcta(jk+1) < htop_thcknlimit) THEN
            zvcta(jk) = zvcta(jk+1)+MIN(max_lay_thckn,dvct(jk))
          ELSE IF (jk2 == 0) THEN
            jk2 = jk+1
            jks = MAX(0,jk1-jk2)  ! shift layers from which thicknesses are taken downward in order to prevent sudden jumps
            zvcta(jk) = zvcta(jk+1)+dvct(jk+jks)
          ELSE
            zvcta(jk) = zvcta(jk+1)+dvct(jk+jks)
          ENDIF
        ENDDO
        IF (jks == 0) THEN ! either jk1 < htop_thcknlimit, which means that the thickness limiter has nothing to do,
                           ! or htop_thcknlimit is larger than the provisional model top;
          stretchfac = 1   ! in the latter case, the model top height is reset, overriding the top_height parameter
                           !
        ELSE               ! stretch remaining model levels such as to retain the original model top height
          stretchfac = (vct_a(1)-(zvcta(jk2)+REAL(jk2-1,wp)*max_lay_thckn))/&
                       (zvcta(1)-(zvcta(jk2)+REAL(jk2-1,wp)*max_lay_thckn))
        ENDIF
        DO jk = nlev, 1, -1
          IF (vct_a(jk+1) < htop_thcknlimit) THEN
            vct_a(jk) = vct_a(jk+1)+MIN(max_lay_thckn,dvct(jk))
          ELSE
            vct_a(jk) = vct_a(jk+1)+max_lay_thckn+(dvct(jk+jks)-max_lay_thckn)*stretchfac
          ENDIF
        ENDDO

        ! Try to apply additional smoothing on the stretching factor above the constant-thickness layer
        IF (stretchfac /= 1._wp .AND. jk1 < nlev-3) THEN
          DO jk = nlev, 1, -1
            IF (zvcta(jk+1) < htop_thcknlimit) THEN
              zvcta(jk) = vct_a(jk)
            ELSE
              zdvct = MIN(1.025_wp*(vct_a(jk)-vct_a(jk+1)), 1.025_wp*(zvcta(jk1+1)-zvcta(jk1+2))/ &
                         (zvcta(jk1+2)-zvcta(jk1+3))*(zvcta(jk+1)-zvcta(jk+2)) )
              zvcta(jk) = MIN(vct_a(jk),zvcta(jk+1)+zdvct)
            ENDIF
          ENDDO
          IF (zvcta(1) == vct_a(1)) THEN
            vct_a(1:2) = zvcta(1:2)
            vct_a(jk2+1:nlev) = zvcta(jk2+1:nlev)
            DO jk = 3, jk2
              vct_a(jk) = 0.5_wp*(zvcta(jk-1)+zvcta(jk+1))
            ENDDO
          ENDIF
        ENDIF

      ENDIF
    ELSE
     ! Use constant layer thicknesses determined by nlev and top_height
      DO jk = 1, nlevp1
        vct_a(jk) = top_height*(REAL(nlevp1,wp)-REAL(jk,wp))/REAL(nlev,wp)
        vct_b(jk) = EXP(-vct_a(jk)/5000._wp)
      ENDDO
    ENDIF

  END SUBROUTINE init_sleve_coord


  !---------------------------------------------------------------------------
  !> Utility routine: computes some values based on vct_a, vct_b
  !
  SUBROUTINE prepare_sleve_coord(nlev, vct_a, vct_b, vct, nflatlev)

    INTEGER,  INTENT(IN)  :: nlev  !< number of full levels
    REAL(wp), INTENT(IN)  :: vct_a(:), vct_b(:)
    REAL(wp), INTENT(OUT) :: vct(:)
    INTEGER,  INTENT(OUT) :: nflatlev(:)

    INTEGER  :: jk, nflat
    INTEGER  :: nlevp1        !< number of full and half levels

    ! number of vertical levels
    nlevp1 = nlev+1

    !------------------------------
    ! derived quantities: input: vct_a, vct_b -> output: vct, nflat, nflatlev

    IF (min_lay_thckn > 0.01_wp) THEN
      DO jk = 1, nlevp1
        vct(jk)        = vct_a(jk)
        vct(jk+nlevp1) = vct_b(jk)
      ENDDO
    ELSE
      ! Use constant layer thicknesses determined by nlev and top_height
      DO jk = 1, nlevp1
        vct(jk)        = vct_a(jk)
        vct(jk+nlevp1) = vct_b(jk)
      ENDDO
    ENDIF

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

    IF (msg_level >= 7) THEN
      CALL print_vcoord_info(nlev, vct_a, "init_sleve_coord")
    ENDIF

!!$    CALL message('mo_nh_init_utils: init_sleve_coord', '')
!!$
!!$    IF (msg_level >= 7) THEN
!!$     WRITE(message_text,'(a)') 'Nominal heights of coordinate half levels and layer thicknesses (m):'
!!$        CALL message('', TRIM(message_text))
!!$
!!$      DO jk = 1, nlevp1
!!$       WRITE(message_text,'(a,i4,2F12.3)') 'jk, vct_a, dvct: ',jk, vct_a(jk), vct_a(jk)-vct_a(MIN(jk+1,nlevp1))
!!$        CALL message('', TRIM(message_text))
!!$      ENDDO
!!$    ENDIF

  END SUBROUTINE prepare_sleve_coord


  !----------------------------------------------------------------------------
  !>
  !! Print some information about the vertical coordinate
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2011-07-01)
  !! Modification by Daniel Reinert, DWD (2017-01-12)
  !! - encapsulate into subroutine
  !!
  SUBROUTINE print_vcoord_info(nlev, vct_a, vcoord_type)
    !
    INTEGER,  INTENT(IN)  :: nlev  !< number of full levels
    REAL(wp), INTENT(IN)  :: vct_a(:)
    CHARACTER(len=MAX_CHAR_LENGTH) :: vcoord_type
    !
    ! local
    INTEGER :: jk
    INTEGER :: nlevp1        !< number of full and half levels

    ! number of vertical levels
    nlevp1 = nlev+1

    CALL message(TRIM(vcoord_type), 'Coordinate setup finished')

    WRITE(message_text,'(a)') 'Nominal heights of coordinate half levels and layer thicknesses (m):'
      CALL message('', TRIM(message_text))

     DO jk = 1, nlevp1
       WRITE(message_text,'(a,i4,2F12.3)') 'jk, vct_a, dvct: ',jk, vct_a(jk), vct_a(jk)-vct_a(MIN(jk+1,nlevp1))
       CALL message('', TRIM(message_text))
     ENDDO

  END SUBROUTINE print_vcoord_info

  !----------------------------------------------------------------------------
  !>
  !! Computes the 3D vertical coordinate fields for the nonhydrostatic model.
  !! (was originally included in subroutine set_nh_metrics but has been
  !! encapsulated because IFS2ICON needs the coordinate fields as input.
  !! Note: this routine is supposed to be used for cells and vertices.
  !! Therefore, field dimensions are passed.
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2011-07-01)
  !!
  SUBROUTINE init_vert_coord(vct_a, vct_b, topo, topo_smt, z3d_i,     &
                             nlev, nblks, npromz, nshift, nflat)

    ! Input parameters:
    REAL(wp), INTENT(IN) :: vct_a(:), vct_b(:)

    INTEGER, INTENT(IN) :: nlev, nblks, & ! field dimensions
                           npromz,      & ! length of last block
                           nshift,      & ! shift parameter for vertical nesting
                           nflat          ! index below which model levels are flat

    ! Input fields: "Normal" and smooth topography
    REAL(wp),  INTENT(IN) :: topo    (nproma,nblks), &
                             topo_smt(nproma,nblks)


    ! Output fields: 3D coordinate fields at interface and main levels
    REAL(wp),  INTENT(OUT) :: z3d_i(nproma,nlev+1,nblks)

    INTEGER :: jc, jk, jk1, jb, nlen, nlevp1, ierr(nblks), nerr, ktop_thicklimit(nproma)
    REAL(wp) :: z_fac1, z_fac2, z_topo_dev(nproma), min_lay_spacing, &
                dvct, dvct1, dvct2, minrat1, minrat2, wfac, dz1, dz2, dz3, dzr
    !-------------------------------------------------------------------------

    nlevp1  = nlev+1
    dvct1   = 100._wp
    minrat1 = 1._wp/3._wp ! minimum relative layer thickness for nominal thicknesses <= dvct1 (in m)
    dvct2   = 500._wp
    minrat2 = 0.5_wp  ! minimum relative layer thickness for a nominal thickness of dvct2
    ierr(:) = 0

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, jk1, z_fac1, z_fac2, z_topo_dev, min_lay_spacing,&
!$OMP dvct, wfac, ktop_thicklimit, dz1, dz2, dz3, dzr) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1,nblks

      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        z3d_i(nlen+1:nproma,:,jb) = 0._wp
     ENDIF

     z3d_i(1:nlen,nlevp1,jb) = topo(1:nlen,jb)
     ktop_thicklimit(:) = nlevp1

     ! vertical interface height
     IF (ivctype == 1 .OR. decay_scale_1 >= 0.5_wp*top_height) THEN ! hybrid Gal-Chen coordinate
       DO jk = 1, nlev
         jk1 = jk + nshift
         z3d_i(1:nlen,jk,jb) = vct_a(jk1) + vct_b(jk1)*z3d_i(1:nlen,nlevp1,jb)
       ENDDO
     ELSE IF (ivctype == 2) THEN ! SLEVE coordinate (Leuenberger et al. MWR 2010)
       DO jk = 1, nflat
         jk1 = jk + nshift
         z3d_i(1:nlen,jk,jb) = vct_a(jk1)
       ENDDO
       DO jk = nflat + 1, nlev
         jk1 = jk + nshift
         ! Scaling factors for large-scale and small-scale topography
         z_fac1 = SINH((top_height/decay_scale_1)**decay_exp - &
                       (vct_a(jk1)/decay_scale_1)**decay_exp)/ &
                  SINH((top_height/decay_scale_1)**decay_exp)
         z_fac2 = SINH((top_height/decay_scale_2)**decay_exp - &
                       (vct_a(jk1)/decay_scale_2)**decay_exp)/ &
                  SINH((top_height/decay_scale_2)**decay_exp)

         ! Small-scale topography (i.e. full topo - smooth topo)
         z_topo_dev(1:nlen) = topo(1:nlen,jb) - topo_smt(1:nlen,jb)

         z3d_i(1:nlen,jk,jb) = vct_a(jk1) + topo_smt(1:nlen,jb)*z_fac1 + &
            z_topo_dev(1:nlen)*z_fac2
       ENDDO
       ! Ensure that layer thicknesses are not too small; this would potentially cause
       ! instabilities in vertical advection
       DO jk = nlev, 1, -1
         jk1 = jk + nshift
         dvct = vct_a(jk1) - vct_a(jk1+1)
         IF (dvct < dvct1) THEN ! limit layer thickness to minrat1 times its nominal value
           min_lay_spacing = minrat1*dvct
         ELSE IF (dvct < dvct2) THEN ! limitation factor changes from minrat1 to minrat2
           wfac = ((dvct2-dvct)/(dvct2-dvct1))**2
           min_lay_spacing = (minrat1*wfac + minrat2*(1._wp-wfac))*dvct
         ELSE ! limitation factor decreases again
           min_lay_spacing = minrat2*dvct2*(dvct/dvct2)**(1._wp/3._wp)
         ENDIF
         min_lay_spacing = MAX(min_lay_spacing,MIN(50._wp,min_lay_thckn))
         DO jc = 1, nlen
           IF (z3d_i(jc,jk+1,jb)+min_lay_spacing > z3d_i(jc,jk,jb)) THEN
             z3d_i(jc,jk,jb) = z3d_i(jc,jk+1,jb) + min_lay_spacing
             ktop_thicklimit(jc) = jk
           ENDIF
         ENDDO
       ENDDO
       ! Smooth layer thickness ratios in the transition layer of columns where the thickness limiter has been active
       DO jc = 1, nlen
         jk = ktop_thicklimit(jc)
         IF (jk <= nlev-2 .AND. jk >= 4) THEN
           ! TODO : array access with subscript (jk-3)=0
           dz1 = z3d_i(jc,jk+1,jb)-z3d_i(jc,jk+2,jb)
           dz2 = z3d_i(jc,jk-3,jb)-z3d_i(jc,jk-2,jb)
           dzr = (dz2/dz1)**0.25_wp ! stretching factor
           dz3 = (z3d_i(jc,jk-2,jb)-z3d_i(jc,jk+1,jb))/(dzr*(1._wp+dzr*(1._wp+dzr)))
           z3d_i(jc,jk  ,jb) = MAX(z3d_i(jc,jk  ,jb), z3d_i(jc,jk+1,jb) + dz3*dzr)
           z3d_i(jc,jk-1,jb) = MAX(z3d_i(jc,jk-1,jb), z3d_i(jc,jk  ,jb) + dz3*dzr*dzr)
         ENDIF
       ENDDO
       ! Check if level nflat is still flat
       IF (ANY(z3d_i(1:nlen,nflat,jb) /= vct_a(nflat+nshift))) ierr(jb) = 1
       ! Check also if ktop_thicklimit is sufficiently far away from the model top
       IF (nlev > 6 .AND. ANY(ktop_thicklimit(1:nlen) <= 3)) ierr(jb) = ierr(jb) + 1
     ENDIF

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   nerr = SUM(ierr(1:nblks))
   IF (nerr > 0) CALL finish ('init_vert_coord: ', &
      'flat_height in sleve_nml or model top is too low')

  END SUBROUTINE init_vert_coord


END MODULE mo_init_vgrid

