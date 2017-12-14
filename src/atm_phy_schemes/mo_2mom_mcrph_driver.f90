!!==============================================================================
!!
!! Two-moment mixed-phase bulk microphysics
!!
!! original version by Axel Seifert, May 2003
!! with modifications by Ulrich Blahak, August 2007
!!
!!==============================================================================
!!
!!
!! @par Revision History
!! Ported into ICON from UCLA-LES by Anurag Dipankar (2013-12-15)
!!
!!==============================================================================
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!==============================================================================

MODULE mo_mcrph_sb

!------------------------------------------------------------------------------
!>
!! Description:
!!
!!   The subroutines in the module "gscp" calculate the rates of change of
!!   temperature, cloud condensate and water vapor due to cloud microphysical
!!   processes related to the formation of grid scale clouds and precipitation.
!!   In the COSMO model the microphysical subroutines are either
!!   called from "organize_gscp" or from "organize_physics" itself.
!!
!==============================================================================
!
! Declarations:
!
! Modules used:

!------------------------------------------------------------------------------
! Microphysical constants and variables
!------------------------------------------------------------------------------

!$ser verbatim USE mo_ser_2mom,           ONLY: serialize_input, &
!$ser&                                          serialize_output
USE mo_kind,                 ONLY: wp
USE mo_math_constants,       ONLY: pi
USE mo_vertical_coord_table, ONLY: vct_a
USE mo_physical_constants,   ONLY: &
    rhoh2o,        & ! density of liquid water
    alv,           & ! latent heat of vaporization
    als,           & ! latent heat of sublimation
    cpdr  => rcpd, & ! (spec. heat of dry air at constant press)^-1
    cvdr  => rcvd    ! (spec. heat of dry air at const vol)^-1

USE mo_exception,      ONLY: finish, message, message_text

USE mo_2mom_mcrph_main,     ONLY:                                &
     &                       clouds_twomoment,                   &
     &                       atmosphere, particle, particle_frozen, particle_lwf, &
     &                       rain_coeffs, ice_coeffs, snow_coeffs, graupel_coeffs, hail_coeffs, &
     &                       ccn_coeffs, in_coeffs,                     &
     &                       init_2mom_scheme, init_2mom_scheme_once,   &
     &                       qnc_const

USE mo_2mom_mcrph_processes, ONLY:                                &
     &                       sedi_vel_rain, sedi_vel_sphere,     &
     &                       sedi_icon_rain, sedi_icon_sphere, sedi_icon_sphere_lwf, &
     &                       particle_meanmass, &
     &                       q_crit

USE mo_2mom_mcrph_util, ONLY:                            &
     &                       gfct,                       &  ! Gamma function (becomes intrinsic in Fortran2008)
     &                       ltabdminwgg,                &
     &                       init_dmin_wetgrowth,        &
     &                       init_dmin_wg_gr_ltab_equi

USE mo_2mom_prepare, ONLY: prepare_twomoment, post_twomoment

!==============================================================================

IMPLICIT NONE
PUBLIC

CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_mcrph_driver'
INTEGER,          PARAMETER :: dbg_level = 25                   ! level for debug prints

  ! .. exponents for simple height dependency of terminal fall velocity
  REAL(wp), PARAMETER :: rho_vel    = 0.4e0_wp    !..exponent for density correction
  REAL(wp), PARAMETER :: rho_vel_c  = 0.2e0_wp    !..for cloud droplets
  REAL(wp), PARAMETER :: rho0    = 1.225_wp     !..Norm-Luftdichte

INTEGER :: cloud_type, ccn_type

INTEGER, PARAMETER :: cloud_type_default_gscp4 = 2103, ccn_type_gscp4 = 1
INTEGER, PARAMETER :: cloud_type_default_gscp5 = 2603, ccn_type_gscp5 = 8

! AS: For gscp=4 use 2103 with ccn_type = 1 (HDCP2 IN and CCN schemes)
!     For gscp=5 use 2603 with ccn_type = 8 (PDA ice nucleation and Segal&Khain CCN activation)

! AS: Runs without hail, e.g, 1503 are buggy and give a segmentation fault.
!     So far I was not able to identify the problem, needs more detailed debugging.

CONTAINS

  !==============================================================================
  !
  ! Two-moment mixed-phase bulk microphysics
  !
  ! original version by Axel Seifert, May 2003
  ! with modifications by Ulrich Blahak, August 2007

  !==============================================================================
  SUBROUTINE two_moment_mcrph(            &
                       isize,             & ! in: array size
                       ke,                & ! in: end level/array size
                       is,                & ! in: start index, optional
                       ie,                & ! in: end index, optional
                       ks,                & ! in: start index vertical , optional
                       dt,                & ! in: time step
                       dz,                & ! in: vertical layer thickness
                       rho,               & ! in: density
                       pres,              & ! in: pressure
                       qv,                & ! inout: specific humidity
                       qc, qnc,           & ! inout: cloud water
                       qr, qnr,           & ! inout: rain
                       qi, qni,           & ! inout: ice
                       qs, qns,           & ! inout: snow
                       qg, qng, qgl,      & ! inout: graupel
                       qh, qnh, qhl,      & ! inout: hail
                       nccn,              & ! inout: ccn
                       ninpot,            & ! inout: potential ice nuclei
                       ninact,            & ! inout: activated ice nuclei
                       tk,                & ! inout: temp
                       w,                 & ! inout: w
                       prec_r,            & ! inout: precip rate rain
                       prec_i,            & ! inout: precip rate ice
                       prec_s,            & ! inout: precip rate snow
                       prec_g,            & ! inout: precip rate graupel
                       prec_h,            & ! inout: precip rate hail
                       dtemp,             & ! inout: opt. temp increment
                       msg_level,         & ! in: msg_level
                       l_cv          )      ! in: switch for cv/cp


    ! Declare variables in argument list

    INTEGER,            INTENT (IN)  :: isize, ke    ! grid sizes
    INTEGER,  OPTIONAL, INTENT (IN)  :: is, ie, ks   ! start/end indices

    REAL(wp), INTENT (IN)            :: dt           ! time step

    ! Dynamical core variables
    REAL(wp), DIMENSION(:,:), INTENT(IN), TARGET :: dz, rho, pres, w

    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET :: tk

    ! Microphysics variables
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) , TARGET :: &
         qv, qc, qnc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh, ninact

    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
         &               qgl, qhl

    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
         &               nccn, ninpot

    ! Precip rates, vertical profiles
    REAL(wp), DIMENSION(:), INTENT (INOUT)  :: &
         &               prec_r, prec_i, prec_s, prec_g, prec_h

    REAL(wp), OPTIONAL, INTENT (INOUT)  :: dtemp(:,:)

    INTEGER,  INTENT (IN)             :: msg_level
    LOGICAL,  OPTIONAL,  INTENT (IN)  :: l_cv

    ! ... Variables which are global in module_2mom_mcrph_main

    REAL(wp), TARGET, DIMENSION(isize,ke) ::        &
         &  qnc_dummy,     & ! cloud droplet number
         &  rhocorr,       & ! density dependency of particle fall speed
         &  rhocld           ! density dependency of particle fall speed for cloud droplets

    REAL(wp) :: q_liq_old(isize,ke), q_vap_old(isize,ke)  ! to store old values for latent heat calc

    INTEGER  :: its,ite,kts,kte
    INTEGER  :: ii,kk
    INTEGER  :: ntsedi     ! for sedimentation sub stepping

    REAL(wp) :: q_liq_new,q_vap_new
    REAL(wp) :: zf,hlp,dtemp_loc
    REAL(wp) :: convliq,convice
    REAL(wp), PARAMETER :: tau_inact =  600.  ! relaxation time scale for activated IN number density
    REAL(wp), PARAMETER :: tau_inpot = 1800.  ! relaxation time scale for potential IN number density
    REAL(wp) :: in_bgrd            ! background profile of IN number density
    REAL(wp) :: z_heat_cap_r       ! reciprocal of cpdr or cvdr (depending on l_cv)
    REAL(wp) :: rdz(isize,ke), rho_r(isize,ke)

    LOGICAL :: lprogccn, lprogin, lprogmelt

    LOGICAL, PARAMETER :: debug     = .false.       !
    LOGICAL, PARAMETER :: clipping  = .true.        ! not really necessary, just for cleanup

    CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_mcrph_driver'

    LOGICAL, PARAMETER :: explicit_solver = .true.  ! explicit or semi-implicit solver

    ! These structures include the pointers to the model arrays (which are automatic arrays
    ! of this driver subroutine). These structures live only for one time step and are
    ! different for the OpenMP threads. In contrast, the types like rain_coeffs, ice_coeffs,
    ! etc. that are declared in mo_2mom_mcrph_main live for the whole runtime and may include
    ! coefficients that are calculated once during initialization (and there is only one per
    ! mpi thread).
    TYPE(atmosphere)           :: atmo

    ! These are the fundamental hydrometeor particle variables for the two-moment scheme
    TYPE(particle), target          :: cloud_hyd, rain_hyd
    TYPE(particle_frozen), target   :: ice_frz, snow_frz, graupel_frz, hail_frz
    TYPE(particle_lwf), target      :: graupel_lwf, hail_lwf

    ! Pointers to the derived types that are actually needed
    CLASS(particle), pointer        :: cloud, rain
    CLASS(particle_frozen), pointer :: ice, snow, graupel, hail

    INTEGER :: ik_slice(4)

    !$ser verbatim CALL serialize_input(dz, rho, pres, w)
    !$ser verbatim CALL serialize_output(dz, rho, pres, w)

    IF (PRESENT(nccn)) THEN
       cloud_type = cloud_type_default_gscp5 + 10 * ccn_type
    ELSE
       cloud_type = cloud_type_default_gscp4 + 10 * ccn_type
    END IF

    cloud => cloud_hyd
    rain  => rain_hyd
    ice   => ice_frz
    snow  => snow_frz

    IF (PRESENT(qgl)) THEN
       graupel => graupel_lwf   ! gscp=4,5,6
       hail => hail_lwf
    ELSE
       graupel => graupel_frz   ! gscp=7
       hail => hail_frz
    END IF

    ! start/end indices
    IF (PRESENT(is)) THEN
      its = is
    ELSE
      its = 1
    END IF

    IF (PRESENT(ie)) THEN
      ite = ie
    ELSE
      ite = isize
    END IF

    IF (PRESENT(ks)) THEN
      kts = ks
    ELSE
      kts = 1
    END IF
    kte = ke

    ! inverse of vertical layer thickness
    rdz(its:ite,kts:kte) = 1._wp / dz(its:ite,kts:kte)

    lprogccn  = PRESENT(nccn)
    lprogin   = PRESENT(ninpot)
    lprogmelt = PRESENT(qgl)

    IF (clipping) THEN
       WHERE(qr(its:ite,kts:kte) < 0.0_wp) qr(its:ite,kts:kte) = 0.0_wp
       WHERE(qi(its:ite,kts:kte) < 0.0_wp) qi(its:ite,kts:kte) = 0.0_wp
       WHERE(qs(its:ite,kts:kte) < 0.0_wp) qs(its:ite,kts:kte) = 0.0_wp
       WHERE(qg(its:ite,kts:kte) < 0.0_wp) qg(its:ite,kts:kte) = 0.0_wp
       WHERE(qh(its:ite,kts:kte) < 0.0_wp) qh(its:ite,kts:kte) = 0.0_wp
       IF (lprogmelt) THEN
          WHERE(qgl(its:ite,kts:kte) < 0.0_wp) qgl(its:ite,kts:kte) = 0.0_wp
          WHERE(qhl(its:ite,kts:kte) < 0.0_wp) qhl(its:ite,kts:kte) = 0.0_wp
       END IF
    END IF

    ! indices as used in two-moment scheme
    ik_slice(1) = its
    ik_slice(2) = ite
    ik_slice(3) = kts
    ik_slice(4) = kte

    IF (PRESENT(l_cv)) THEN
      IF (l_cv) THEN
        z_heat_cap_r = cvdr
      ELSE
        z_heat_cap_r = cpdr
      ENDIF
    ELSE
       z_heat_cap_r = cpdr
    ENDIF

    IF (msg_level>dbg_level) CALL message(TRIM(routine),'')

    IF (msg_level>dbg_level)THEN
       WRITE (message_text,'(1X,A,I4,3(A,L2))') &
            & "cloud_type = ",cloud_type,", lprogccn = ",lprogccn,", lprogin = ",lprogin,", lprogmelt = ",lprogmelt
       CALL message(TRIM(routine),TRIM(message_text))
    END IF

    ! time step for two-moment microphysics is the same as all other fast physics
    IF (msg_level>dbg_level) CALL message(TRIM(routine), "prepare variable for 2mom")

    DO kk = kts, kte
       DO ii = its, ite

          ! ... 1/rho is used quite often
          rho_r(ii,kk) = 1.0 / rho(ii,kk)

          ! ... height dependency of terminal fall velocities
          hlp = log(rho(ii,kk)/rho0)
          rhocorr(ii,kk) = exp(-rho_vel*hlp)
          rhocld(ii,kk)  = exp(-rho_vel_c*hlp)

       END DO
    END DO

    ! .. set the particle types, but no calculations
    CALL init_2mom_scheme(cloud,rain,ice,snow,graupel,hail)

    ! .. convert to densities and set pointerns to two-moment module
    !    (pointers are used to avoid passing everything explicitly by argument and
    !     to avoid local allocates within the OpenMP-loop, and keep everything on stack)

    CALL prepare_twomoment(atmo, cloud, rain, ice, snow, graupel, hail, &
         rho, rhocorr, rhocld, pres, w, tk, &
         nccn, ninpot, ninact, &
         qv, qc, qnc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh, qgl, qhl, &
         lprogccn, lprogin, lprogmelt, its, ite, kts, kte)

    IF (msg_level>dbg_level) CALL message(TRIM(routine)," calling clouds_twomoment")

    IF (explicit_solver) then

       ! ... save old variables for latent heat calculation
       q_vap_old(its:ite,kts:kte) = qv(its:ite,kts:kte)

       if (lprogmelt) then
          q_liq_old(its:ite,kts:kte) = qc(its:ite,kts:kte) + qr(its:ite,kts:kte)  &
               &                     + qgl(its:ite,kts:kte) + qhl(its:ite,kts:kte)
       else
          q_liq_old(its:ite,kts:kte) = qc(its:ite,kts:kte) + qr(its:ite,kts:kte)
       end if

       ! .. this subroutine calculates all the microphysical sources and sinks
       CALL clouds_twomoment(ik_slice, dt, lprogin, &
            atmo, cloud, rain, ice, snow, graupel, hail, ninact, nccn, ninpot)

       IF (lprogccn) THEN
          !WHERE(qc > 0.0_wp)  cloud%n = 5000e6_wp
          WHERE(qc(its:ite,kts:kte) == 0.0_wp) cloud%n(its:ite,kts:kte) = 0.0_wp
       END IF

       ! .. latent heat term for temperature equation
       convice = z_heat_cap_r * als
       convliq = z_heat_cap_r * (alv-als)
       DO kk = kts, kte
          DO ii = its, ite

             ! .. new variables
             q_vap_new = qv(ii,kk)
             if (lprogmelt) then
                q_liq_new = qr(ii,kk) + qc(ii,kk) + qgl(ii,kk) + qhl(ii,kk) 
             else
                q_liq_new = qr(ii,kk) + qc(ii,kk)
             end if

             ! .. update temperature
             dtemp_loc  = - convice * rho_r(ii,kk) * (q_vap_new - q_vap_old(ii,kk))  &
                  &       + convliq * rho_r(ii,kk) * (q_liq_new - q_liq_old(ii,kk))

             tk(ii,kk) = tk(ii,kk) + dtemp_loc

             IF(PRESENT(dtemp)) &
               dtemp(ii,kk) = dtemp_loc

          ENDDO
       ENDDO

       ! .. if we solve explicitly, then sedimentation is done here after microphysics
       CALL sedimentation_explicit()

    ELSE

       ! .. semi-implicit solver includes microphysics and sedimentation
       if (lprogin) THEN
          CALL message(TRIM(routine)," ERROR: gscp=5 not implemented for implicit solver")
          CALL finish(TRIM(routine),'Error in two_moment_mcrph')
       ELSE
          CALL clouds_twomoment_implicit ()
       END IF
    END IF

    ! .. check for negative values
    IF (debug) CALL check_clouds()

    ! .. convert back and nullify two-moment pointers
    CALL post_twomoment(atmo, cloud, rain, ice, snow, graupel, hail, &
         rho_r, qnc, nccn, ninpot, ninact, &
         qv, qc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh, qgl, qhl,  &
         lprogccn, lprogin, lprogmelt, its, ite, kts, kte)

    IF (clipping) THEN
       WHERE(qr(its:ite,kts:kte) < 0.0_wp) qr(its:ite,kts:kte) = 0.0_wp
       WHERE(qi(its:ite,kts:kte) < 0.0_wp) qi(its:ite,kts:kte) = 0.0_wp
       WHERE(qs(its:ite,kts:kte) < 0.0_wp) qs(its:ite,kts:kte) = 0.0_wp
       WHERE(qg(its:ite,kts:kte) < 0.0_wp) qg(its:ite,kts:kte) = 0.0_wp
       WHERE(qh(its:ite,kts:kte) < 0.0_wp) qh(its:ite,kts:kte) = 0.0_wp
       WHERE(qnr(its:ite,kts:kte) < 0.0_wp) qnr(its:ite,kts:kte) = 0.0_wp
       WHERE(qni(its:ite,kts:kte) < 0.0_wp) qni(its:ite,kts:kte) = 0.0_wp
       WHERE(qns(its:ite,kts:kte) < 0.0_wp) qns(its:ite,kts:kte) = 0.0_wp
       WHERE(qng(its:ite,kts:kte) < 0.0_wp) qng(its:ite,kts:kte) = 0.0_wp
       WHERE(qnh(its:ite,kts:kte) < 0.0_wp) qnh(its:ite,kts:kte) = 0.0_wp
       IF (lprogmelt) THEN
          WHERE(qgl(its:ite,kts:kte) < 0.0_wp) qgl(its:ite,kts:kte) = 0.0_wp
          WHERE(qhl(its:ite,kts:kte) < 0.0_wp) qhl(its:ite,kts:kte) = 0.0_wp
       END IF
    END IF

    IF (lprogccn) THEN
      DO kk=kts,kte
        zf = 0.5_wp*(vct_a(kk)+vct_a(kk+1))
        DO ii=its,ite
          !..reset nccn for cloud-free grid points to background profile
          IF (qc(ii,kk) .LE. q_crit) THEN
            IF(zf > ccn_coeffs%z0) THEN
              nccn(ii,kk) = MAX(nccn(ii,kk),ccn_coeffs%Ncn0 &
                   * EXP((ccn_coeffs%z0 - zf)*(1._wp/ccn_coeffs%z1e)))
            ELSE
              nccn(ii,kk) = MAX(nccn(ii,kk),ccn_coeffs%Ncn0)
            END IF
          END IF
        END DO
      END DO
    END IF

    DO kk=kts,kte
      DO ii=its,ite
        !..relaxation of activated IN number density to zero
        IF(qi(ii,kk) == 0) THEN
          ninact(ii,kk) = ninact(ii,kk) - ninact(ii,kk)*(1._wp/tau_inact)*dt
        END IF
      END DO
    END DO

    IF (lprogin) THEN
      DO kk=kts,kte
        DO ii=its,ite
          zf = 0.5_wp*(vct_a(kk)+vct_a(kk+1))
          !..relaxation of potential IN number density to background profile
          IF(zf > in_coeffs%z0) THEN
            in_bgrd = in_coeffs%N0*EXP((in_coeffs%z0 - zf)*(1._wp/in_coeffs%z1e))
          ELSE
            in_bgrd = in_coeffs%N0
          END IF
          ninpot(ii,kk) = ninpot(ii,kk) - (ninpot(ii,kk)-in_bgrd)*(1._wp/tau_inpot)*dt
        END DO
      END DO
    END IF

    IF (lprogccn) &
      WHERE(nccn(its:ite,kts:kte) < 35.0_wp) nccn(its:ite,kts:kte) = 35e6_wp

    WHERE(qc(its:ite,kts:kte) < 1.0e-12_wp) qnc(its:ite,kts:kte) = 0.0_wp
    WHERE(qr(its:ite,kts:kte) > 0.02_wp) qr(its:ite,kts:kte) = 0.02_wp
    WHERE(qi(its:ite,kts:kte) > 0.02_wp) qi(its:ite,kts:kte) = 0.02_wp
    WHERE(qs(its:ite,kts:kte) > 0.02_wp) qs(its:ite,kts:kte) = 0.02_wp
    WHERE(qg(its:ite,kts:kte) > 0.02_wp) qg(its:ite,kts:kte) = 0.02_wp
    WHERE(qh(its:ite,kts:kte) > 0.02_wp) qh(its:ite,kts:kte) = 0.02_wp

    IF (msg_level>dbg_level) CALL message(TRIM(routine), "two moment mcrph ends!")

    RETURN
    !
    ! end of driver routine, but many details are below in the contains-part of this subroutine
    !
  CONTAINS

    SUBROUTINE clouds_twomoment_implicit()
      !
      ! semi-implicit solver for sedimentation including microphysics, the same
      ! approach is used in the COSMO microphysics, e.g, hydci_pp
      ! (see COSMO documentation for details)
      !

      ! a few 1d arrays, maybe we can reduce this later or we keep them ...
      real(wp), dimension(isize) :: &
           & qr_flux_now,qr_flux_new,qr_sum,qr_flux_sum,vr_sedq_new,vr_sedq_now,qr_impl,qr_star,xr_now, &
           & nr_flux_now,nr_flux_new,nr_sum,nr_flux_sum,vr_sedn_new,vr_sedn_now,nr_impl,nr_star,xr_new, &
           & qs_flux_now,qs_flux_new,qs_sum,qs_flux_sum,vs_sedq_new,vs_sedq_now,qs_impl,qs_star,xs_now, &
           & ns_flux_now,ns_flux_new,ns_sum,ns_flux_sum,vs_sedn_new,vs_sedn_now,ns_impl,ns_star,xs_new, &
           & qg_flux_now,qg_flux_new,qg_sum,qg_flux_sum,vg_sedq_new,vg_sedq_now,qg_impl,qg_star,xg_now, &
           & ng_flux_now,ng_flux_new,ng_sum,ng_flux_sum,vg_sedn_new,vg_sedn_now,ng_impl,ng_star,xg_new, &
           & qh_flux_now,qh_flux_new,qh_sum,qh_flux_sum,vh_sedq_new,vh_sedq_now,qh_impl,qh_star,xh_now, &
           & nh_flux_now,nh_flux_new,nh_sum,nh_flux_sum,vh_sedn_new,vh_sedn_now,nh_impl,nh_star,xh_new, &
           & qi_flux_now,qi_flux_new,qi_sum,qi_flux_sum,vi_sedq_new,vi_sedq_now,qi_impl,qi_star,xi_now, &
           & ni_flux_now,ni_flux_new,ni_sum,ni_flux_sum,vi_sedn_new,vi_sedn_now,ni_impl,ni_star,xi_new

      REAL(wp), DIMENSION(isize,ke) :: rdzdt
      INTEGER :: i, ii, k, kk

      logical, parameter :: lmicro_impl = .true.  ! microphysics within semi-implicit sedimentation loop?

      convice = z_heat_cap_r * als
      convliq = z_heat_cap_r * (alv-als)

      if (.not.lmicro_impl) then

         ! ... save old variables for latent heat calculation
         q_vap_old(its:ite,kts:kte) = qv(its:ite,kts:kte)
         q_liq_old(its:ite,kts:kte) = qc(its:ite,kts:kte) + qr(its:ite,kts:kte)

         ! .. this subroutine calculates all the microphysical sources and sinks
         CALL clouds_twomoment(ik_slice, dt, lprogin, atmo, cloud, rain, &
              ice, snow, graupel, hail, ninact, nccn, ninpot)

         DO kk=kts,kte
            DO ii = its, ite
               ! .. latent heat term for temperature equation
               q_vap_new = qv(ii,kk)
               q_liq_new = qr(ii,kk) + qc(ii,kk)
               tk(ii,kk) = tk(ii,kk) - convice * rho_r(ii,kk) * (q_vap_new - q_vap_old(ii,kk))  &
                    &                + convliq * rho_r(ii,kk) * (q_liq_new - q_liq_old(ii,kk))
            ENDDO
         ENDDO

      end if

      ! clipping maybe not necessary
      WHERE(qr(its:ite,kts:kte) < 0.0_wp) qr(its:ite,kts:kte) = 0.0_wp
      WHERE(qi(its:ite,kts:kte) < 0.0_wp) qi(its:ite,kts:kte) = 0.0_wp
      WHERE(qs(its:ite,kts:kte) < 0.0_wp) qs(its:ite,kts:kte) = 0.0_wp
      WHERE(qg(its:ite,kts:kte) < 0.0_wp) qg(its:ite,kts:kte) = 0.0_wp
      WHERE(qh(its:ite,kts:kte) < 0.0_wp) qh(its:ite,kts:kte) = 0.0_wp
      WHERE(qnr(its:ite,kts:kte) < 0.0_wp) qnr(its:ite,kts:kte) = 0.0_wp
      WHERE(qni(its:ite,kts:kte) < 0.0_wp) qni(its:ite,kts:kte) = 0.0_wp
      WHERE(qns(its:ite,kts:kte) < 0.0_wp) qns(its:ite,kts:kte) = 0.0_wp
      WHERE(qng(its:ite,kts:kte) < 0.0_wp) qng(its:ite,kts:kte) = 0.0_wp
      WHERE(qnh(its:ite,kts:kte) < 0.0_wp) qnh(its:ite,kts:kte) = 0.0_wp

      rdzdt(its:ite,kts:kte) = 0.5_wp * rdz(its:ite,kts:kte) * dt

      qr_flux_now(:) = 0.0_wp
      nr_flux_now(:) = 0.0_wp
      qr_flux_new(:) = 0.0_wp
      nr_flux_new(:) = 0.0_wp

      qi_flux_now(:) = 0.0_wp
      ni_flux_now(:) = 0.0_wp
      qi_flux_new(:) = 0.0_wp
      ni_flux_new(:) = 0.0_wp

      qs_flux_now(:) = 0.0_wp
      ns_flux_now(:) = 0.0_wp
      qs_flux_new(:) = 0.0_wp
      ns_flux_new(:) = 0.0_wp

      qg_flux_now(:) = 0.0_wp
      ng_flux_now(:) = 0.0_wp
      qg_flux_new(:) = 0.0_wp
      ng_flux_new(:) = 0.0_wp

      qh_flux_now(:) = 0.0_wp
      nh_flux_now(:) = 0.0_wp
      qh_flux_new(:) = 0.0_wp
      nh_flux_new(:) = 0.0_wp

      convice = z_heat_cap_r * als
      convliq = z_heat_cap_r * (alv-als)

      do i=its,ite
         vr_sedn_new(i) = rain%vsedi_min
         vi_sedn_new(i) = ice%vsedi_min
         vs_sedn_new(i) = snow%vsedi_min
         vg_sedn_new(i) = graupel%vsedi_min
         vh_sedn_new(i) = hail%vsedi_min
      end do

      ! here we simply assume that there is no cloud or precip in the uppermost level
      ! i.e. we start from kts+1 going down in physical space

!DIR$ IVDEP
      DO k=kts+1,kte

        do i=its,ite
           xr_now(i) = particle_meanmass(rain, qr(i,k),qnr(i,k))
           xi_now(i) = particle_meanmass(ice, qi(i,k),qni(i,k))
           xs_now(i) = particle_meanmass(snow, qs(i,k),qns(i,k))
           xg_now(i) = particle_meanmass(graupel, qg(i,k),qng(i,k))
           xh_now(i) = particle_meanmass(hail, qh(i,k),qnh(i,k))
        end do

        call sedi_vel_rain(rain,rain_coeffs,qr(:,k),xr_now,vr_sedn_now,vr_sedq_now,its,ite,qc(:,k))
        call sedi_vel_sphere(ice,ice_coeffs,qi(:,k),xi_now,vi_sedn_now,vi_sedq_now,its,ite)
        call sedi_vel_sphere(snow,snow_coeffs,qs(:,k),xs_now,vs_sedn_now,vs_sedq_now,its,ite)
        call sedi_vel_sphere(graupel,graupel_coeffs,qg(:,k),xg_now,vg_sedn_now,vg_sedq_now,its,ite)
        call sedi_vel_sphere(hail,hail_coeffs,qh(:,k),xh_now,vh_sedn_now,vh_sedq_now,its,ite)

        do i=its,ite

           ! .... rain ....
           vr_sedn_new(i) = 0.5 * (vr_sedn_now(i) + vr_sedn_new(i))
           vr_sedq_new(i) = 0.5 * (vr_sedq_now(i) + vr_sedq_new(i))

           ! qflux_new, nflux_new are the updated flux values from the level above
           ! qflux_now, nflux_now are here the old (current time step) flux values from the level above
           ! In COSMO-Docu  {...} =  flux_(k-1),new + flux_(k-1),start
           nr_flux_sum(i) = nr_flux_new(i) + nr_flux_now(i)
           qr_flux_sum(i) = qr_flux_new(i) + qr_flux_now(i)

           ! qflux_now, nflux_now are here overwritten with the current level
           nr_flux_now(i) = min(vr_sedn_now(i) * qnr(i,k), nr_flux_sum(i))  ! this is then passed to the level below
           qr_flux_now(i) = min(vr_sedq_now(i) * qr(i,k),  qr_flux_sum(i))  ! (loop dependency)
           nr_flux_now(i) = max(nr_flux_now(i),0.0_wp) ! maybe not necessary
           qr_flux_now(i) = max(qr_flux_now(i),0.0_wp) ! maybe not necessary

           nr_sum(i) = qnr(i,k) + rdzdt(i,k) * (nr_flux_sum(i) - nr_flux_now(i))
           qr_sum(i) = qr(i,k)  + rdzdt(i,k) * (qr_flux_sum(i) - qr_flux_now(i))

           nr_impl(i) = 1.0_wp/(1.0_wp + vr_sedn_new(i) * rdzdt(i,k))
           qr_impl(i) = 1.0_wp/(1.0_wp + vr_sedq_new(i) * rdzdt(i,k))

           nr_star(i) = nr_impl(i) * nr_sum(i)  ! intermediate values for calculating
           qr_star(i) = qr_impl(i) * qr_sum(i)  ! sources and sinks

           qnr(i,k) = nr_star(i)                ! overwrite array with intermediate
           qr(i,k)  = qr_star(i)                ! values to do micro processes on this level

           nr_sum(i) = nr_sum(i) - nr_star(i)   ! final time integration starts from sum-values
           qr_sum(i) = qr_sum(i) - qr_star(i)   ! but source/sinks work on star-values

           ! .... ice ....
           vi_sedn_new(i) = 0.5 * (vi_sedn_now(i) + vi_sedn_new(i))
           vi_sedq_new(i) = 0.5 * (vi_sedq_now(i) + vi_sedq_new(i))

           ! qflux_new, nflux_new are the updated flux values from the level above
           ! qflux_now, nflux_now are here the old (current time step) flux values from the level above
           ! In COSMO-Docu  {...} =  flux_(k-1),new + flux_(k-1),start
           ni_flux_sum(i) = ni_flux_new(i) + ni_flux_now(i)
           qi_flux_sum(i) = qi_flux_new(i) + qi_flux_now(i)

           ! qflux_now, nflux_now are here overwritten with the current level
           ni_flux_now(i) = min(vi_sedn_now(i) * qni(i,k), ni_flux_sum(i))  ! this is then passed to the level below
           qi_flux_now(i) = min(vi_sedq_now(i) * qi(i,k),  qi_flux_sum(i))  ! (loop dependency)
           ni_flux_now(i) = max(ni_flux_now(i),0.0_wp)
           qi_flux_now(i) = max(qi_flux_now(i),0.0_wp)

           ni_sum(i) = qni(i,k) + rdzdt(i,k) * (ni_flux_sum(i) - ni_flux_now(i))
           qi_sum(i) = qi(i,k)  + rdzdt(i,k) * (qi_flux_sum(i) - qi_flux_now(i))

           ni_impl(i) = 1.0_wp/(1.0_wp + vi_sedn_new(i) * rdzdt(i,k))
           qi_impl(i) = 1.0_wp/(1.0_wp + vi_sedq_new(i) * rdzdt(i,k))

           ni_star(i) = ni_impl(i) * ni_sum(i)  ! intermediate values for calculating
           qi_star(i) = qi_impl(i) * qi_sum(i)  ! sources and sinks

           qni(i,k) = ni_star(i)                ! overwrite array with intermediate
           qi(i,k)  = qi_star(i)                ! values to do micro processes on this level

           ni_sum(i) = ni_sum(i) - ni_star(i)   ! final time integration starts from sum-values
           qi_sum(i) = qi_sum(i) - qi_star(i)   ! but source/sinks work on star-values

           ! .... snow ....
           vs_sedn_new(i) = 0.5 * (vs_sedn_now(i) + vs_sedn_new(i))
           vs_sedq_new(i) = 0.5 * (vs_sedq_now(i) + vs_sedq_new(i))

           ! qflux_new, nflux_new are the updated flux values from the level above
           ! qflux_now, nflux_now are here the old (current time step) flux values from the level above
           ! In COSMO-Docu  {...} =  flux_(k-1),new + flux_(k-1),start
           ns_flux_sum(i) = ns_flux_new(i) + ns_flux_now(i)
           qs_flux_sum(i) = qs_flux_new(i) + qs_flux_now(i)

           ! qflux_now, nflux_now are here overwritten with the current level
           ns_flux_now(i) = min(vs_sedn_now(i) * qns(i,k), ns_flux_sum(i))  ! this is then passed to the level below
           qs_flux_now(i) = min(vs_sedq_now(i) * qs(i,k),  qs_flux_sum(i))  ! (loop dependency)
           ns_flux_now(i) = max(ns_flux_now(i),0.0_wp)
           qs_flux_now(i) = max(qs_flux_now(i),0.0_wp)

           ns_sum(i) = qns(i,k) + rdzdt(i,k) * (ns_flux_sum(i) - ns_flux_now(i))
           qs_sum(i) = qs(i,k)  + rdzdt(i,k) * (qs_flux_sum(i) - qs_flux_now(i))

           ns_impl(i) = 1.0_wp/(1.0_wp + vs_sedn_new(i) * rdzdt(i,k))
           qs_impl(i) = 1.0_wp/(1.0_wp + vs_sedq_new(i) * rdzdt(i,k))

           ns_star(i) = ns_impl(i) * ns_sum(i)  ! intermediate values for calculating
           qs_star(i) = qs_impl(i) * qs_sum(i)  ! sources and sinks

           qns(i,k) = ns_star(i)                ! overwrite array with intermediate
           qs(i,k)  = qs_star(i)                ! values to do micro processes on this level

           ns_sum(i) = ns_sum(i) - ns_star(i)   ! final time integration starts from sum-values
           qs_sum(i) = qs_sum(i) - qs_star(i)   ! but source/sinks work on star-values

           ! .... graupel ....
           vg_sedn_new(i) = 0.5 * (vg_sedn_now(i) + vg_sedn_new(i))
           vg_sedq_new(i) = 0.5 * (vg_sedq_now(i) + vg_sedq_new(i))

           ! qflux_new, nflux_new are the updated flux values from the level above
           ! qflux_now, nflux_now are here the old (current time step) flux values from the level above
           ! In COSMO-Docu  {...} =  flux_(k-1),new + flux_(k-1),start
           ng_flux_sum(i) = ng_flux_new(i) + ng_flux_now(i)
           qg_flux_sum(i) = qg_flux_new(i) + qg_flux_now(i)

           ! qflux_now, nflux_now are here overwritten with the current level
           ng_flux_now(i) = min(vg_sedn_now(i) * qng(i,k), ng_flux_sum(i))  ! this is then passed to the level below
           qg_flux_now(i) = min(vg_sedq_now(i) * qg(i,k),  qg_flux_sum(i))  ! (loop dependency)
           ng_flux_now(i) = max(ng_flux_now(i),0.0_wp)
           qg_flux_now(i) = max(qg_flux_now(i),0.0_wp)

           ng_sum(i) = qng(i,k) + rdzdt(i,k) * (ng_flux_sum(i) - ng_flux_now(i))
           qg_sum(i) = qg(i,k)  + rdzdt(i,k) * (qg_flux_sum(i) - qg_flux_now(i))

           ng_impl(i) = 1.0_wp/(1.0_wp + vg_sedn_new(i) * rdzdt(i,k))
           qg_impl(i) = 1.0_wp/(1.0_wp + vg_sedq_new(i) * rdzdt(i,k))

           ng_star(i) = ng_impl(i) * ng_sum(i)  ! intermediate values for calculating
           qg_star(i) = qg_impl(i) * qg_sum(i)  ! sources and sinks

           qng(i,k) = ng_star(i)                ! overwrite array with intermediate
           qg(i,k)  = qg_star(i)                ! values to do micro processes on this level

           ng_sum(i) = ng_sum(i) - ng_star(i)   ! final time integration starts from sum-values
           qg_sum(i) = qg_sum(i) - qg_star(i)   ! but source/sinks work on star-values

           ! .... hail ....
           vh_sedn_new(i) = 0.5 * (vh_sedn_now(i) + vh_sedn_new(i))
           vh_sedq_new(i) = 0.5 * (vh_sedq_now(i) + vh_sedq_new(i))

           ! qflux_new, nflux_new are the updated flux values from the level above
           ! qflux_now, nflux_now are here the old (current time step) flux values from the level above
           ! In COSMO-Docu  {...} =  flux_(k-1),new + flux_(k-1),start
           nh_flux_sum(i) = nh_flux_new(i) + nh_flux_now(i)
           qh_flux_sum(i) = qh_flux_new(i) + qh_flux_now(i)

           ! qflux_now, nflux_now are here overwritten with the current level
           nh_flux_now(i) = min(vh_sedn_now(i) * qnh(i,k), nh_flux_sum(i))  ! this is then passed to the level below
           qh_flux_now(i) = min(vh_sedq_now(i) * qh(i,k),  qh_flux_sum(i))  ! (loop dependency)
           nh_flux_now(i) = max(nh_flux_now(i),0.0_wp)
           qh_flux_now(i) = max(qh_flux_now(i),0.0_wp)

           nh_sum(i) = qnh(i,k) + rdzdt(i,k) * (nh_flux_sum(i) - nh_flux_now(i))
           qh_sum(i) = qh(i,k)  + rdzdt(i,k) * (qh_flux_sum(i) - qh_flux_now(i))

           nh_impl(i) = 1.0_wp/(1.0_wp + vh_sedn_new(i) * rdzdt(i,k))
           qh_impl(i) = 1.0_wp/(1.0_wp + vh_sedq_new(i) * rdzdt(i,k))

           nh_star(i) = nh_impl(i) * nh_sum(i)  ! intermediate values for calculating
           qh_star(i) = qh_impl(i) * qh_sum(i)  ! sources and sinks

           qnh(i,k) = nh_star(i)                ! overwrite array with intermediate
           qh(i,k)  = qh_star(i)                ! values to do micro processes on this level

           nh_sum(i) = nh_sum(i) - nh_star(i)   ! final time integration starts from sum-values
           qh_sum(i) = qh_sum(i) - qh_star(i)   ! but source/sinks work on star-values

        end do

        ! do microphysics on this k-level only (using the star-values)
        IF (lmicro_impl) THEN

           ! .. save old variables for latent heat calculation
           DO ii = its, ite
              q_vap_old(ii,k) = qv(ii,k)
              q_liq_old(ii,k) = qc(ii,k) + qr(ii,k)
           END DO

           ik_slice(3) = k
           ik_slice(4) = k
           CALL clouds_twomoment(ik_slice, dt, lprogin, &
                atmo, cloud, rain, ice, snow, graupel, hail, &
                ninact, nccn, ninpot)

           ! .. latent heat term for temperature equation
           DO ii = its, ite
              q_vap_new  = qv(ii,k)
              q_liq_new  = qr(ii,k) + qc(ii,k)
              tk(ii,k)   = tk(ii,k) - convice * rho_r(ii,k) * (q_vap_new - q_vap_old(ii,k))  &
                    &               + convliq * rho_r(ii,k) * (q_liq_new - q_liq_old(ii,k))
           ENDDO

        END IF

        do i=its,ite

           ! time integration
           qnr(i,k) = max( 0.0_wp, nr_impl(i)*(nr_sum(i) + qnr(i,k)))
           qni(i,k) = max( 0.0_wp, ni_impl(i)*(ni_sum(i) + qni(i,k)))
           qns(i,k) = max( 0.0_wp, ns_impl(i)*(ns_sum(i) + qns(i,k)))
           qng(i,k) = max( 0.0_wp, ng_impl(i)*(ng_sum(i) + qng(i,k)))
           qnh(i,k) = max( 0.0_wp, nh_impl(i)*(nh_sum(i) + qnh(i,k)))
           qr(i,k)  = max( 0.0_wp, qr_impl(i)*(qr_sum(i) + qr(i,k)))
           qi(i,k)  = max( 0.0_wp, qi_impl(i)*(qi_sum(i) + qi(i,k)))
           qs(i,k)  = max( 0.0_wp, qs_impl(i)*(qs_sum(i) + qs(i,k)))
           qg(i,k)  = max( 0.0_wp, qg_impl(i)*(qg_sum(i) + qg(i,k)))
           qh(i,k)  = max( 0.0_wp, qh_impl(i)*(qh_sum(i) + qh(i,k)))

           ! prepare for next level
           nr_flux_new(i) = qnr(i,k) * vr_sedn_new(i)     ! flux_(k),new
           ni_flux_new(i) = qni(i,k) * vi_sedn_new(i)     ! for next level (loop dependency)
           ns_flux_new(i) = qns(i,k) * vs_sedn_new(i)     !
           ng_flux_new(i) = qng(i,k) * vg_sedn_new(i)     !
           nh_flux_new(i) = qnh(i,k) * vh_sedn_new(i)     !
           qr_flux_new(i) = qr(i,k)  * vr_sedq_new(i)     !
           qi_flux_new(i) = qi(i,k)  * vi_sedq_new(i)     !
           qs_flux_new(i) = qs(i,k)  * vs_sedq_new(i)     !
           qg_flux_new(i) = qg(i,k)  * vg_sedq_new(i)     !
           qh_flux_new(i) = qh(i,k)  * vh_sedq_new(i)     !

           vr_sedn_new(i) = vr_sedn_now(i)
           vi_sedn_new(i) = vi_sedn_now(i)
           vs_sedn_new(i) = vs_sedn_now(i)
           vg_sedn_new(i) = vg_sedn_now(i)
           vh_sedn_new(i) = vh_sedn_now(i)
           vr_sedq_new(i) = vr_sedq_now(i)
           vi_sedq_new(i) = vi_sedq_now(i)
           vs_sedq_new(i) = vs_sedq_now(i)
           vg_sedq_new(i) = vg_sedq_now(i)
           vh_sedq_new(i) = vh_sedq_now(i)

        end do

     END DO

   END SUBROUTINE clouds_twomoment_implicit
   
   !
   ! sedimentation for explicit solver, i.e., sedimentation is done with an explicit
   ! flux-form semi-lagrangian scheme after the microphysics.
   !
   SUBROUTINE sedimentation_explicit()
     ! D.Rieger: the parameter lfullyexplicit needs to be set false, otherwise the nproma/mpi tests of buildbot are not passed
     logical, parameter :: lfullyexplicit = .false.
     REAL(wp) :: cmax, rdzmaxdt

     cmax = 0.0_wp

     IF (lfullyexplicit) THEN
       rdzmaxdt = maxval(rdz(its:ite,kts:kte)) * dt
       ntsedi = ceiling(rain%vsedi_max*rdzmaxdt)
     ELSE
       ntsedi = 1       
     END IF

     prec_r(:) = 0.0_wp
     prec_i(:) = 0.0_wp
     prec_s(:) = 0.0_wp
     prec_g(:) = 0.0_wp
     prec_h(:) = 0.0_wp

     DO ii=1,ntsedi       
       CALL sedi_icon_rain (rain,rain_coeffs,qr,qnr,prec_r,qc,rhocorr,rdz,dt/ntsedi,its,ite,kts,kte,cmax)
     END DO
          
     IF (ANY(qi(its:ite,kts:kte)>0._wp)) &
          call sedi_icon_sphere (ice,ice_coeffs,qi,qni,prec_i,rhocorr,rdz,dt,its,ite,kts,kte)
     
     IF (ANY(qs(its:ite,kts:kte)>0._wp)) &
          call sedi_icon_sphere (snow,snow_coeffs,qs,qns,prec_s,rhocorr,rdz,dt,its,ite,kts,kte)
     
     IF (ANY(qg(its:ite,kts:kte)>0._wp)) THEN
       if (lfullyexplicit) then
         ntsedi = ceiling(graupel%vsedi_max*rdzmaxdt)
       else
         ntsedi = 1
       end if       
       IF (lprogmelt) THEN
         DO ii=1,ntsedi
           call sedi_icon_sphere_lwf(graupel_lwf,graupel_coeffs,qg,qng,qgl,&
                &                    prec_g,rhocorr,rdz,dt/ntsedi,its,ite,kts,kte,cmax)
         END DO
       ELSE
         DO ii=1,ntsedi
           call sedi_icon_sphere (graupel,graupel_coeffs,qg,qng,prec_g,rhocorr,rdz,dt/ntsedi,its,ite,kts,kte,cmax)
         END DO
       END IF
     END IF

     IF (ANY(qh(its:ite,kts:kte)>0._wp)) THEN
       if (lfullyexplicit) then
         ntsedi = ceiling(hail%vsedi_max*rdzmaxdt)
       else
         ntsedi = 1
       end if       
       IF (lprogmelt) THEN
         DO ii=1,ntsedi
           call sedi_icon_sphere_lwf(hail_lwf,hail_coeffs,qh,qnh,qhl,&
                &                    prec_h,rhocorr,rdz,dt/ntsedi,its,ite,kts,kte,cmax)
         END DO
       ELSE
         DO ii=1,ntsedi
           call sedi_icon_sphere (hail,hail_coeffs,qh,qnh,prec_h,rhocorr,rdz,dt/ntsedi,its,ite,kts,kte,cmax)
         END DO
       END IF
     END IF
     
     IF (msg_level > 100)THEN
       WRITE (message_text,'(1X,A,f8.2)') ' sedimentation_explicit  cmax = ',cmax
       CALL message(routine,TRIM(message_text))
     END IF

   END SUBROUTINE sedimentation_explicit

    !
    ! check for negative values after microphysics
    !
    SUBROUTINE check_clouds()

      REAL(wp), PARAMETER :: meps = -1e-12

      IF (cloud_type.lt.2000) THEN
         IF (ANY(qh(its:ite,kts:kte)>0._wp)) THEN
            qh(its:ite,kts:kte)  = 0.0_wp
            WRITE (message_text,'(1X,A)') '  qh > 0, after cloud_twomoment for cloud_type < 2000'
            CALL message(routine,TRIM(message_text))
            CALL finish(TRIM(routine),'Error in two_moment_mcrph')
         END IF
         IF (ANY(qnh(its:ite,kts:kte)>0._wp)) THEN
            qnh(its:ite,kts:kte)  = 0.0_wp
            WRITE (message_text,'(1X,A)') '  qnh > 0, after cloud_twomoment for cloud_type < 2000'
            CALL message(routine,TRIM(message_text))
            CALL finish(TRIM(routine),'Error in two_moment_mcrph')
         END IF
      END IF
      IF (msg_level>dbg_level) CALL message(TRIM(routine), " test for negative values")
      IF (MINVAL(cloud%q(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, cloud%q < 0')
      ENDIF
      IF (MINVAL(rain%q(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, rain%q < 0')
      ENDIF
      IF (MINVAL(ice%q(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, ice%q < 0,')
      ENDIF
      IF (MINVAL(snow%q(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, snow%q < 0')
      ENDIF
      IF (MINVAL(graupel%q(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, graupel%q < 0')
      ENDIF
      IF (MINVAL(hail%q(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, hail%q < 0')
      ENDIF
!      IF (MINVAL(cloud%n) < meps) THEN
!         CALL finish(TRIM(routine),'Error in two_moment_mcrph, cloud%n < 0')
!      ENDIF
      IF (MINVAL(rain%n(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, rain%n < 0')
      ENDIF
      IF (MINVAL(ice%n(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, ice%n < 0')
      ENDIF
      IF (MINVAL(snow%n(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, snow%n < 0')
      ENDIF
      IF (MINVAL(graupel%n(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, graupel%n < 0')
      ENDIF
      IF (MINVAL(hail%n(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, hail%n < 0')
      ENDIF
    END subroutine check_clouds

  END SUBROUTINE two_moment_mcrph

!===========================================================================================

  SUBROUTINE two_moment_mcrph_init(igscp,N_cn0,z0_nccn,z1e_nccn,N_in0,z0_nin,z1e_nin,msg_level)

    INTEGER, INTENT(IN) :: igscp, msg_level

    REAL(wp), OPTIONAL, INTENT(OUT) ::             & ! for CCN and IN in case of gscp=5
         & N_cn0,z0_nccn,z1e_nccn,    &
         & N_in0,z0_nin,z1e_nin

    TYPE(particle)        :: cloud, rain
    TYPE(particle_frozen) :: ice, snow, graupel, hail
    TYPE(particle_lwf)    :: graupel_lwf, hail_lwf

    INTEGER        :: unitnr

    IF (msg_level>5) CALL message (TRIM(routine), " Initialization of two-moment microphysics scheme")

    unitnr = 11
    CALL init_dmin_wetgrowth('dmin_wetgrowth_lookup.dat', unitnr)

    CALL init_dmin_wg_gr_ltab_equi('dmin_wetgrowth_lookup.dat', unitnr, 61, ltabdminwgg)

    IF (msg_level>dbg_level) CALL message (TRIM(routine), " finished init_dmin_wetgrowth")

    IF (PRESENT(N_cn0)) THEN
       ccn_type   = ccn_type_gscp5
       cloud_type = cloud_type_default_gscp5 + 10 * ccn_type
    ELSE
       ccn_type   = ccn_type_gscp4
       cloud_type = cloud_type_default_gscp4 + 10 * ccn_type
    END IF

    ! .. set the particle types, and calculate some coefficients
    IF (igscp.eq.7) THEN
       CALL init_2mom_scheme_once(cloud,rain,ice,snow,graupel_lwf,hail_lwf,cloud_type)
    ELSE
       CALL init_2mom_scheme_once(cloud,rain,ice,snow,graupel,hail,cloud_type)
    END IF

    IF (present(N_cn0)) THEN

       !..parameters for CCN and IN are set here. The 3D fields are then
       !  initialized in mo_nwp_phy_init.

       !..parameters for exponential decrease of N_ccn with height
       !  z0:  up to this height (m) constant unchanged value
       !  z1e: height interval at which N_ccn decreases by factor 1/e above z0_nccn

       ccn_coeffs%z0  = 4000.0d0
       ccn_coeffs%z1e = 2000.0d0

       ! characteristics of different kinds of CN
       ! (copied from COSMO 5.0 Segal & Khain nucleation subroutine)
       SELECT CASE(ccn_type)
       CASE(6)
          !... maritime case
          ccn_coeffs%Ncn0 = 100.0d06   ! CN concentration at ground
          ccn_coeffs%Nmin =  35.0d06
          ccn_coeffs%lsigs = 0.4d0      ! log(sigma_s)
          ccn_coeffs%R2    = 0.03d0     ! in mum
          ccn_coeffs%etas  = 0.9        ! soluble fraction
       CASE(7)
          !... intermediate case
          ccn_coeffs%Ncn0 = 500.0d06
          ccn_coeffs%Nmin =  35.0d06
          ccn_coeffs%lsigs = 0.4d0
          ccn_coeffs%R2    = 0.03d0       ! in mum
          ccn_coeffs%etas  = 0.8          ! soluble fraction
       CASE(8)
          !... continental case
          ccn_coeffs%Ncn0 = 1700.0d06
          ccn_coeffs%Nmin =   35.0d06
          ccn_coeffs%lsigs = 0.2d0
          ccn_coeffs%R2    = 0.03d0       ! in mum
          ccn_coeffs%etas  = 0.7          ! soluble fraction
       CASE(9)
          !... "polluted" continental
          ccn_coeffs%Ncn0 = 3200.0d06
          ccn_coeffs%Nmin =   35.0d06
          ccn_coeffs%lsigs = 0.2d0
          ccn_coeffs%R2    = 0.03d0       ! in mum
          ccn_coeffs%etas  = 0.7          ! soluble fraction
       CASE(1)
          !... dummy values
          ccn_coeffs%Ncn0 =  200.0d06
          ccn_coeffs%Nmin =   10.0d06
          ccn_coeffs%lsigs = 0.0
          ccn_coeffs%R2    = 0.0
          ccn_coeffs%etas  = 0.0
       CASE DEFAULT
          CALL finish(TRIM(routine),'Error in two_moment_mcrph_init: Invalid value for ccn_type')
       END SELECT

       z0_nccn  = ccn_coeffs%z0
       z1e_nccn = ccn_coeffs%z1e
       N_cn0    = ccn_coeffs%Ncn0

       WRITE(message_text,'(A)') "  CN properties:" ; CALL message(TRIM(routine),TRIM(message_text))
       WRITE(message_text,'(A,D10.3)') "    Ncn0 = ",ccn_coeffs%Ncn0 ; CALL message(TRIM(routine),TRIM(message_text))
       WRITE(message_text,'(A,D10.3)') "    z0   = ",ccn_coeffs%z0  ; CALL message(TRIM(routine),TRIM(message_text))
       WRITE(message_text,'(A,D10.3)') "    z1e  = ",ccn_coeffs%z1e ; CALL message(TRIM(routine),TRIM(message_text))
    END IF

    IF (present(N_in0)) THEN

       in_coeffs%N0  = 200.0e6_wp ! this is currently just a scaling factor for the PDA scheme
       in_coeffs%z0  = 3000.0d0
       in_coeffs%z1e = 1000.0d0

       N_in0   = in_coeffs%N0
       z0_nin  = in_coeffs%z0
       z1e_nin = in_coeffs%z1e

       WRITE(message_text,'(A)') "  IN properties:" ; CALL message(TRIM(routine),TRIM(message_text))
       WRITE(message_text,'(A,D10.3)') "    Ncn0 = ",in_coeffs%N0  ; CALL message(TRIM(routine),TRIM(message_text))
       WRITE(message_text,'(A,D10.3)') "    z0   = ",in_coeffs%z0  ; CALL message(TRIM(routine),TRIM(message_text))
       WRITE(message_text,'(A,D10.3)') "    z1e  = ",in_coeffs%z1e ; CALL message(TRIM(routine),TRIM(message_text))
    END IF

  END SUBROUTINE two_moment_mcrph_init

!==============================================================================

  FUNCTION set_qnc(qc)

    REAL(wp), INTENT(in)  :: qc
    REAL(wp) :: set_qnc
    REAL(wp), PARAMETER   :: Dmean = 10e-6_wp    ! Diameter of mean particle mass:

    set_qnc = qc * 6.0_wp / (pi * rhoh2o * Dmean**3.0_wp)

  END FUNCTION set_qnc

  FUNCTION set_qni(qi)
    REAL(wp) :: set_qni
    REAL(wp), INTENT(in)  :: qi
    !  REAL(wp), PARAMETER   :: Dmean = 100e-6_wp  ! Diameter of mean particle mass:

    set_qni  = qi / 1e-10   !  qiin / ( ( Dmean / ageo) ** (1.0_wp / bgeo) )

  END FUNCTION set_qni

  FUNCTION set_qnr(qr)
    REAL(wp) :: set_qnr
    REAL(wp), INTENT(in)  :: qr
    REAL(wp), PARAMETER   :: N0r = 8000.0e3_wp ! intercept of MP distribution

    set_qnr = N0r * ( qr * 6.0_wp / (pi * rhoh2o * N0r * gfct(4.0_wp)))**(0.25_wp)

  END FUNCTION set_qnr

  FUNCTION set_qns(qs)
    REAL(wp) :: set_qns
    REAL(wp), INTENT(in)  :: qs
    REAL(wp), PARAMETER :: N0s = 800.0e3_wp

    REAL(wp), PARAMETER :: ams = 0.038_wp  ! needs to be connected to snow-type
    REAL(wp), PARAMETER :: bms = 2.0_wp

    set_qns = N0s * ( qs / ( ams * N0s * gfct(bms+1.0_wp)))**( 1.0_wp/(1.0_wp+bms) )

  END FUNCTION set_qns

  FUNCTION set_qng(qg)
    REAL(wp) :: set_qng
    REAL(wp), INTENT(in)  :: qg

    REAL(wp), PARAMETER   :: N0g = 4000.0e3_wp
    REAL(wp), PARAMETER   :: amg = 169.6_wp     ! needs to be connected to graupel-type
    REAL(wp), PARAMETER   :: bmg = 3.1_wp

    set_qng = N0g * ( qg / ( amg * N0g * gfct(bmg+1.0_wp)))**( 1.0_wp/(1.0_wp+bmg) )

  END FUNCTION set_qng

END MODULE mo_mcrph_sb
