!>
!! cloud microphysics
!!
!!==============================================================================
!!
!! Two-moment mixed-phase bulk microphysics
!!
!! original version by Axel Seifert, May 2003
!! with modifications by Ulrich Blahak, August 2007
!!
!!==============================================================================
!!
!! $Id: n/a$
!!
!! @par Revision History
!! Ported into ICON from UCLA-LES by Anurag Dipankar (2013-12-15) 
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

USE mo_kind              , ONLY: wp 
USE mo_physical_constants, ONLY: &
    alv,           & !! latent heat of vaporization
    als,           & !! latent heat of sublimation
    cpdr  => rcpd, & !! (spec. heat of dry air at constant press)^-1
    cvdr  => rcvd    !! (spec. heat of dry air at const vol)^-1

USE mo_exception,      ONLY: finish, message, message_text

#ifndef __SX__
USE mo_2mom_mcrph_main,     ONLY:                              &
     &                       istart,iend,kstart,kend,          &
     &                       clouds_twomoment,                 &
     &                       sedi_vel_rain, sedi_vel_sphere,   &
     &                       sedi_icon_rain,sedi_icon_sphere,  &
     &                       dt_twomoment => dt,               &
     &                       w_p, T_p, p_p,                    &
     &                       ptr_rho => rho_p, ptr_qv => qv,            &
     &                       rrho_04, rrho_c, rho_vel, rho_vel_c, rho0, &
     &                       q_cloud, q_ice, q_rain, q_snow, q_graupel, &
     &                       n_cloud, n_ice, n_rain, n_snow, n_graupel, &
     &                       cloud, rain, ice, snow, graupel, hail,     &
     &                       rain_coeffs, ice_coeffs, snow_coeffs, graupel_coeffs, hail_coeffs, &
     &                       n_hail, q_hail,                            &
     &                       ltabdminwgg,                               &
     &                       init_2mom_scheme,                          &
     &                       qnc_const
USE mo_2mom_mcrph_util,     ONLY:                              &
     &                       init_dmin_wetgrowth,              &
     &                       init_dmin_wg_gr_ltab_equi
#endif

!==============================================================================

IMPLICIT NONE
PUBLIC 

INTEGER :: cloud_type = 2503
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
                       qc,                & ! inout: cloud water
                       qr, qnr,           & ! inout: rain
                       qi, qni,           & ! inout: ice
                       qs, qns,           & ! inout: snow
                       qg, qng,           & ! inout: graupel
                       qh, qnh,           & ! inout: hail
                       tk,                & ! inout: temp 
                       w,                 & ! inout: w
                       prec_r,            & ! inout: precip rate rain
                       prec_i,            & ! inout: precip rate ice
                       prec_s,            & ! inout: precip rate snow
                       prec_g,            & ! inout: precip rate graupel
                       prec_h,            & ! inout: precip rate hail
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
         &               qv, qc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh

    ! Precip rates, vertical profiles
    REAL(wp), DIMENSION(:), INTENT (INOUT)  :: &
         &               prec_r, prec_i, prec_s, prec_g, prec_h

    INTEGER,  INTENT (IN)             :: msg_level 
    LOGICAL,  OPTIONAL,  INTENT (IN)  :: l_cv

    ! ... Variables which are global in module_2mom_mcrph_main

    REAL(wp), TARGET, DIMENSION(isize,ke) ::        &
         &  qnc,           & ! cloud droplet number
         &  rhocorr,       & ! density dependency of particle fall speed
         &  rhocld           ! density dependency of particle fall speed for cloud droplets

    REAL(wp) :: q_liq_old(isize,ke), q_vap_old(isize,ke)  ! to store old values for latent heat calc

    INTEGER  :: its,ite,kts,kte
    INTEGER  :: ii,kk
    INTEGER  :: ntsedi     ! for sedimentation sub stepping
    
    REAL(wp) :: q_liq_new,q_vap_new
    REAL(wp) :: hlp
    REAL(wp) :: convliq,convice
    REAL(wp) :: z_heat_cap_r !! reciprocal of cpdr or cvdr (depending on l_cv)
    REAL(wp) :: rdz(isize,ke), rho_r(isize,ke)

    INTEGER, PARAMETER :: dbg_level = 25            ! level for debug prints
    LOGICAL, PARAMETER :: debug = .false.    
    CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_mcrph_driver'

    LOGICAL, PARAMETER :: explicit_solver = .true.  ! explicit or semi-implicit solver

#ifndef __SX__

    ! inverse of vertical layer thickness
    rdz = 1._wp / dz

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

    ! indices as used in two-moment scheme
    istart = its
    iend   = ite
    kstart = kts
    kend   = kte
     
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

    ! ... constant cloud droplet number
    qnc_const = 200.0e6_wp

    IF (msg_level>dbg_level)THEN
       WRITE (message_text,'(1X,A,I4,A,E12.4)') "cloud_type = ",cloud_type,", qnc_const = ",qnc_const
       CALL message(TRIM(routine),TRIM(message_text))
    END IF

    ! time step for two-moment microphysics is the same as all other fast physics
    dt_twomoment = dt

    IF (msg_level>dbg_level) CALL message(TRIM(routine), "prepare variable for 2mom")

    DO kk = kts, kte
       DO ii = its, ite

          rho_r(ii,kk) = 1.0 / rho(ii,kk)

          ! ... save old variables for latent heat calculation
          q_vap_old(ii,kk) = qv(ii,kk)
          q_liq_old(ii,kk) = qc(ii,kk) + qr(ii,kk)

          ! ... height dependency of terminal fall velocities
          hlp = log(rho(ii,kk)/rho0)
          rhocorr(ii,kk) = exp(-rho_vel*hlp)
          rhocld(ii,kk)  = exp(-rho_vel_c*hlp)

       END DO
    END DO

    ! .. convert to densities and set pointerns to two-moment module
    !    (pointers are used to avoid passing everything explicitly by argument and
    !     to avoid local allocates within the OpenMP-loop, and keep everything on stack)

    CALL prepare_twomoment()
                     
    IF (msg_level>dbg_level) CALL message(TRIM(routine)," calling clouds_twomoment")

    IF (explicit_solver) then
       ! .. this subroutine calculates all the microphysical sources and sinks
       CALL clouds_twomoment (isize,ke)
    ELSE
       ! .. semi-implicit solver doing microphysics and sedimentation
       
       ! TO-DO: Latent heat with implicit solver is wrong
       !        Budget has to be done on each k-level

       CALL clouds_twomoment_implicit ()
    END IF

    ! .. check for negative values
    IF (debug) CALL check_clouds()

    ! .. latent heat term for temperature equation
    convice = z_heat_cap_r * als
    convliq = z_heat_cap_r * (alv-als)
    DO kk = kts, kte
       DO ii = its, ite

          ! .. new variables
          q_vap_new = qv(ii,kk)
          q_liq_new = qr(ii,kk) + qc(ii,kk)

          ! .. update temperature
          tk(ii,kk) = tk(ii,kk) - convice * (rho_r(ii,kk)*q_vap_new - q_vap_old(ii,kk))  &
               &                + convliq * (rho_r(ii,kk)*q_liq_new - q_liq_old(ii,kk))
       ENDDO
    ENDDO

    IF (explicit_solver) THEN
       IF (msg_level>dbg_level) CALL message(TRIM(routine),"sedimentation")

       ! .. if we solve explicitly, then sedimentation is done here after microphysics
       CALL sedimentation_explicit()

    END IF

    ! .. convert back and nullify two-moment pointers
    CALL post_twomoment()

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
      real(wp), dimension(isize,ke) :: rdzdt 

      integer :: i,k

      logical, parameter :: lmicro_impl = .false.  ! microphysics within semi-implicit sedimentation loop?

      if (.not.lmicro_impl) then
         CALL clouds_twomoment (isize,ke)
      end if

      rdzdt = 0.5_wp * rdz * dt

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

      ! here we simply assume that there is no cloud or precip in the uppermost level
      ! i.e. we start from kts+1 going down in physical space
     
!DIR$ IVDEP 
      DO k=kts+1,kte

        do i=its,ite
           xr_now(i) = rain%meanmass(qr(i,k),qnr(i,k))                       
           xr_new(i) = 0.5*(rain%meanmass(qr(i,k-1),qnr(i,k-1)) + xr_now(i))

           xi_now(i) = ice%meanmass(qi(i,k),qni(i,k))                       
           xi_new(i) = 0.5*(ice%meanmass(qi(i,k-1),qni(i,k-1)) + xi_now(i))

           xs_now(i) = ice%meanmass(qs(i,k),qns(i,k))                       
           xs_new(i) = 0.5*(ice%meanmass(qs(i,k-1),qns(i,k-1)) + xs_now(i))

           xg_now(i) = ice%meanmass(qg(i,k),qng(i,k))                       
           xg_new(i) = 0.5*(ice%meanmass(qg(i,k-1),qng(i,k-1)) + xg_now(i))

           xh_now(i) = ice%meanmass(qh(i,k),qnh(i,k))                       
           xh_new(i) = 0.5*(ice%meanmass(qh(i,k-1),qnh(i,k-1)) + xh_now(i))
        end do

        call sedi_vel_rain(rain,rain_coeffs,xr_new,vr_sedn_new,vr_sedq_new,its,ite,qc(:,k))
        call sedi_vel_rain(rain,rain_coeffs,xr_now,vr_sedn_now,vr_sedq_now,its,ite,qc(:,k))

        call sedi_vel_sphere(ice,ice_coeffs,xi_new,vi_sedn_new,vi_sedq_new,its,ite)
        call sedi_vel_sphere(ice,ice_coeffs,xi_now,vi_sedn_now,vi_sedq_now,its,ite)

        call sedi_vel_sphere(snow,snow_coeffs,xs_new,vs_sedn_new,vs_sedq_new,its,ite)
        call sedi_vel_sphere(snow,snow_coeffs,xs_now,vs_sedn_now,vs_sedq_now,its,ite)

        call sedi_vel_sphere(graupel,graupel_coeffs,xg_new,vg_sedn_new,vg_sedq_new,its,ite)
        call sedi_vel_sphere(graupel,graupel_coeffs,xg_now,vg_sedn_now,vg_sedq_now,its,ite)

        call sedi_vel_sphere(hail,hail_coeffs,xh_new,vh_sedn_new,vh_sedq_new,its,ite)
        call sedi_vel_sphere(hail,hail_coeffs,xh_now,vh_sedn_now,vh_sedq_now,its,ite)

        do i=its,ite

           ! .... rain ....

           ! qflux_new, nflux_new are the updated flux values from the level above
           ! qflux_now, nflux_now are here the old (current time step) flux values from the level above
           ! In COSMO-Docu  {...} =  flux_(k-1),new + flux_(k-1),start
           nr_flux_sum(i) = nr_flux_new(i) + nr_flux_now(i)
           qr_flux_sum(i) = qr_flux_new(i) + qr_flux_now(i)

           ! qflux_now, nflux_now are here overwritten with the current level
           nr_flux_now(i) = min(vr_sedn_now(i) * qnr(i,k), nr_flux_sum(i))  ! this is then passed to the level below
           qr_flux_now(i) = min(vr_sedq_now(i) * qr(i,k),  qr_flux_sum(i))  ! (loop dependency)

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

           ! qflux_new, nflux_new are the updated flux values from the level above
           ! qflux_now, nflux_now are here the old (current time step) flux values from the level above
           ! In COSMO-Docu  {...} =  flux_(k-1),new + flux_(k-1),start
           ni_flux_sum(i) = ni_flux_new(i) + ni_flux_now(i)
           qi_flux_sum(i) = qi_flux_new(i) + qi_flux_now(i)

           ! qflux_now, nflux_now are here overwritten with the current level
           ni_flux_now(i) = min(vi_sedn_now(i) * qni(i,k), ni_flux_sum(i))  ! this is then passed to the level below
           qi_flux_now(i) = min(vi_sedq_now(i) * qi(i,k),  qi_flux_sum(i))  ! (loop dependency)

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

           ! qflux_new, nflux_new are the updated flux values from the level above
           ! qflux_now, nflux_now are here the old (current time step) flux values from the level above
           ! In COSMO-Docu  {...} =  flux_(k-1),new + flux_(k-1),start
           ns_flux_sum(i) = ns_flux_new(i) + ns_flux_now(i)
           qs_flux_sum(i) = qs_flux_new(i) + qs_flux_now(i)

           ! qflux_now, nflux_now are here overwritten with the current level
           ns_flux_now(i) = min(vs_sedn_now(i) * qns(i,k), ns_flux_sum(i))  ! this is then passed to the level below
           qs_flux_now(i) = min(vs_sedq_now(i) * qs(i,k),  qs_flux_sum(i))  ! (loop dependency)

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

           ! qflux_new, nflux_new are the updated flux values from the level above
           ! qflux_now, nflux_now are here the old (current time step) flux values from the level above
           ! In COSMO-Docu  {...} =  flux_(k-1),new + flux_(k-1),start
           ng_flux_sum(i) = ng_flux_new(i) + ng_flux_now(i)
           qg_flux_sum(i) = qg_flux_new(i) + qg_flux_now(i)

           ! qflux_now, nflux_now are here overwritten with the current level
           ng_flux_now(i) = min(vg_sedn_now(i) * qng(i,k), ng_flux_sum(i))  ! this is then passed to the level below
           qg_flux_now(i) = min(vg_sedq_now(i) * qg(i,k),  qg_flux_sum(i))  ! (loop dependency)

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

           ! qflux_new, nflux_new are the updated flux values from the level above
           ! qflux_now, nflux_now are here the old (current time step) flux values from the level above
           ! In COSMO-Docu  {...} =  flux_(k-1),new + flux_(k-1),start
           nh_flux_sum(i) = nh_flux_new(i) + nh_flux_now(i)
           qh_flux_sum(i) = qh_flux_new(i) + qh_flux_now(i)

           ! qflux_now, nflux_now are here overwritten with the current level
           nh_flux_now(i) = min(vh_sedn_now(i) * qnh(i,k), nh_flux_sum(i))  ! this is then passed to the level below
           qh_flux_now(i) = min(vh_sedq_now(i) * qh(i,k),  qh_flux_sum(i))  ! (loop dependency)

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

        if (lmicro_impl) then
           ! do microphysics on this k-level only (using the star-values)
           kstart = k  
           kend   = k
           call clouds_twomoment(isize,ke)
        end if

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
        end do

     END DO

   END SUBROUTINE clouds_twomoment_implicit

   SUBROUTINE prepare_twomoment()
     
     ! ... Transformation of microphysics variables to densities
     DO kk = kts, kte
        DO ii = its, ite
           
           ! ... concentrations --> number densities
           qnr(ii,kk) = rho(ii,kk) * qnr(ii,kk) 
           qni(ii,kk) = rho(ii,kk) * qni(ii,kk) 
           qns(ii,kk) = rho(ii,kk) * qns(ii,kk)
           qng(ii,kk) = rho(ii,kk) * qng(ii,kk)
           qnh(ii,kk) = rho(ii,kk) * qnh(ii,kk)
           
           ! ... mixing ratios -> mass densities
           qv(ii,kk) = rho(ii,kk) * qv(ii,kk) 
           qc(ii,kk) = rho(ii,kk) * qc(ii,kk) 
           qr(ii,kk) = rho(ii,kk) * qr(ii,kk) 
           qi(ii,kk) = rho(ii,kk) * qi(ii,kk) 
           qs(ii,kk) = rho(ii,kk) * qs(ii,kk) 
           qg(ii,kk) = rho(ii,kk) * qg(ii,kk) 
           qh(ii,kk) = rho(ii,kk) * qh(ii,kk) 
           
        END DO
     END DO
     
     ! set pointers
     w_p => w
     t_p => tk
     p_p => pres
     ptr_qv => qv
     ptr_rho => rho
     
     rrho_04 => rhocorr
     rrho_c  => rhocld

     q_cloud   => qc
     n_cloud   => qnc
     q_rain    => qr
     n_rain    => qnr
     q_ice     => qi
     n_ice     => qni
     q_snow    => qs
     n_snow    => qns
     q_graupel => qg
     n_graupel => qng
     q_hail    => qh
     n_hail    => qnh
     
   END SUBROUTINE prepare_twomoment

   SUBROUTINE post_twomoment()

     ! nullify pointers
     w_p => null()
     t_p => null()
     p_p => null()
     ptr_qv  => null()
     ptr_rho => null()
     
     rrho_04 => null()
     rrho_c  => null()
     
     q_cloud   => null()
     n_cloud   => null()
     q_rain    => null()
     n_rain    => null()
     q_ice     => null()
     n_ice     => null()
     q_snow    => null()
     n_snow    => null()
     q_graupel => null()
     n_graupel => null()
     q_hail    => null()
     n_hail    => null()

     ! ... Transformation of variables back to ICON standard variables
     DO kk = kts, kte
        DO ii = its, ite
           
           hlp = rho_r(ii,kk)

           ! ... from mass densities back to mixing ratios
           qv(ii,kk) = hlp * qv(ii,kk)
           qc(ii,kk) = hlp * qc(ii,kk)
           qr(ii,kk) = hlp * qr(ii,kk)
           qi(ii,kk) = hlp * qi(ii,kk)
           qs(ii,kk) = hlp * qs(ii,kk)
           qg(ii,kk) = hlp * qg(ii,kk)
           qh(ii,kk) = hlp * qh(ii,kk)

           ! ... number concentrations
           qnr(ii,kk) = hlp * qnr(ii,kk)
           qni(ii,kk) = hlp * qni(ii,kk)
           qns(ii,kk) = hlp * qns(ii,kk)
           qng(ii,kk) = hlp * qng(ii,kk)
           qnh(ii,kk) = hlp * qnh(ii,kk)

        ENDDO
     ENDDO

   END SUBROUTINE post_twomoment

   !
   ! sedimentation for explicit solver, i.e., sedimentation is done with an explicit
   ! flux-form semi-lagrangian scheme after the microphysics.
   !
   SUBROUTINE sedimentation_explicit()

     REAL(wp) :: cmax
     
     cmax = 0.0_wp

     prec_r  = 0._wp
     CALL sedi_icon_rain (qr,qnr,prec_r,qc,rhocorr,rdz,dt,its,ite,kts,kte,cmax)
      
      IF (cloud_type.ge.1000) THEN
         prec_i(:) = 0._wp
         prec_s(:) = 0._wp
         prec_g(:) = 0._wp
         
         IF (ANY(qi(:,:)>0._wp)) &
              call sedi_icon_sphere (ice,ice_coeffs,qi,qni,prec_i,rhocorr,rdz,dt,its,ite,kts,kte)
         
         IF (ANY(qs(:,:)>0._wp)) &
              call sedi_icon_sphere (snow,snow_coeffs,qs,qns,prec_s,rhocorr,rdz,dt,its,ite,kts,kte)
         
         IF (ANY(qg(:,:)>0._wp)) THEN
            ntsedi = 1
            DO ii=1,ntsedi
               call sedi_icon_sphere (graupel,graupel_coeffs,qg,qng,prec_g,rhocorr,rdz,dt/ntsedi,its,ite,kts,kte,cmax)
            END DO
         END IF
      END IF

      IF (cloud_type.ge.2000) THEN
         prec_h(:) = 0.0
         IF (ANY(qh(:,:)>0._wp)) THEN
            ntsedi = 1
            DO ii=1,ntsedi
               call sedi_icon_sphere (hail,hail_coeffs,qh,qnh,prec_h,rhocorr,rdz,dt/ntsedi,its,ite,kts,kte,cmax)
            END DO
         ENDIF
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

      IF (cloud_type.lt.2000) THEN
         IF (ANY(qh(:,:)>0._wp)) THEN
            qh(:,:)  = 0.0_wp
            WRITE (message_text,'(1X,A)') '  qh > 0, after cloud_twomoment for cloud_type < 2000'
            CALL message(routine,TRIM(message_text))
            CALL finish(TRIM(routine),'Error in two_moment_mcrph')
         END IF
         IF (ANY(qnh(:,:)>0._wp)) THEN
            qnh(:,:)  = 0.0_wp
            WRITE (message_text,'(1X,A)') '  qnh > 0, after cloud_twomoment for cloud_type < 2000'
            CALL message(routine,TRIM(message_text))
            CALL finish(TRIM(routine),'Error in two_moment_mcrph')
         END IF
      END IF
      IF (msg_level>dbg_level) CALL message(TRIM(routine), " test for negative values")
      IF (MINVAL(q_cloud) < 0.0) THEN
         WRITE (message_text,'(1X,A)') '  q_cloud < 0, after cloud_twomoment'
         CALL message(routine,TRIM(message_text))
         CALL finish(TRIM(routine),'Error in two_moment_mcrph')
      ENDIF
      IF (MINVAL(q_rain) < 0.0) THEN
         WRITE (message_text,'(1X,A)') '  q_rain < 0, after cloud_twomoment'
         CALL message(routine,TRIM(message_text))
         CALL finish(TRIM(routine),'Error in two_moment_mcrph')
      ENDIF
      IF (MINVAL(q_ice) < 0.0) THEN
         WRITE (message_text,'(1X,A)') '  q_ice < 0, after cloud_twomoment'
         CALL message(routine,TRIM(message_text))
         CALL finish(TRIM(routine),'Error in two_moment_mcrph')
      ENDIF
      IF (MINVAL(q_snow) < 0.0) THEN
         WRITE (message_text,'(1X,A)') '  q_snow < 0, after cloud_twomoment'
         CALL message(routine,TRIM(message_text))
         CALL finish(TRIM(routine),'Error in two_moment_mcrph')
      ENDIF
      IF (MINVAL(q_graupel) < 0.0) THEN
         WRITE (MESSAGE_TEXT,'(1X,A)') '  q_graupel < 0, after cloud_twomoment'
         CALL message(routine,TRIM(message_text))
         CALL finish(TRIM(routine),'Error in two_moment_mcrph')
      ENDIF
      IF (MINVAL(q_hail) < 0.0) THEN
         WRITE (MESSAGE_TEXT,'(1X,A)') '  q_hail < 0, after cloud_twomoment'
         CALL message(routine,TRIM(message_text))
         CALL finish(TRIM(routine),'Error in two_moment_mcrph')
      ENDIF
      IF (MINVAL(n_cloud) < 0.0) THEN
         WRITE (message_text,'(1X,A)') '  n_cloud < 0, after cloud_twomoment'
         CALL message(routine,TRIM(message_text))
         CALL finish(TRIM(routine),'Error in two_moment_mcrph')
      ENDIF
      IF (MINVAL(n_rain) < 0.0) THEN
         WRITE (message_text,'(1X,A)') '  n_rain < 0, after cloud_twomoment'
         CALL message(routine,TRIM(message_text))
         CALL finish(TRIM(routine),'Error in two_moment_mcrph')
      ENDIF
      IF (MINVAL(n_ice) < 0.0) THEN
         WRITE (message_text,'(1X,A)') '  n_ice < 0, after cloud_twomoment'
         CALL message(routine,TRIM(message_text))
         CALL finish(TRIM(routine),'Error in two_moment_mcrph')
      ENDIF
      IF (MINVAL(n_snow) < 0.0) THEN
         WRITE (message_text,'(1X,A)') '  n_snow < 0, after cloud_twomoment'
         CALL message(routine,TRIM(message_text))
         CALL finish(TRIM(routine),'Error in two_moment_mcrph')
      ENDIF
      IF (MINVAL(n_graupel) < 0.0) THEN
         WRITE (message_text,'(1X,A)') '  n_graupel < 0, after cloud_twomoment'
         CALL message(routine,TRIM(message_text))
         CALL finish(TRIM(routine),'Error in two_moment_mcrph')
      ENDIF
      IF (MINVAL(n_hail) < 0.0) THEN
         WRITE (message_text,'(1X,A)') '  n_hail < 0, after cloud_twomoment'
         CALL message(routine,TRIM(message_text))
         CALL finish(TRIM(routine),'Error in two_moment_mcrph')
      ENDIF
    END subroutine check_clouds
    
#endif

  END SUBROUTINE two_moment_mcrph

!===========================================================================================

  SUBROUTINE two_moment_mcrph_init()
     
    INTEGER :: unitnr

#ifndef __SX__

    CALL init_2mom_scheme(cloud_type)

    unitnr = 11
    CALL init_dmin_wetgrowth('dmin_wetgrowth_lookup.dat', unitnr)

    CALL init_dmin_wg_gr_ltab_equi(&
         'dmin_wetgrowth_lookup.dat', &
         unitnr, 61, ltabdminwgg)

    CALL message ("", " finished init_dmin_wetgrowth")

#endif

  END SUBROUTINE two_moment_mcrph_init

!==============================================================================

END MODULE mo_mcrph_sb
