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
!! Description:
!!
!! The subroutine is the interface to the original KAMM2 modules,
!! which are provided by src_seifert.f90. A major difference between UCLA-LES and
!! KAMM2 is that KAMM2 uses mass densities instead of mixing ratios, thus
!! q_cloud=rho*qc.
!! Temporary allocation of memory to the KAMM2 variables is done by the
!! subroutines ALLOC_DRIVER and ALLOC_WOLKEN.
!! All microphysical source terms e.g. nucleation, condensation, coagulation,
!! freezing and melting are then calculated and time integrated within the
!! subroutine CLOUDS.
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
    alv    , & !! latent heat of vaporization
    als    , & !! latent heat of sublimation
    cpdr  => rcpd   , & !! (spec. heat of dry air at constant press)^-1
    cvdr  => rcvd       !! (spec. heat of dry air at const vol)^-1

USE mo_satad          , ONLY: satad_v_3d     , & !! new saturation adjustment
                              sat_pres_water , & !! saturation vapor pressure w.r.t. water
                              sat_pres_ice   , & !! saturation vapor pressure w.r.t. ice
                              spec_humi          !! Specific humidity 

USE mo_exception      , ONLY: message, message_text

USE wolken_driver,     ONLY: loc_ix, loc_iy, loc_iz,           &
     &                       dt_kamm2 => dt,                   &
     &                       dz_kamm2 => dz3d,                 &
     &                       w_kamm2 => w,                     &
     &                       T_kamm2 => T,                     &
     &                       p_0, T_0, rho_k=>rho_0, q,        &
     &                       q_cloud, n_cloud,                 &
     &                       q_ice, q_rain, q_snow, q_graupel, &
     &                       n_ice, n_rain, n_snow, n_graupel, &
     &                       prec_ice, prec_rain,              &
     &                       prec_snow, prec_graupel, prec,    &
     &                       n_hail, q_hail, prec_hail,        &
     &                       alloc_driver, dealloc_driver, cp, &
     &                       S_w, S_i,dSwdz, dSidz, dT0dz,     &
     &                       zml_k,                            &
     &                       w_cb,                             &
     &                       speichere_umwandlungsraten,       &
     &                       speichere_precipstat,             &
     &                       dmin_wg_g, pvec_wg_g, Tvec_wg_g,  &
     &                       qwvec_wg_g, qivec_wg_g, &
     &                       anzp_wg, anzT_wg, anzi_wg, anzw_wg, &
     &                       ltabdminwgg, cloud_type => wolke_typ


USE wolken,            ONLY: alloc_wolken, dealloc_wolken,     &
     &                       clouds,                           &
     &                       rrho_04, rrho_c

USE wolken_eis, ONLY : init_dmin_wetgrowth, &
                       init_dmin_wg_gr_ltab_equi

USE wolken_konstanten, ONLY: init_seifert,                     &
                             L_ew,L_ed,R_d,R_l,                &
                             myparticle=>particle,             &
                             mycloud=>cloud,myrain=>rain,myice=>ice,&
                             mysnow=>snow,mygraupel=>graupel,myhail=>hail, &
                             e_ws,e_es,                        &
                             e_ws_vec,e_es_vec,                &
                             T_nuc, T_f, satad_nach_mikrophysik, & 
                             graupel_shedding, hail_shedding,    &
                             qnc_const

 USE mo_sync,           ONLY: global_max
 USE wolken_sedi


!==============================================================================

IMPLICIT NONE
PUBLIC 

CHARACTER(len=*), PARAMETER, PRIVATE :: &
  &  version = '$Id$'

CONTAINS

  !==============================================================================
  !
  ! Two-moment mixed-phase bulk microphysics
  !
  ! original version by Axel Seifert, May 2003
  ! with modifications by Ulrich Blahak, August 2007

  !==============================================================================
  SUBROUTINE two_moment_mcrph(            &
                       isize,             &!in:array size
                       ke,                &!in: end level/array size
                       is,                &!in: start index, optional
                       ie,                &!in: end index, optional
                       ks,                &!in: start index vertical , optional 
                       dt,                &!in: time step
                       dz,                &!in: vertical layer thickness
                       rho,               &!in: density
                       pres,              &!in: pressure
                       qv,                &!inout: sp humidity
                       qc,                &!inout: cloud water
                       qr, qnr,           &!inout: rain
                       qi, qni,           &!inout: ice
                       qs, qns,           &!inout: snow
                       qg, qng,           &!inout: graupel
                       qh, qnh,           &!inout: hail
                       tk,                &!inout: temp 
                       w,                 &!inout: w
                       prec_r,            &!inout precp rate rain
                       prec_i,            &!inout precp rate ice
                       prec_s,            &!inout precp rate snow
                       prec_g,            &!inout precp rate graupel
                       prec_h,            &!inout precp rate hail
                       msg_level,         &!in msg_level
                       l_cv          )    
                

    ! Declare variables in argument list

    INTEGER,            INTENT (IN)  :: isize, ke    !grid sizes
    INTEGER,  OPTIONAL, INTENT (IN)  :: is, ie, ks   !start/end indices

    REAL(wp), INTENT (IN)  :: dt        !time step
    REAL(wp), DIMENSION(:,:), INTENT(IN) :: dz, rho, pres
     
    !Tracers
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: qv, qc, qr, qnr, qi,   &
                                               qni, qs, qns, qg, qng, &
                                               qh, qnh

    !Dynamic variables
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: tk, w 


    !Precp rates
    REAL(wp), DIMENSION(:), INTENT (INOUT)  :: prec_r, prec_i, &
                                               prec_s, prec_g, prec_h

    INTEGER,  INTENT (IN)             :: msg_level 
    LOGICAL,  OPTIONAL,  INTENT (IN)  :: l_cv

    ! ... Local Variables

    INTEGER    :: its,ite,kts,kte
    INTEGER    :: i,j,k,ii,kk,n,nt_kamm2,ntsedi,igridpoints,izstat,kkk
    
    REAL(wp) :: q_vap_new,q_vap_old
    REAL(wp) :: q_liq_new,q_liq_old
    REAL(wp) :: q_ice_new,q_ice_old
    REAL(wp) :: q_v, hlp, hlp2, hlp3, dt_sedi,e_v,T_a,qerr, zxi, zxc
    REAL(wp) :: convliq,convice
    REAL(wp) :: z_heat_cap_r !! reciprocal of cpdr or cvdr (depending on l_cv)
    REAL(wp) :: rdz(isize,ke), rhocorr(isize,ke)

    INTEGER, DIMENSION(:), AllOCATABLE    :: ilm,klm
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ssi 

    ! using eps=0 increases runtime, but may remove artifacts
    REAL,    PARAMETER :: eps = 0._wp
    INTEGER, PARAMETER :: dbg_level = 25 !level for debug prints
    LOGICAL, PARAMETER :: CGP_SEARCH = .true.
    LOGICAL, PARAMETER :: debug      = .false.    
    CHARACTER(len=*), PARAMETER :: routine = 'mo_mcrph_sb:two_moment_mcrph'

    !inverse of vertical layer thickness
    rdz = 1._wp / dz

    !start/end indices
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

    !..constant cloud droplet number
    qnc_const = 200.0e6_wp

    loc_ix = (kte-kts+1)*(ite-its+1)

    IF (msg_level>dbg_level)THEN
       WRITE (message_text,'(1X,A,I4,A,E12.4)') "mcrph_sb: cloud_type = ",cloud_type,", CCN = ",qnc_const
       CALL message(TRIM(routine),TRIM(message_text))
    END IF


    ALLOCATE( ilm(0:loc_ix), klm(0:loc_ix))
    ALLOCATE( ssi(its:ite,kts:kte))

    !..calculate supersaturations after dynamics/advection
    DO ii = its, ite
       DO kk = kts, kte
           ! supersaturation w.r.t. ice
           ssi(ii,kk) = R_d  * rho(ii,kk) * qv(ii,kk) * &
                        tk(ii,kk) / sat_pres_ice(tk(ii,kk)) - 1.0
       ENDDO
    ENDDO

    !..search for cloudy grid points and store locations
    i = -1
    IF (msg_level>dbg_level) CALL message(TRIM(routine), "mcrph_sb: search for cloudy grid points")

    IF (CGP_SEARCH) THEN
          DO ii = its, ite
             DO kk = kts, kte
                IF (  ssi(ii,kk)  > eps .or. &
                     & qc(ii,kk)  > eps .or. &
                     & qr(ii,kk)  > eps .or. &
                     & qi(ii,kk)  > eps .or. &
                     & qs(ii,kk)  > eps .or. &
                     & qg(ii,kk)  > eps .or. &
                     & qh(ii,kk)  > eps ) THEN
                   i = i+1
                   ilm(i) = ii     ! they start with i=0
                   klm(i) = kk
                ENDIF
             ENDDO
          ENDDO
    ELSE
          DO ii = its, ite
             DO kk = kts, kte
                i = i+1
                ilm(i) = ii
                klm(i) = kk
             ENDDO
          ENDDO
    ENDIF
    igridpoints = i

    !..return now, if no clouds are found
    IF (igridpoints == -1) THEN

       IF(msg_level>dbg_level)CALL message(TRIM(routine), "mcrph_sb: no clouds found!")

       DEALLOCATE(ssi)
       DEALLOCATE(ilm,klm)

    ELSE ! cloudy points have been found

       !IF (debug.and.isIO()) WRITE (0,'(a,i5,a,i10,a,f6.2)') &
       !     & " mcrph_sb: myid = ",myid,", number of cloudy grid points = ",i+1, &
       !     & ", percent of domain = ", 100.0*(i+1.0)/loc_ix

       loc_ix = igridpoints
       loc_iy = 1
       loc_iz = 1

       !IF (debug.and.isIO()) THEN
       !   WRITE (*,*) "mcrph_sb: grid"
       !   WRITE (*,*) "      its = ",its
       !   WRITE (*,*) "      ite = ",ite
       !   WRITE (*,*) "      jts = ",jts
       !   WRITE (*,*) "      jte = ",jte
       !   WRITE (*,*) "      kts = ",kts
       !   WRITE (*,*) "      kte = ",kte
       !   WRITE (*,*) "      loc_ix = ",loc_ix
       !   WRITE (*,*) "      loc_iy = ",loc_iy
       !   WRITE (*,*) "      loc_iz = ",loc_iz
       !ENDIF

       IF (msg_level>dbg_level) CALL message(TRIM(routine), "mcrph_sb: calling allocation")

       ! ... Allocate memory to temporary KAMM2 variables
       CALL alloc_driver()
       CALL alloc_wolken()

       IF (msg_level>dbg_level) CALL message(TRIM(routine), "start copy to one-dimensional array")

       ! ... transpose to one-dimensional array and variables used in cloud module
       j = 1
       k = 1
       do i=0,loc_ix
          ! ... grid points
          ii = ilm(i)
          kk = klm(i)
           
          ! ... dynamics
          T_0(i,j,k)      = tk(ii,kk)
          p_0(i,j,k)      = pres(ii,kk)
          rho_k(i,j,k)    = rho(ii,kk)

          ! .. the ice supersaturation
          S_i(i,j,k)       = ssi(ii,kk)

          ! ... concentrations --> number densities
          n_rain(i,j,k)    = rho_k(i,j,k) * qnr(ii,kk) 
          n_ice(i,j,k)     = rho_k(i,j,k) * qni(ii,kk) 
          n_snow(i,j,k)    = rho_k(i,j,k) * qns(ii,kk)
          n_graupel(i,j,k) = rho_k(i,j,k) * qng(ii,kk)
          n_hail(i,j,k)    = rho_k(i,j,k) * qnh(ii,kk)

          ! ... mixing ratios -> mass densities
          q(i,j,k)         = rho_k(i,j,k) * qv(ii,kk) 
          q_cloud(i,j,k)   = rho_k(i,j,k) * qc(ii,kk) 
          q_rain(i,j,k)    = rho_k(i,j,k) * qr(ii,kk) 
          q_ice(i,j,k)     = rho_k(i,j,k) * qi(ii,kk) 
          q_snow(i,j,k)    = rho_k(i,j,k) * qs(ii,kk) 
          q_graupel(i,j,k) = rho_k(i,j,k) * qg(ii,kk) 
          q_hail(i,j,k)    = rho_k(i,j,k) * qh(ii,kk) 

          w_kamm2(i,j,k)   = w(ii,kk) !+ 0.5 * tke(ii,kk)

       enddo

       IF (msg_level>dbg_level) CALL message(TRIM(routine), "finished copy to one-dimensional array")

       ! ... timestep
       dt_kamm2 = dt

       IF (msg_level>dbg_level) CALL message(TRIM(routine),"mcrph_sb: calling clouds")

       !======================================================================
       ! .. this subroutine calculates all the microphysical sources and sinks

       CALL clouds ()

       !======================================================================

       IF (debug) THEN
          IF (MINVAL(q) < 0.0) THEN
             WRITE (*,*) ' mcrph_sb: q < 0 ENCOUNTERED, FILLED WITH 0.0'
             WHERE (q < 0.0) q = 0.0
          ENDIF
          IF (MINVAL(q_cloud) < 0.0) THEN
             write (*,*) ' mcrph_sb: q_cloud < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(q_rain) < 0.0) THEN
             write (*,*) ' mcrph_sb: q_rain < 0, STOPPED AFTER CLOUDS 1', MINVAL(q_rain)
             stop
          ENDIF
          IF (MINVAL(q_ice) < 0.0) THEN
             write (*,*) ' mcrph_sb: q_ice < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(q_snow) < 0.0) THEN
             write (*,*) ' mcrph_sb: q_snow < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(q_graupel) < 0.0) THEN
             WRITE (*,*) ' mcrph_sb: q_graupel < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(q_hail) < 0.0) THEN
             WRITE (*,*) ' mcrph_sb: q_hail < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(n_cloud) < 0.0) THEN
             write (*,*) ' mcrph_sb: n_cloud < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(n_rain) < 0.0) THEN
             write (*,*) ' mcrph_sb: n_rain < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(n_ice) < 0.0) THEN
             write (*,*) ' mcrph_sb: n_ice < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(n_snow) < 0.0) THEN
             write (*,*) ' mcrph_sb: n_snow < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(n_graupel) < 0.0) THEN
             write (*,*) ' mcrph_sb: n_graupel < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(n_hail) < 0.0) THEN
             write (*,*) ' mcrph_sb: n_hail < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
       ENDIF

       ! ... Transformation of variables back to driving model and latent heat equation
       IF (msg_level>dbg_level) CALL message(TRIM(routine), "mcrph_sb: trafo back")

       j = 1
       k = 1
       DO i=0,loc_ix
          ii = ilm(i)
          kk = klm(i)

          hlp = 1.0 / rho_k(i,j,k)

          ! ... latent heat for temperature equation
          !     in UCLA-LES the qc is part of theta_l, but qr is not

          ! ... ... new variables
          q_vap_new = hlp * q(i,j,k)
          q_liq_new = hlp *  ( q_rain(i,j,k) + q_cloud(i,j,k) )
          q_ice_new = hlp * ( q_ice(i,j,k)+q_snow(i,j,k)+q_graupel(i,j,k)+q_hail(i,j,k) )

          ! ... ... old variables
          q_vap_old = qv(ii,kk)
          q_liq_old = qc(ii,kk) + qr(ii,kk)
          q_ice_old = qi(ii,kk) + qs(ii,kk) + qg(ii,kk) + qh(ii,kk)

          !...update temperature here
          convice = z_heat_cap_r * als
          convliq = z_heat_cap_r * (alv-als)
          tk(ii,kk) = tk(ii,kk) &
               &      - convice * (q_vap_new - q_vap_old)  &
               &      + convliq * (q_liq_new - q_liq_old)

          ! ... mass densities to mixing ratios with actual density:
          qv(ii,kk) = hlp * q(i,j,k)
          qc(ii,kk) = hlp * q_cloud(i,j,k)
          qr(ii,kk) = hlp * q_rain(i,j,k)
          qi(ii,kk) = hlp * q_ice(i,j,k)
          qs(ii,kk) = hlp * q_snow(i,j,k)
          qg(ii,kk) = hlp * q_graupel(i,j,k)
          qh(ii,kk) = hlp * q_hail(i,j,k)

          ! ... number concentrations
          qnr(ii,kk) = hlp * n_rain(i,j,k)
          qni(ii,kk) = hlp * n_ice(i,j,k)
          qns(ii,kk) = hlp * n_snow(i,j,k)
          qng(ii,kk) = hlp * n_graupel(i,j,k)
          qnh(ii,kk) = hlp * n_hail(i,j,k)
       ENDDO

       DEALLOCATE(ilm,klm)

       ! ... deallocation
       IF (msg_level>dbg_level) CALL message(TRIM(routine),"calling deallocation")
       CALL dealloc_driver()
       CALL dealloc_wolken()

       DEALLOCATE(ssi)

       IF (msg_level>dbg_level) CALL message(TRIM(routine),"cloud microphysics end")

    END IF ! ... This ends the loooong if-block 'cloudy points present'


!===========================================================================================
    IF (msg_level>dbg_level) CALL message(TRIM(routine),"sedimentation")
!===========================================================================================

    prec_r  = 0._wp
    rhocorr = SQRT(1.2_wp/rho) 
    CALL rain_sedi_icon (qr, qnr, qc, prec_r, rho, rhocorr, rdz, dt, its, ite, kts, kte)


    IF (cloud_type.ge.1000) THEN

      prec_i = 0._wp
      prec_s = 0._wp
      prec_g = 0._wp

      IF (ANY(qi(:,:)>0._wp)) &
        call ice_sedi_icon (qi, qni, prec_i, rho, rhocorr, rdz, dt, its, ite, kts, kte)

      IF (ANY(qs(:,:)>0._wp)) &
        call snow_sedi_icon (qs, qns, prec_s, rho, rhocorr, rdz, dt, its, ite, kts, kte)

      IF (ANY(qg(:,:)>0._wp)) THEN
        ntsedi = 3
        DO ii=1,ntsedi
          call graupel_sedi_icon (qg, qng, prec_g, rho, rhocorr, rdz, dt/ntsedi, its, ite, kts, kte)
        END DO
      END IF

    END IF

    IF (cloud_type.ge.2000) THEN

      prec_h = 0.0
      IF (ANY(qh(:,:)>0._wp)) THEN
        ntsedi = 3
        DO ii=1,ntsedi
          call hail_sedi_icon (qh, qnh, prec_h, rho, rhocorr, rdz, dt/ntsedi, its, ite, kts, kte)
        END DO
      ENDIF

    END IF

    IF (msg_level>dbg_level) CALL message(TRIM(routine), "two moment mcrph ends!")

    RETURN

  END SUBROUTINE two_moment_mcrph

!===========================================================================================

  SUBROUTINE two_moment_mcrph_init( )
     
    INTEGER :: unitnr

    LOGICAL, PARAMETER :: lgshed_ascm = .false.
    LOGICAL, PARAMETER :: lhshed_ascm = .false.

     
    CALL init_seifert( cloud_type )

    unitnr = 11
    CALL init_dmin_wetgrowth('dmin_wetgrowth_lookup.dat', unitnr)

    CALL init_dmin_wg_gr_ltab_equi(&
         'dmin_wetgrowth_lookup.dat', &
         unitnr, 61, ltabdminwgg)

    CALL message ("", "mcrph_sb: finished init_dmin_wetgrowth")

    graupel_shedding = lgshed_ascm
    hail_shedding    = lhshed_ascm
    speichere_precipstat = .FALSE.
    speichere_umwandlungsraten = .FALSE.


  END SUBROUTINE two_moment_mcrph_init

!==============================================================================

END MODULE mo_mcrph_sb
