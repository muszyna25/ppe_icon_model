!>
!!=======================================================================
!!
!!  Function for setting up idealized tropical cyclone initial conditions
!!
!!  Given a point specified by: 
!!
!!  	longitude (radians) 
!! 	latitude (radians) 
!! 	height
!!
!!  the functions will return:
!!	u	zonal wind (m s^-1)
!!	v	meridional wind (m s^-1)
!!	w	vertical velocity (m s^-1)
!!	t	temperature (K)
!!	phis	surface geopotential (m^2 s^-2)
!!	ps	surface pressure (Pa)
!!	rho	density (kj m^-3)
!!	q	specific humidity (kg/kg)
!!	qi	tracers (kg/kg)
!!      p       pressure (Pa)  
!!
!!  Initial data are currently identical to:
!!
!!                 Reed, K. A., and C. Jablonowski, 2011: An analytic
!!                 vortex initialization technique for idealized tropical
!!                 cyclone studies in AGCMs. Mon. Wea. Rev., 139, 689-710. 
!!
!!  Author: Kevin A. Reed (University of Michigan, kareed@umich.edu)
!!          version 1
!!          5/21/2012
!!
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Kevin A. Reed, University of Michigan, kareed@umich.edu
!!
!!
!! @par Revision History
!! Adaptation for ICON by Marco Giorgetta (2012-07-02)
!!
!! @par Copyright
!! 2002-2012 by DWD and MPI-M
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
MODULE mo_nh_dcmip_tc

  USE mo_kind,                ONLY: wp

  USE mo_math_constants,      ONLY: deg2rad                   ! convert degree to radian

  USE mo_physical_constants,  ONLY: rd                     ,& ! ideal gas const dry air (J kg^-1 K^1)
    &                               rd_o_cpd               ,& ! rd/cpd
    &                               p0ref                  ,& ! reference pressure for Exner function
    &                               g => grav                 ! gravity (m s^2)
  USE mo_grid_config,         ONLY: grid_rescale_factor    ,& ! factor for small planet scaling (=1/x)
    &                               a => grid_sphere_radius   ! Earth radius/x (m)

  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_intp_data_strc,      ONLY: t_int_state

  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: iqv

  USE mo_intp,                ONLY: cells2edges_scalar


  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: init_nh_dcmip_tc

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !-----------------------------------------------------------------------
  !     Tropical Cyclone Test Case Tuning Parameters
  !-----------------------------------------------------------------------

  REAL(wp), PARAMETER :: rp1        = 282000.0_wp       ,& ! (m)   Radius for calculation of PS
    &                    dp         = 1115.0_wp         ,& ! (Pa)  Delta P for calculation of PS
    &                    zp         = 7000.0_wp         ,& ! (m)   Height for calculation of P
    &                    q0         = 0.0210_wp         ,& !       q at surface from Jordan
    &                    gamma      = 0.0070_wp         ,& ! (K/m) Lapse rate
    &                    ts0        = 302.15_wp         ,& ! (K)   Surface temperature (SST)
    &                    p00        = 101500.0_wp       ,& ! (Pa)  Global mean surface pressure
    &                    p0         = 100000.0_wp       ,& ! (Pa)  p for model level calculation
    &                    cen_lat    =  10.0_wp*deg2rad  ,& ! (rad) Center latitude of initial vortex
    &                    cen_lon    = 180.0_wp*deg2rad  ,& ! (rad) Center longitufe of initial vortex
    &                    zq1        = 3000.0_wp         ,& ! (m)   Height 1 for q calculation
    &                    zq2        = 8000.0_wp         ,& ! (m)   Height 2 for q calculation
    &                    exppr      = 1.5_wp            ,& !       Exponent for r dependence of p
    &                    exppz      = 2.0_wp            ,& !       Exponent for z dependence of p
    &                    ztrop      = 15000.0_wp        ,& ! (m)   Tropopause Height
    &                    qtrop      = 1.e-11_wp         ,& !       Tropopause specific humidity
    &                    consttv    = 0.6080_wp         ,& !       Constant for Virtual Temp Conversion
    &                    epsilon    = 1.e-25_wp            !       Small number to aviod dividing by zero in wind calc
                   
  REAL(wp), PARAMETER :: exponent   = rd*gamma/g                        ,& !       Exponent
    &                    t0         = ts0*(1.0_wp+consttv*q0)           ,& ! (K)   Surface    virtual temp (eq.93)
    &                    ttrop      = t0 - gamma*ztrop                  ,& ! (K)   Tropopause virtual temp (eq.92)
    &                    ptrop      = p00*(ttrop/t0)**(1.0_wp/exponent)    ! (Pa)  Tropopause pressure     (eq.95)

CONTAINS

  !>
  !! <Short description of the subroutine for listings and indices>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !! Adaptation for ICON by Marco Giorgetta (2012-07-02)
  !!
  SUBROUTINE init_nh_dcmip_tc( p_patch, p_nh_prog, p_nh_diag, p_metrics, p_int )

    TYPE(t_patch),         INTENT(inout), TARGET :: p_patch    !< patch on which computation is performed
    TYPE(t_nh_prog),       INTENT(inout)         :: p_nh_prog  !< prognostic state vector
    TYPE(t_nh_diag),       INTENT(inout)         :: p_nh_diag  !< diagnostic state vector
    TYPE(t_nh_metrics),    INTENT(in)            :: p_metrics  !< NH metrics state
    TYPE(t_int_state),     INTENT(in)            :: p_int      !< interpolation state

    ! block dimensions and indices

    INTEGER  :: nblks_c   !< number of blocks of cells
    INTEGER  :: nblks_e   !< number of blocks of edges
    INTEGER  :: jb        !< index of block

    INTEGER  :: npromz_c  !< length of last row in block of cells      (<= nproma)
    INTEGER  :: npromz_e  !< length of last row in block of edges      (<= nproma)
    INTEGER  :: nlen      !< length of loop in row           ( = nproma or npromz)

    INTEGER  :: nlev      !< number of full levels
    INTEGER  :: jk        !< index of level


    ! Arguments of original DCMIP subroutine
    !
    REAL(wp), POINTER     :: lon(:,:)   !< Longitude (radians)
    REAL(wp), POINTER     :: lat(:,:)   !< Latitude (radians)
    REAL(wp), POINTER     :: z(:,:,:)   !< Height (m)
    !
    REAL(wp), POINTER     :: p(:,:,:)   !< Pressure  (Pa)
    !
!!$    REAL(wp), POINTER     :: u(:,:,:)   !< Zonal wind (m s^-1)
!!$    REAL(wp), POINTER     :: v(:,:,:)   !< Meridional wind (m s^-1)
    !
    REAL(wp), POINTER     :: w(:,:,:)   !< Vertical velocity (m s^-1)
    !
    REAL(wp), POINTER     :: t(:,:,:)   !< Temperature (K)
    !
!!$    REAL(wp), POINTER     :: phis(:,:)  !< Surface Geopotential (m^2 s^-2)
!!$    REAL(wp), POINTER     :: ps(:,:)    !< Surface Pressure (Pa)
    !
    REAL(wp), POINTER     :: rho(:,:,:) !< density (kg m^-3)
    REAL(wp), POINTER     :: q(:,:,:)   !< Specific Humidity (kg/kg)


    ! fields to be computed for ICON as a function of the fields above
    !
    REAL(wp), POINTER     :: exner(:,:,:)      !< exner pressure ()
    REAL(wp), POINTER     :: theta_v(:,:,:)    !< virtual potential temperature (K)
    REAL(wp), POINTER     :: rhotheta_v(:,:,:) !< virtual potential temperature (K)
    REAL(wp), POINTER     :: vn(:,:,:)         !< normal wind (m s^-1)


    ! additional fields
    !
    REAL(wp)              :: rp                        !< rp1/x, radius for PS scaled by x
    !
    REAL(wp), POINTER     :: f(:,:)                    !< scaled Coriolis parameter f*x
    REAL(wp), ALLOCATABLE :: gr(:,:)                   !< great dcircle distance to center of vortex for eq.97
    REAL(wp), ALLOCATABLE :: d1(:,:), d2(:,:), d(:,:)  !< factors for computing ufac anv vfac for eqs.109,110,111
    REAL(wp), ALLOCATABLE :: ufac(:,:), vfac(:,:)      !< factors for computing u anv v from tangential wind for eq.112,113
!!$    REAL(wp), ALLOCATABLE :: u(:,:,:)                  !< zonal      wind of cyclone for eq.108
!!$    REAL(wp), ALLOCATABLE :: v(:,:,:)                  !< meridional wind of cyclone for eq.108
    REAL(wp), ALLOCATABLE :: vtc_c(:,:,:)              !< tangential wind of cyclone at cell center for eq.108
    REAL(wp), ALLOCATABLE :: vtc_e(:,:,:)              !< tangential wind of cyclone at edge center for eq.112,113


    ! number of blocks (nblks) and length of loop in last block (npromz)
    !
    ! - for cell variables
    nblks_c  = p_patch%nblks_c
    npromz_c = p_patch%npromz_c
    !
    ! - for edge variables
    nblks_e  = p_patch%nblks_e
    npromz_e = p_patch%npromz_e

    ! number of full levels
    nlev = p_patch%nlev

    ! allocate additional variables
    !
    ! - in cell centers
    ALLOCATE(gr(nproma,nblks_c))
    ALLOCATE(vtc_c(nproma,nlev,nblks_c))
    !
    ! - in edge centers
    ALLOCATE(d1(nproma,nblks_e), d2(nproma,nblks_e), d(nproma,nblks_e))
    ALLOCATE(ufac(nproma,nblks_e), vfac(nproma,nblks_e))
    ALLOCATE(vtc_e(nproma,nlev,nblks_e))

    ! initialize
    !
    gr(:,:)      = 0.0_wp
    !
    d1(:,:)      = 0.0_wp
    d2(:,:)      = 0.0_wp
    d(:,:)       = 0.0_wp
    !
    ufac(:,:)    = 0.0_wp
    vfac(:,:)    = 0.0_wp
    !
    vtc_c(:,:,:) = 0.0_wp
    vtc_e(:,:,:) = 0.0_wp

    ! 1. initialize variables defined at cell centers
    ! ===============================================

    ! connect link single varaibles to components of data structures
    !
    ! horizontal coordinates, heights and geopotential
    lon        => p_patch%cells%center%lon      ! (rad)
    lat        => p_patch%cells%center%lat      ! (rad)
    z          => p_metrics%z_mc                ! (m)
    !
    ! Coriolis parameter
    f          => p_patch%cells%f_c             ! (1/s)
    !
    ! thermodynamics
!!$    ps         => p_nh_diag%pres_sfc            ! (Pa)
    p          => p_nh_diag%pres                ! (Pa)
    exner      => p_nh_prog%exner               ! (),        prognostic
    rho        => p_nh_prog%rho                 ! (kg/m3),   prognostic
    theta_v    => p_nh_prog%theta_v             ! (K),       prognostic, to be computed from p, T, and q
    rhotheta_v => p_nh_prog%rhotheta_v          ! (K*kg/m3), prognostic
    q          => p_nh_prog%tracer(:,:,:,iqv)   ! (kg/kg),   prognostic
    t          => p_nh_diag%temp                ! (K)
    !
    ! vertical wind
    w          => p_nh_prog%w                   ! (m/s),     prognostic


    ! scale radius for pressure field equally as the Earth radius of the grid
    rp = rp1 * grid_rescale_factor

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF

      ! great circle distance (gr) to vortex center                    (eq.97)
      !-----------------------------------------------------------------------

      gr(1:nlen,jb) = a*ACOS( SIN(cen_lat)*SIN(lat(1:nlen,jb))                               &
        &                    +COS(cen_lat)*COS(lat(1:nlen,jb))*COS(lon(1:nlen,jb)-cen_lon) )


      ! initialize atmosphere ...
      !-----------------------------------------------------------------------
      DO jk = 1, nlev

        WHERE (z(1:nlen,jk,jb) > ztrop) ! ... above tropopause ...

          ! pressure                                                   (eqs.94+98)
          !-----------------------------------------------------------------------
          p(1:nlen,jk,jb) = ptrop*EXP(-(g*(z(1:nlen,jk,jb)-ztrop))/(rd*ttrop))

          ! virtual potential temperature                         (eqs.93,104,105)
          !-----------------------------------------------------------------------
          theta_v(1:nlen,jk,jb) = ttrop

          ! specific humidity                                              (eq.91)
          !-----------------------------------------------------------------------
          q(1:nlen,jk,jb) = qtrop

          ! tangential wind of cyclone                                    (eq.108)
          !-----------------------------------------------------------------------
          vtc_c(1:nlen,jk,jb) = 0.0_wp

        ELSEWHERE ! ... below tropopause ...

          ! pressure                                                   (eqs.94+98)
          !-----------------------------------------------------------------------
          p(1:nlen,jk,jb) = (p00-dp*EXP(-(gr(1:nlen,jb)/rp)**exppr-(z(1:nlen,jk,jb)/zp)**exppz)) &
            &              *((t0-gamma*z(1:nlen,jk,jb))/t0)**(1/exponent)

          ! virtual potential temperature                         (eqs.93,104,105)
          !-----------------------------------------------------------------------
          theta_v(1:nlen,jk,jb) =  (t0-gamma*z(1:nlen,jk,jb))                                                 &
            &                     /(1.0_wp + (exppz*rd*(t0-gamma*z(1:nlen,jk,jb))*z(1:nlen,jk,jb))            &
            &                               /(g*zp**exppz*(1.0_wp-p00/dp*EXP( (gr(1:nlen,jb)/rp)**exppr       &
            &                                                                +(z(1:nlen,jk,jb)/zp)**exppz))))

          ! specific humidity                                              (eq.91)
          !-----------------------------------------------------------------------
          q(1:nlen,jk,jb) = q0*EXP(-(z(1:nlen,jk,jb)/zq1)-(z(1:nlen,jk,jb)/zq2)**exppz)

          ! tangential wind of cyclone                                    (eq.108)
          !-----------------------------------------------------------------------
          vtc_c(1:nlen,jk,jb) = -f(1:nlen,jb)*gr(1:nlen,jb)/2.0_wp                                                    &
            &                   +SQRT( (f(1:nlen,jb)*gr(1:nlen,jb)/2.0_wp)**(2.0_wp)                                  &
            &                         - ( exppr*(gr(1:nlen,jb)/rp)**exppr*rd*(t0-gamma*z(1:nlen,jk,jb)))              &
            &                          /( 1.0_wp+exppz*rd*(t0-gamma*z(1:nlen,jk,jb))*z(1:nlen,jk,jb)/(g*zp**exppz)    &
            &                            -p00/dp*EXP((gr(1:nlen,jb)/rp)**exppr+(z(1:nlen,jk,jb)/zp)**exppz)))

        END WHERE

        ! ... and at all levels

        ! Exner pressure
        !-----------------------------------------------------------------------
        exner(1:nlen,jk,jb) = (p(1:nlen,jk,jb)/p0ref)**rd_o_cpd

        ! density of moist air                                          (eq.106)
        !-----------------------------------------------------------------------
        rho(1:nlen,jk,jb) = p(1:nlen,jk,jb)/(rd*theta_v(1:nlen,jk,jb) )

        ! density*virtual potential temparature
        !-----------------------------------------------------------------------
        rhotheta_v(1:nlen,jk,jb) = rho(1:nlen,jk,jb)*theta_v(1:nlen,jk,jb)

        ! temperature                                           (eqs.93,104,105)
        !-----------------------------------------------------------------------
        t(1:nlen,jk,jb) = theta_v(1:nlen,jk,jb)/(1.0_wp+consttv*q(1:nlen,jk,jb))

        ! vertical velocity
        !-----------------------------------------------------------------------
        w(1:nlen,jk,jb) = 0.0_wp

!!$        ! surface geopotential
!!$        !-----------------------------------------------------------------------
!!$        phis(1:nlen,jb) = 0.0_wp

      END DO !jk
    END DO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! 2. interpolate from cell centers to edge centers
    ! ================================================

    CALL cells2edges_scalar(vtc_c,                 &
      &                     p_patch,p_int%c_lin_e, &
      &                     vtc_e                  )


    ! 3. initialize variables defined at edge centers
    ! ===============================================

    ! connect link single varaibles to components of data structures
    !
    ! horizontal coordinates, heights and geopotential
    lon        => p_patch%edges%center%lon      ! (rad)
    lat        => p_patch%edges%center%lat      ! (rad)
    !
    ! horizontal wind normal to edges
    vn         => p_nh_prog%vn                  ! (m/s),     prognostic

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF

      ! factors for u and v components of the tangential wind    (eqs.109-111)
      !-----------------------------------------------------------------------

      d1(1:nlen,jb) =  SIN(cen_lat)*COS(lat(1:nlen,jb))                             &
        &             -COS(cen_lat)*SIN(lat(1:nlen,jb))*COS(lon(1:nlen,jb)-cen_lon)
      d2(1:nlen,jb) =  COS(cen_lat)*SIN(lon(1:nlen,jb)-cen_lon)
      d (1:nlen,jb) =  MAX(epsilon, SQRT(d1(1:nlen,jb)**2.0_wp + d2(1:nlen,jb)**2.0_wp))
      
      ufac(1:nlen,jb) = d1(1:nlen,jb)/d(1:nlen,jb)
      vfac(1:nlen,jb) = d2(1:nlen,jb)/d(1:nlen,jb)


      DO jk = 1, nlev

        ! compute wind normal to edges from
        ! tangential wind of cyclone at edges                  (eqs.108,112,113)
        !-----------------------------------------------------------------------

        vn(1:nlen,jk,jb) = vtc_e(1:nlen,jk,jb)*ufac(1:nlen,jb) * p_patch%edges%primal_normal(1:nlen,jb)%v1   &
          &               +vtc_e(1:nlen,jk,jb)*vfac(1:nlen,jb) * p_patch%edges%primal_normal(1:nlen,jb)%v2

      END DO !jk
    END DO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DEALLOCATE(gr)
    DEALLOCATE(d1, d2, d)
    DEALLOCATE(ufac, vfac)
    DEALLOCATE(vtc_c, vtc_e)

  END SUBROUTINE init_nh_dcmip_tc

END MODULE mo_nh_dcmip_tc

