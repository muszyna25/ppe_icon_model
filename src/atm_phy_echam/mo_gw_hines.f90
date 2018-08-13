#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
#endif
#include "fsel.inc"
!>
!! Hines parameterization for the vertical transport and dissipation
!! of unresolved gravity wave spectra originating from the troposphere.
!!
!!
!! Authors:
!!
!! @author   c. mclandress   ists   august 1995
!! @author   n. mcfarlane    cccma  may 1995
!! @author   m. charron      mpi-m  2000-2001
!! @author   e. manzini      mpi-m  february 2002 (re-write, based on cccgwd)
!! @author   h. schmidt      mpi-m  march 2003
!! @author   h. schmidt      mpi-m  april 2010 (add latitude dependent source function)
!!
!!
!! @par Revision History
!! - implementation in ICON based on echam-6.0.10, r2480
!!   by Marco Giorgetta, MPI-M (2011-07-28)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_gw_hines

  USE mo_exception,            ONLY: message_text, message, finish

  USE mo_kind,                 ONLY: wp
  USE mo_physical_constants,   ONLY: grav, rd, cpd !!$, re, rhoh2o

  USE mo_echam_gwd_config,     ONLY: echam_gwd_config

  USE mo_math_constants,       ONLY: cos45, one_third
  USE mo_fast_math_lib,        ONLY: vec_cbrt ! cube root

!!$  USE mo_geoloc,               ONLY: ilat
!!$  USE mo_vertical_coord_table, ONLY: vct_a, vct_b
!!$  USE mo_gaussgrid,            ONLY: gl_twomu, gl_sqcst
!!$  USE mo_memory_g1a,           ONLY: dalpslm1, dalpsmm1, vom1, dm1, alpsm1, tm1
!!$  USE mo_memory_g2a,           ONLY: dtlm1, dtmm1, um1, vm1, dudlm1, dvdlm1

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: gw_hines

  !----------------------------------
  ! Internal switches and constants !
  ! ---------------------------------

  INTEGER  :: naz        = 8

  REAL(wp) :: slope      = 1.0_wp

  REAL(wp) :: f1         = 1.5_wp
  REAL(wp) :: f2         = 0.3_wp
  REAL(wp) :: f3         = 1.0_wp
  REAL(wp) :: f5         = 1.0_wp
  REAL(wp) :: f6         = 0.5_wp

  INTEGER  :: icutoff    = 0
  REAL(wp) :: alt_cutoff = 105.e3_wp

  REAL(wp) :: smco       = 2.0_wp      !  (test value: smco = 1.0)
  INTEGER  :: nsmax      = 5           !  (test value: nsmax = 2)

CONTAINS

  SUBROUTINE gw_hines ( jg         ,&! in,  grid level/domain index
    &                   nbdim      ,&! in,  dimension of block of cells/columns
    &                   jcs        ,&! in,  start index of loops over cells/columns
    &                   jce        ,&! in,  end   index ...
    &                   nc         ,&! in,  number of cells/columns in loop (jce-jcs+1)
    &                   nlev       ,&! in,  number of levels
    &                   paphm1     ,&! in,  p at half levels       [Pa]
    &                   papm1      ,&! in,  p at full levels       [Pa]
    &                   pzh        ,&! in,  half level height asl. [m]
    &                   prho       ,&! in,  full level density     [kg/m3]
    &                   pmair      ,&! in,  air mass in layer      [kg/m2]
    &                   ptm1       ,&! in,  T                      [K]
    &                   pum1       ,&! in,  u                      [m/s]
    &                   pvm1       ,&! in,  v                      [m/s]
!!$    &                   lat_deg    ,&! in,  latitude               [degN]
!!$    &                   paprflux   ,&! in, precipitation flux at surface
    &                   dissip_gwd ,&! out,     Q|Hines            [W/kg]
    &                   tend_u_gwd ,&! out, du/dt|Hines            [m/s2]
    &                   tend_v_gwd ) ! out, dv/dt|Hines            [m/s2]


    !
    ! Description:
    !
    !   Hines parameterization from ccc/mam (Hines, 1997a,b):
    !   physical tendencies of the prognostic variables u,v
    !   due to vertical transports by a broad band spectrum
    !   of gravity waves.
    !
    !   Note that diffusion coefficient and  heating rate
    !   only calculated if lheatcal = .TRUE.
    !
    !        *gw_hines* is called from *physc*.
    !

    IMPLICIT NONE

    ! scalar argument with intent(IN)
    INTEGER  ,INTENT(in)  :: jg
    INTEGER  ,INTENT(in)  :: nbdim
    INTEGER  ,INTENT(in)  :: jcs
    INTEGER  ,INTENT(in)  :: jce
    INTEGER  ,INTENT(in)  :: nc
    INTEGER  ,INTENT(in)  :: nlev

    !  Array arguments with intent(IN):
    ! Input 1D
!!$    REAL(wp) ,INTENT(in)  :: paprflux(nbdim)         ! precipitation flux
!!$    REAL(wp) ,INTENT(in)  :: lat_deg(nbdim)          ! latitude in deg N
    ! Input 2D
    REAL(wp) ,INTENT(in)  :: paphm1(nbdim,nlev+1)    ! half level pressure (t-dt)
    REAL(wp) ,INTENT(in)  :: papm1(nbdim,nlev)       ! full level pressure (t-dt)
    REAL(wp) ,INTENT(in)  :: pzh(nbdim,nlev+1)       ! half level height asl. (m)
    REAL(wp) ,INTENT(in)  :: prho(nbdim,nlev)        ! full level density (kg/m3)
    REAL(wp) ,INTENT(in)  :: pmair(nbdim,nlev)       ! full level air mass (kg/m2)
    REAL(wp) ,INTENT(in)  :: ptm1(nbdim,nlev)        ! temperature (t-dt)
    REAL(wp) ,INTENT(in)  :: pum1(nbdim,nlev)        ! zonal wind (t-dt)
    REAL(wp) ,INTENT(in)  :: pvm1(nbdim,nlev)        ! meridional wind (t-dt)

    !  Array arguments with intent(OUT):
    ! - input/output 2d
    REAL(wp) ,INTENT(out) :: dissip_gwd(nbdim,nlev)  ! gw energy dissipation
    REAL(wp) ,INTENT(out) :: tend_u_gwd(nbdim,nlev)  ! tendency of zonal wind
    REAL(wp) ,INTENT(out) :: tend_v_gwd(nbdim,nlev)  ! tendency of meridional wind

    !  Local arrays for ccc/mam hines gwd scheme:

    ! Important local parameter (passed to all subroutines):
    INTEGER, PARAMETER :: nazmth = 8     ! max azimuth array dimension size

    REAL(wp) :: pressg(nc)               ! Surface pressure (pascal)
!!$    REAL(wp) :: zpr(nc)                  ! precipitation (check dims: echam5 change)

    ! * Vertical positioning arrays and work arrays:
    REAL(wp) :: sgj(nc,nlev)
    REAL(wp) :: shj(nc,nlev)
    REAL(wp) :: shxkj(nc,nlev)
    REAL(wp) :: dttdsf(nc)

    REAL(wp) :: diffco(nc,nlev)          ! diffusion coefficient (m^2/s)

    REAL(wp) :: flux_u(nc,nlev)          ! zonal momentum flux (pascals)
    REAL(wp) :: flux_v(nc,nlev)          ! meridional momentum flux (pascals)

    REAL(wp) :: uhs(nc,nlev)             ! zonal wind (m/s), input for hines param
    REAL(wp) :: vhs(nc,nlev)             ! merid wind (m/s), input for hines param
    REAL(wp) :: bvfreq(nc,nlev)          ! background brunt vassala frequency (rad/s)
    REAL(wp) :: density(nc,nlev)         ! background density (kg/m^3)
    REAL(wp) :: mair(nc,nlev)            ! background airmass (kg/m^2)
    REAL(wp) :: visc_mol(nc,nlev)        ! molecular viscosity (m^2/s)
    REAL(wp) :: alt(nc,nlev)             ! background altitude above ground of lower half level of a layer (m)

    REAL(wp) :: rmswind(nc)              ! rms gravity wave  wind, lowest level (m/s)
    REAL(wp) :: anis(nc,nazmth)          ! anisotropy factor (sum over azimuths = 1)
    REAL(wp) :: k_alpha(nc,nazmth)       ! horizontal wavenumber of each azimuth (1/m)
    LOGICAL  :: lorms(nc)                ! .true. for rmswind /=0 at launching level

    REAL(wp) :: m_alpha(nc,nlev,nazmth)  ! cutoff vertical wavenumber (1/m)
    REAL(wp) :: mmin_alpha(nc,nazmth)    ! minumum value of m_alpha
    REAL(wp) :: sigma_t(nc,nlev)         ! total rms gw wind (m/s)

    ! gw variances from orographic sources (for coupling to a orogwd)
    REAL(wp) :: sigsqmcw(nc,nlev,nazmth)
    REAL(wp) :: sigmatm(nc,nlev)

    !
    ! Local scalars:
    INTEGER  :: jk, jl
    INTEGER  :: levbot     ! gravity wave spectrum lowest level
    REAL(wp) :: rgocp, ratio, pressg_inv
!!$    REAL(wp) :: zpcons

    CHARACTER(len=*), PARAMETER :: routine = 'mo_gw_hines:gw_hines'

  ! Shortcuts to components of echam_gwd_config
  !
  INTEGER , POINTER :: emiss_lev     !< number of levels above the ground at which gw are emitted
  REAL(wp), POINTER :: rmscon        !< [m/s] root mean square gravity wave wind at emission level
  REAL(wp), POINTER :: kstar         !< [1/m] typical gravity wave horizontal wavenumber
  !
!!$  LOGICAL , POINTER :: lfront        !< true: compute gw sources emerging from fronts and background
!!$                                     !< (Charron and Manzini, 2002)
!!$                                     !< for which gravity waves are emitted from fronts
!!$  !
!!$  LOGICAL , POINTER :: lozpr         !< true: for background enhancement associated with precipitation
!!$                                     !< (Manzini et al., 1997)
!!$  REAL(wp), POINTER :: pcrit         !< [mm/d] critical precipitation value, above which
!!$                                     !< gravity wave rms wind enhancement is applied
!!$  REAL(wp), POINTER :: pcons         !< [] adimensional factor for background enhancement
!!$                                     !< associated with precipitation
!!$  !
!!$  LOGICAL , POINTER :: lrmscon_lat   !< true:  use latitude dependent rmscon
!!$                                     !< - |latitude| >= lat_rmscon:
!!$                                     !<      use rmscon
!!$                                     !< - |latitude| <= lat_rmscon_eq:
!!$                                     !<      use rmscon_eq
!!$                                     !< - lat_rmscon_eq < |latitude| < lat_rmscon:
!!$                                     !<      use linear interpolation between rmscon_eq and rmscon
!!$                                     !< false: use rmscon for all latitudes
!!$                                     !< attention: may be overwritten if lfront or lozpr is true
!!$  REAL(wp), POINTER :: lat_rmscon_eq !< [degN] rmscon_eq is used equatorward of this latitude
!!$  REAL(wp), POINTER :: lat_rmscon    !< [degN] rmscon is used poleward of this latitude
!!$  REAL(wp), POINTER :: rmscon_eq     !< [m/s]  rms constant used equatorward of lat_rmscon_eq
    !
    emiss_lev     => echam_gwd_config(jg)% emiss_lev
    rmscon        => echam_gwd_config(jg)% rmscon
    kstar         => echam_gwd_config(jg)% kstar
    !
!!$    lfront        => echam_gwd_config(jg)% lfront
!!$    !
!!$    lozpr         => echam_gwd_config(jg)% lozpr
!!$    pcrit         => echam_gwd_config(jg)% pcrit
!!$    pcons         => echam_gwd_config(jg)% pcons
!!$    !
!!$    lrmscon_lat   => echam_gwd_config(jg)% lrmscon_lat
!!$    lat_rmscon_eq => echam_gwd_config(jg)% lat_rmscon_eq
!!$    lat_rmscon    => echam_gwd_config(jg)% lat_rmscon
!!$    rmscon_eq     => echam_gwd_config(jg)% rmscon_eq

    dissip_gwd(:,:) = 0.0_wp
    tend_u_gwd(:,:) = 0.0_wp
    tend_v_gwd(:,:) = 0.0_wp

    !--  Check consistency of nc, jcs and jce
    !
    IF ( nc /= jce-jcs+1 ) CALL finish(TRIM(routine),'nc /= jce-jcs+1')

    !--  Initialize the ccc/mam hines gwd scheme
    !

    diffco(:,:) = 0.0_wp

    flux_u(:,:) = 0.0_wp
    flux_v(:,:) = 0.0_wp

    uhs(:,:) = 0.0_wp
    vhs(:,:) = 0.0_wp

    ! Wind variances form orographic gravity waves
    ! Note: the code is NOT fully implemeted for this case!

    sigsqmcw(:,:,:) = 0.0_wp
    sigmatm(:,:)    = 0.0_wp

!!$    ! precipitation (check the units!):
!!$    zpcons = (1000.0_wp*86400.0_wp)/rhoh2o
!!$    zpr(1:nc)=zpcons*paprflux(1:nc)

    rgocp=rd/cpd

    ! Surface pressure:
    pressg(1:nc)=paphm1(jcs:jce,nlev+1)

    ! Vertical positioning arrays:
    DO jk=1,nlev
!IBM* novector
      DO jl=1,nc
        pressg_inv = 1._wp/pressg(jl)
        shj(jl,jk)=papm1(jl+jcs-1,jk)*pressg_inv
        sgj(jl,jk)=papm1(jl+jcs-1,jk)*pressg_inv
      END DO
      shxkj(1:nc,jk) = sgj(1:nc,jk)**rgocp
    END DO

!     sgj(1:nc,1:nlev)=shj(1:nc,1:nlev)

    !
    !     * calculate b v frequency at all points
    !     * and smooth bvfreq.

    DO jk=2,nlev
!IBM* novector
      DO jl=1,nc
        dttdsf(jl)=(ptm1(jl+jcs-1,jk)/shxkj(jl,jk)-ptm1(jl+jcs-1,jk-1)/shxkj(jl,jk-1)) &
          /(shj(jl,jk)-shj(jl,jk-1))
        dttdsf(jl)=MIN(dttdsf(jl), -5.0_wp/sgj(jl,jk))
        dttdsf(jl)=dttdsf(jl)*shxkj(jl,jk)

        bvfreq(jl,jk)=-dttdsf(jl)*sgj(jl,jk)/rd
      END DO
      bvfreq(1:nc,jk) = SQRT(bvfreq(1:nc,jk))

!!HW!! !IBM* novector
      bvfreq(1:nc,jk) = bvfreq(1:nc,jk)*grav/ptm1(jcs:jce,jk)
    END DO

    bvfreq(:,1) = bvfreq(:,2)

    DO jk=2,nlev
      DO jl=1,nc
        ratio=5.0_wp*LOG(sgj(jl,jk)/sgj(jl,jk-1))
        bvfreq(jl,jk) = (bvfreq(jl,jk-1) + ratio*bvfreq(jl,jk))/(1.0_wp+ratio)
      END DO
    END DO

    !     * altitude above ground, density and air mass.

    DO jk=1,nlev
       alt    (1:nc,jk) = pzh  (jcs:jce,jk+1) - pzh(jcs:jce,nlev+1)
       density(1:nc,jk) = prho (jcs:jce,jk)
       mair   (1:nc,jk) = pmair(jcs:jce,jk)
    END DO

    !
    !     * set molecular viscosity to a very small value.
    !     * if the model top is greater than 100 km then the actual
    !     * viscosity coefficient could be specified here.


    visc_mol(:,:) = 1.5e-5_wp


    ! use single value for azimuthal-dependent horizontal wavenumber:
    ! kstar = (old latitudinal dependence, introduce here if necessary)

    k_alpha(:,:) = kstar

    !     * defile bottom launch level (emission level of gws)

    levbot = nlev-emiss_lev

    !     * initialize switch for column calculation

    lorms(:) = .FALSE.

    !     * background wind minus value at bottom launch level.

    DO jk=1,levbot
      uhs(1:nc,jk) = pum1(jcs:jce,jk) - pum1(jcs:jce,levbot)
      vhs(1:nc,jk) = pvm1(jcs:jce,jk) - pvm1(jcs:jce,levbot)
    END DO

    !     * specify root mean square wind at bottom launch level.

    anis(1:nc,1:naz) = 1.0_wp/REAL(naz,wp)

!!$    IF (lrmscon_lat) THEN
!!$      ! latitude dependent gravity wave source
!!$      ! - poleward of lat_rmscon               : rmscon
!!$      ! - equatorward of lat_rmscon_eq         : rmscon_eq
!!$      ! - between lat_rmscon_eq and lat_rmscon : linear interpolation between rmscon and rmscon_eq
!!$      rmswind(1:nc) = ( MAX(MIN((ABS(lat_deg(1:nc))-lat_rmscon_eq),lat_rmscon-lat_rmscon_eq),0.0_wp) * rmscon      &
!!$        &              +MAX(MIN((lat_rmscon-ABS(lat_deg(1:nc)))   ,lat_rmscon-lat_rmscon_eq),0.0_wp) * rmscon_eq ) &
!!$        &            /(lat_rmscon-lat_rmscon_eq)
!!$    ELSE
      rmswind(1:nc) = rmscon
!!$    ENDIF

!!$       !     * gravity waves from fronts:
!!$       IF (lfront) THEN
!!$          CALL gw_fronts(krow, nc, nlev, nazmth, rmswind, anis)
!!$       ENDIF

!!$       !     * modulation by precipitation:
!!$       IF (lozpr) THEN
!!$          DO jl=1,nc
!!$             IF (zpr(jl) > pcrit) THEN
!!$                rmswind(jl) = rmscon + ( (zpr(jl)-pcrit)/zpr(jl) )*pcons
!!$             ENDIF
!!$          END DO
!!$       ENDIF

    DO jl=1,nc
      IF (rmswind(jl) > 0.0_wp) THEN
        lorms(jl) = .TRUE.
      ENDIF
    END DO

    !
    !     * calculate gw tendencies (note that diffusion coefficient and
    !     * heating rate only calculated if lheatcal = .TRUE.).
    !
    CALL hines_extro ( jg,                                        &
      &                nc, nlev, nazmth,                          &
      &                tend_u_gwd(jcs:jce,:),                     &
      &                tend_v_gwd(jcs:jce,:),                     &
      &                dissip_gwd(jcs:jce,:),                     &
      &                diffco,                                    &
      &                flux_u,                                    &
      &                flux_v,                                    &
      &                uhs, vhs,                                  &
      &                bvfreq,                                    &
      &                density,                                   &
      &                mair,                                      &
      &                visc_mol,                                  &
      &                alt,                                       &
      &                rmswind, anis, k_alpha, sigsqmcw,          &
      &                m_alpha,  mmin_alpha ,sigma_t, sigmatm,    &
      &                levbot, lorms)

  END SUBROUTINE gw_hines
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE hines_extro ( jg,                                        &
    &                      nlons, nlevs, nazmth,                      &
    &                      drag_u,                                    &
    &                      drag_v,                                    &
    &                      heat,                                      &
    &                      diffco,                                    &
    &                      flux_u,                                    &
    &                      flux_v,                                    &
    &                      vel_u,                                     &
    &                      vel_v,                                     &
    &                      bvfreq,                                    &
    &                      density,                                   &
    &                      mair,                                      &
    &                      visc_mol,                                  &
    &                      alt,                                       &
    &                      rmswind, anis, k_alpha, sigsqmcw,          &
    &                      m_alpha,  mmin_alpha, sigma_t, sigmatm,    &
    &                      lev2,  lorms)
    !
    !  main routine for hines' "extrowave" gravity wave parameterization based
    !  on hines' doppler spread theory. this routine calculates zonal
    !  and meridional components of gravity wave drag, heating rates
    !  and diffusion coefficient on a longitude by altitude grid.
    !  no "mythical" lower boundary region calculation is made.
    !
    !  aug. 13/95 - c. mclandress
    !  sept. /95  - n. mcfarlane
    !  1995- 2002 - e. manzini
    !
    !  modifications:
    !
    !  output arguements:
    !
    !     * drag_u = zonal component of gravity wave drag (m/s^2).
    !     * drag_v = meridional component of gravity wave drag (m/s^2).
    !     * heat   = gravity wave heating (J/kg/sec).
    !     * diffco = diffusion coefficient (m^2/sec)
    !     * flux_u = zonal component of vertical momentum flux (pascals)
    !     * flux_v = meridional component of vertical momentum flux (pascals)
    !
    !  input arguements:
    !
    !     * vel_u      = background zonal wind component (m/s).
    !     * vel_v      = background meridional wind component (m/s).
    !     * bvfreq     = background brunt vassala frequency (radians/sec).
    !     * density    = background density (kg/m^3)
    !     * mair       = background airmass (kg/m^2)
    !     * visc_mol   = molecular viscosity (m^2/s)
    !     * alt        = altitude of momentum, density, buoyancy levels (m)
    !     *              (note: levels ordered so that alt(i,1) > alt(i,2), etc.)
    !     * rmswind   = root mean square gravity wave wind at lowest level (m/s).
    !     * anis      = anisotropy factor (sum over azimuths is one)
    !     * lorms     = .true. for drag computation (column selector)
    !     * k_alpha    = horizontal wavenumber of each azimuth (1/m).
    !     * lev2       = index of last level (eg bottom) for drag calculation
    !     *              (i.e., lev1 < lev2 <= nlevs).
    !     * nlons      = number of longitudes.
    !     * nlevs      = number of vertical levels.
    !     * nazmth     = azimuthal array dimension (nazmth >= naz).
    !
    !  ouput diagnostics:
    !     * m_alpha      = cutoff vertical wavenumber (1/m).
    !     * mmin_alpha   = minimum value of cutoff wavenumber.
    !     * sigma_t      = total rms horizontal wind (m/s).

    !  work arrays:
    !     * v_alpha      = wind component at each azimuth (m/s) and if lheatcal=.TRUE.
    !     *                holds vertical derivative of cutoff wavenumber.
    !     * sigma_alpha  = total rms wind in each azimuth (m/s).
    !     * ak_alpha     = spectral amplitude factor at each azimuth
    !     *                (i.e.,{ajkj}) in m^4/s^2.
    !     * densb        = background density at bottom level.
    !     * bvfb         = buoyancy frequency at bottom level and
    !     *                work array for icutoff = 1.
    !     * losigma_t    =  .true. for total sigma not zero

    IMPLICIT NONE

    INTEGER  :: jg
    INTEGER  :: nlons, nlevs, nazmth, lev2

    REAL(wp) :: drag_u(nlons,nlevs),   drag_v(nlons,nlevs)
    REAL(wp) :: heat(nlons,nlevs),     diffco(nlons,nlevs)
    REAL(wp) :: flux_u(nlons,nlevs),   flux_v(nlons,nlevs)
    REAL(wp) :: flux(nlons,nlevs,nazmth)
    REAL(wp) :: vel_u(nlons,nlevs),    vel_v(nlons,nlevs)
    REAL(wp) :: bvfreq(nlons,nlevs),   density(nlons,nlevs), mair(nlons,nlevs)
    REAL(wp) :: visc_mol(nlons,nlevs), alt(nlons,nlevs)
    REAL(wp) :: rmswind(nlons),      bvfb(nlons),   densb(nlons)
    REAL(wp) :: anis(nlons,nazmth)
    REAL(wp) :: sigma_t(nlons,nlevs), sigsqmcw(nlons,nlevs,nazmth)
    REAL(wp) :: sigma_alpha(nlons,nlevs,nazmth), sigmatm(nlons,nlevs)

    REAL(wp) :: m_alpha(nlons,nlevs,nazmth), v_alpha(nlons,nlevs,nazmth)
    REAL(wp) :: ak_alpha(nlons,nazmth),      k_alpha(nlons,nazmth)
    REAL(wp) :: mmin_alpha(nlons,nazmth)
    REAL(wp) :: smoothr1(nlons,nlevs), smoothr2(nlons,nlevs)

    LOGICAL  :: lorms(nlons), losigma_t(nlons,nlevs)
    !
    !  internal variables.
    !
    INTEGER  :: i, n, l, lev1, il1, il2
!!$    INTEGER :: iprint

    ! Shortcuts to components of echam_gwd_config
    !
    LOGICAL , POINTER :: lheatcal      !< true : compute momentum flux dep., heating and diffusion coefficient
    !                                  !< false: compute only momentum flux deposition
    REAL(wp), POINTER :: kstar         !< [1/m] typical gravity wave horizontal wavenumber
    REAL(wp), POINTER :: m_min         !< [1/m] minimum bound in  vertical wavenumber
    !
    lheatcal      => echam_gwd_config(jg)% lheatcal
    kstar         => echam_gwd_config(jg)% kstar
    m_min         => echam_gwd_config(jg)% m_min

    !-----------------------------------------------------------------------
    !

    ! range of longitude index:
    il1 = 1
    il2 = nlons

    lev1=1              ! top level index

!!$    iprint = 0       !     * iprint     = 1 to print out various arrays.

    !
    !  buoyancy and density at bottom level.
    !
    bvfb(il1:il2)  = bvfreq(il1:il2,lev2)
    densb(il1:il2) = density(il1:il2,lev2)
    !
    !  initialize some variables
    !
    DO n = 1,naz
       DO l=lev1,lev2
          DO i=il1,il2
             m_alpha(i,l,n) =  m_min
          END DO
       END DO
    END DO
    !
    !  compute azimuthal wind components from zonal and meridional winds.
    !
    CALL hines_wind ( v_alpha,                                     &
      &               vel_u, vel_v, naz,                           &
      &               il1, il2, lev1, lev2, nlons, nlevs, nazmth )
    !
    !  calculate cutoff vertical wavenumber and velocity variances.
    !
    CALL hines_wavnum ( jg,                                        &
      &                 m_alpha, sigma_t, sigma_alpha, ak_alpha,   &
      &                 mmin_alpha, losigma_t,                     &
      &                 v_alpha, visc_mol, density, densb,         &
      &                 bvfreq, bvfb, rmswind, anis, lorms,        &
      &                 sigsqmcw, sigmatm,                         &
      &                 il1, il2, lev1, lev2, nlons, nlevs, nazmth)
    !
    !  smooth cutoff wavenumbers and total rms velocity in the vertical
    !  direction nsmax times, using flux_u as temporary work array.
    !
    IF (nsmax>0)  THEN
       DO n = 1,naz
          DO l=lev1,lev2
             DO i=il1,il2
                smoothr1(i,l) = m_alpha(i,l,n)
             END DO
          END DO
          CALL vert_smooth (smoothr1,smoothr2, smco, nsmax,      &
            &               il1, il2, lev1, lev2, nlons, nlevs )
          DO l=lev1,lev2
             DO i=il1,il2
                m_alpha(i,l,n) = smoothr1(i,l)
             END DO
          END DO
       END DO
       CALL vert_smooth ( sigma_t, smoothr2, smco, nsmax,      &
         &                il1, il2, lev1, lev2, nlons, nlevs )
    END IF
    !
    !  calculate zonal and meridional components of the
    !  momentum flux and drag.
    !
    CALL hines_flux ( flux_u, flux_v, flux, drag_u, drag_v,       &
      &               densb, mair,                                &
      &               m_alpha,  ak_alpha, k_alpha,                &
      &               m_min, naz,                                 &
      &               il1, il2, lev1, lev2, nlons, nlevs, nazmth, &
      &               lorms )
    !
    !  cutoff drag above alt_cutoff, using bvfb as temporary work array.
    !
    IF (icutoff==1)  THEN
       CALL hines_exp ( drag_u, bvfb, alt, alt_cutoff,  &
         &              il1, il2, lev1, lev2, nlons, nlevs )
       CALL hines_exp ( drag_v, bvfb, alt, alt_cutoff,  &
         &              il1, il2, lev1, lev2, nlons, nlevs )
    END IF
    !
!!$    !  print out various arrays for diagnostic purposes.
!!$    !
!!$    IF (iprint==1)  THEN
!!$       CALL hines_print ( flux_u, flux_v, drag_u, drag_v, alt,     &
!!$         &                sigma_t, sigma_alpha, v_alpha, m_alpha,  &
!!$         &                1, 1, il1, il2, lev1, lev2,              &
!!$         &                naz, nlons, nlevs, nazmth)
!!$    END IF
    !
    !  if not calculating heating rate and diffusion coefficient then finished.
    !
    !  heating rate and diffusion coefficient.

    IF (lheatcal) THEN
       !
       !
       CALL hines_heat ( heat, diffco,                    &
         &  bvfreq, mair,                                 &
         &  sigma_t, sigma_alpha,                         &
         &  flux, visc_mol, kstar, f1, f2, f3, f5, f6,    &
         &  naz, il1, il2, lev1, lev2, nlons, nlevs,      &
         &  nazmth, losigma_t )
    END IF

    !  finished.
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_extro

  SUBROUTINE hines_wavnum ( jg,                                          &
    &                       m_alpha, sigma_t, sigma_alpha, ak_alpha,     &
    &                       mmin_alpha, losigma_t,                       &
    &                       v_alpha, visc_mol, density, densb,           &
    &                       bvfreq, bvfb, rms_wind, anis, lorms,         &
    &                       sigsqmcw, sigmatm,                           &
    &                       il1, il2, levtop, levbot, nlons, nlevs, nazmth)
    !
    !  this routine calculates the cutoff vertical wavenumber and velocity
    !  variances on a longitude by altitude grid for the hines' doppler
    !  spread gravity wave drag parameterization scheme.
    !  note: (1) only values of four or eight can be used for # azimuths (naz).
    !        (2) only values of 1.0, 1.5 or 2.0 can be used for slope (slope).
    !        (3) if m_min not zero, only slope=1. can be used.
    !
    !  aug. 10/95 - c. mclandress
    !  2000-2001  - m. charron
    !  2002       - e. manzini
    !
    !  output arguements:
    !
    !     * m_alpha      = cutoff wavenumber at each azimuth (1/m).
    !     * sigma_t      = total rms horizontal wind (m/s).
    !     * sigma_alpha  = total rms wind in each azimuth (m/s).
    !     * ak_alpha     = spectral amplitude factor at each azimuth
    !     *                (i.e.,{ajkj}) in m^4/s^2.
    !     * losigma_t    =  .true. for total sigma not zero
    !     * mmin_alpha   = minimum value of cutoff wavenumber.


    !  input arguements:
    !
    !     * v_alpha  = wind component at each azimuth (m/s).
    !     * visc_mol = molecular viscosity (m^2/s)
    !     * density  = background density (kg/m^3).
    !     * densb    = background density at model bottom (kg/m^3).
    !     * bvfreq   = background brunt vassala frequency (radians/sec).
    !     * bvfb     = background brunt vassala frequency at model bottom.
    !     * rms_wind = root mean square gravity wave wind at lowest level (m/s).
    !     * anis      = anisotropy factor (sum over azimuths is one)
    !     * lorms       = .true. for drag computation at lowest level

    !     * levbot   = index of lowest vertical level.
    !     * levtop   = index of highest vertical level
    !     *            (note: if levtop < levbot then level index
    !     *             increases from top down).
    !     * il1      = first longitudinal index to use (il1 >= 1).
    !     * il2      = last longitudinal index to use (il1 <= il2 <= nlons).

    !     * nlons    = number of longitudes.
    !     * nlevs    = number of vertical levels.
    !     * nazmth   = azimuthal array dimension (nazmth >= naz).
    !

    !
    !  work arrays:
    !
    !     * i_alpha    = hines' integral at a single level.
    !
    !     * do_alpha = .true. for the azimuths and longitudes for
    !                   which to continue to compute the drag above
    !                   the lowest level


    IMPLICIT NONE

    INTEGER  :: jg
    INTEGER  :: il1, il2, levtop, levbot, nlons, nlevs, nazmth
    REAL(wp) :: m_alpha(nlons,nlevs,nazmth)
    REAL(wp) :: sigma_alpha(nlons,nlevs,nazmth)
    REAL(wp) :: sigalpmc(nlons,nlevs,nazmth)
    REAL(wp) :: sigsqh_alpha(nlons,nlevs,nazmth)
    REAL(wp) :: sigsqmcw(nlons,nlevs,nazmth)
    REAL(wp) :: sigma_t(nlons,nlevs)
    REAL(wp) :: sigmatm(nlons,nlevs)
    REAL(wp) :: ak_alpha(nlons,nazmth)
    REAL(wp) :: v_alpha(nlons,nlevs,nazmth)
    REAL(wp) :: visc_mol(nlons,nlevs)
!     REAL(wp) :: f2mod(nlons,nlevs)
    REAL(wp) :: density(nlons,nlevs),  densb(nlons)
    REAL(wp) :: bvfreq(nlons,nlevs),   bvfb(nlons),  rms_wind(nlons)
    REAL(wp) :: anis(nlons,nazmth)
    REAL(wp) :: i_alpha(nlons,nazmth), mmin_alpha(nlons,nazmth)

    LOGICAL  :: lorms(nlons), losigma_t(nlons,nlevs), do_alpha(nlons,nazmth)
    !
    ! internal variables.
    !
    REAL(wp) :: sqr_rms_wind(nlons)
    INTEGER  :: ilorms(nlons),ialpha(nlons)
    INTEGER  :: i, j, l, n, istart, lend, lincr, lbelow, nlorms, nalpha

    REAL(wp) :: m_sub_m_turb, m_sub_m_mol, m_trial, mmsq
    REAL(wp) :: visc, visc_min, f2mfac

    REAL(wp) :: n_over_m(nlons), sigfac(nlons), vtmp1(nlons), vtmp2(nlons)! , vtmp3(nlons), maxdiff

    ! Here for optimization purposes distinguish:
    ! - the standard slope = 1 -> slope+1=2
    INTEGER,  PARAMETER :: s1p1_par_integer = 2      ! for use as exponential
    REAL(wp), PARAMETER :: s1p1_par_real    = 2.0_wp ! for use as factor
    ! - the generic  slope = 1., 1.5 or 2.
    REAL(wp) :: sp1_var_real                         ! for use as exponential and factor

    CHARACTER(len=*), PARAMETER :: routine = 'mo_gw_hines:hines_wavnum'

    !-----------------------------------------------------------------------
    !
    ! Shortcuts to components of echam_gwd_config
    !
    REAL(wp), POINTER :: kstar         !< [1/m] typical gravity wave horizontal wavenumber
    REAL(wp), POINTER :: m_min         !< [1/m] minimum bound in  vertical wavenumber
    !
    kstar => echam_gwd_config(jg)% kstar
    m_min => echam_gwd_config(jg)% m_min

    visc_min = 1.e-10_wp

    sp1_var_real = slope + 1.0_wp

    mmsq = m_min**2

    !
    !  indices of levels to process.
    !
    IF (levbot > levtop)  THEN
       istart = levbot - 1
       lend   = levtop
       lincr  = -1
    ELSE
       CALL finish(TRIM(routine),'level index not increasing downward')
    END IF

    !   initialize logical flags and arrays
    DO l=1,nlevs
       losigma_t(:,l) = lorms(:)
    ENDDO
    DO n=1,nazmth
       do_alpha(:,n) = lorms(:)
    ENDDO

    nlorms = 0
    DO i = il1,il2
       IF (lorms(i)) THEN
          nlorms = nlorms + 1
          ilorms(nlorms) = i
       END IF
    END DO

    sigsqh_alpha(:,:,:) = 0.0_wp ! mpuetz: is this really necessary -> 15% of wavnum
    i_alpha(:,:) = 0.0_wp

    !
    ! calculate azimuthal variances at bottom level using anisotropy factor
    !
    sqr_rms_wind(il1:il2) = rms_wind(il1:il2)**2
    DO n = 1,naz
       DO i = il1,il2
          sigsqh_alpha(i,levbot,n) = anis(i,n)* sqr_rms_wind(i)
       END DO
    END DO
    !
    !  velocity variances at bottom level.
    !
    CALL hines_sigma ( sigma_t, sigma_alpha,          &
      &                sigsqh_alpha, naz, levbot,     &
      &                il1, il2, nlons, nlevs, nazmth)

    CALL hines_sigma ( sigmatm, sigalpmc,             &
      &                sigsqmcw, naz, levbot,         &
      &                il1, il2, nlons, nlevs, nazmth)
    !
    !  calculate cutoff wavenumber and spectral amplitude factor
    !  at bottom level where it is assumed that background winds vanish
    !  and also initialize minimum value of cutoff wavnumber.
    !

     IF ( ABS(slope-1.0_wp) < EPSILON(1.0_wp) ) THEN
       ! here slope=1 -> use parameters s1p1_par_integer for exponential and s1p1_par_real for factor
       DO n = 1,naz
!IBM* NOVECTOR
!IBM* ASSERT(NODEPS)
          DO j = 1,nlorms
             i = ilorms(j)
                m_alpha(i,levbot,n) = bvfb(i)/(f1*sigma_alpha(i,levbot,n) + f2*sigma_t(i,levbot))
                ak_alpha(i,n)   = s1p1_par_real * sigsqh_alpha(i,levbot,n)/(m_alpha(i,levbot,n)**s1p1_par_integer - mmsq)
                mmin_alpha(i,n) = m_alpha(i,levbot,n)
          END DO
       END DO
     ELSE
       ! here slope/=1 -> use variable real for exponential and factor
       DO n = 1,naz
!IBM* ASSERT(NODEPS)
          DO j = 1,nlorms
             i = ilorms(j)
             vtmp1(j) = m_alpha(i,levbot,n)
          END DO

          vtmp1(1:nlorms) = vtmp1(1:nlorms)**(-sp1_var_real)

!IBM* ASSERT(NODEPS)
          DO j = 1,nlorms
             i = ilorms(j)
                m_alpha(i,levbot,n) = bvfb(i)/(f1*sigma_alpha(i,levbot,n) + f2*sigma_t(i,levbot))
                ak_alpha(i,n)       = sp1_var_real * sigsqh_alpha(i,levbot,n) * vtmp1(j)
                mmin_alpha(i,n)     = m_alpha(i,levbot,n)
          END DO
       END DO
    ENDIF
    !
    !  calculate quantities from the bottom upwards,
    !  starting one level above bottom.
    !
    DO l = istart,lend,lincr
       !
       !  level beneath present level.
       !
       lbelow = l - lincr
       !
       !  calculate n/m_m where m_m is maximum permissible value of the vertical
       !  wavenumber (i.e., m > m_m are obliterated) and n is buoyancy frequency.
       !  m_m is taken as the smaller of the instability-induced
       !  wavenumber (m_sub_m_turb) and that imposed by molecular viscosity
       !  (m_sub_m_mol). since variance at this level is not yet known
       !  use value at level below.
       !
       nlorms = 0
       DO i = il1,il2
          IF (losigma_t(i,lbelow)) THEN
             nlorms = nlorms + 1
             ilorms(nlorms) = i
          END IF
       END DO

!IBM* ASSERT(NODEPS)
       DO j = 1,nlorms
          i = ilorms(j)

          f2mfac=sigmatm(i,lbelow)**2
!           f2mod(i,lbelow) =1.0_wp+ 2.0_wp*f2mfac / (f2mfac+sigma_t(i,lbelow)**2 )

          visc = MAX ( visc_mol(i,l), visc_min )
          n_over_m(i) = f2 * &
            & (1.0_wp+ 2.0_wp*f2mfac / (f2mfac+sigma_t(i,lbelow)**2 )) & !f2mod(i,lbelow)
            * sigma_t(i,lbelow)
          vtmp1(j) = bvfreq(i,l) / n_over_m(i)
          vtmp2(j) = bvfreq(i,l)*kstar/visc
       END DO

       CALL vec_cbrt(vtmp2, vtmp2, n=nlorms)

       DO j = 1,nlorms
          i = ilorms(j)
          m_sub_m_turb = vtmp1(j)
          m_sub_m_mol  = vtmp2(j)/f3
          IF (m_sub_m_turb>=m_sub_m_mol) n_over_m(i) = bvfreq(i,l) / m_sub_m_mol
       END DO
       !
       !  calculate cutoff wavenumber at this level.
       !
       DO n = 1,naz
          ! set the default
          m_alpha(il1:il2,l,n) = m_min

          nalpha = 0
          DO j = 1,nlorms
             i = ilorms(j)
             IF (do_alpha(i,n)) THEN
                nalpha = nalpha + 1
                ialpha(nalpha) = i
             END IF
          END DO

!IBM* ASSERT(NODEPS)
          DO j = 1,nalpha
             i = ialpha(j)
             !
             !  calculate trial value (variance at this level is not yet known:
             !  use value at level below). if trial value negative or larger
             !  minimum value (not permitted) then set it to minimum value.
             !
!IBM* mem_delay(v_alpha(i,l,n),20)
!IBM* mem_delay(sigma_alpha(i,lbelow,n),20)
!IBM* mem_delay(sigalpmc(i,lbelow,n),20)
             m_trial = bvfb(i) / ( f1 * ( sigma_alpha(i,lbelow,n) + sigalpmc(i,lbelow,n)) &
               &       + n_over_m(i) + v_alpha(i,l,n) )

             m_trial = MIN(m_trial,mmin_alpha(i,n))
             IF (m_trial <= 0.0_wp) m_trial = mmin_alpha(i,n)

             m_alpha(i,l,n) = m_trial

             !  do not permit cutoff wavenumber to be less than minimum  value.

             m_alpha(i,l,n) = MAX(m_alpha(i,l,n),m_min)
             !
             !  reset minimum value of cutoff wavenumber if necessary.
             !
             mmin_alpha(i,n) = MIN(mmin_alpha(i,n),m_alpha(i,l,n))
          END DO
       END DO
       !
       !  calculate the hines integral at this level.
       !
       CALL hines_intgrl ( i_alpha,                                     &
         &                 v_alpha, m_alpha, bvfb, m_min, naz,   &
         &                 l, il1, il2, nlons, nlevs, nazmth,           &
         &                 lorms, do_alpha )

       !
       !  calculate the velocity variances at this level.
       !
!IBM* novector
       DO i = il1,il2
          sigfac(i) = densb(i) / density(i,l) * bvfreq(i,l) / bvfb(i)
       END DO
       DO n = 1,naz
          DO i = il1,il2
             sigsqh_alpha(i,l,n) = sigfac(i) * ak_alpha(i,n) * i_alpha(i,n)
          END DO
       END DO
       CALL hines_sigma ( sigma_t, sigma_alpha, sigsqh_alpha, naz, l, &
         &                il1, il2, nlons, nlevs, nazmth )

       CALL hines_sigma ( sigmatm, sigalpmc, sigsqmcw, naz, l,   &
         &                il1, il2, nlons, nlevs, nazmth )

       !
       !  if total rms wind zero (no more drag) then set drag to false
       !
       DO i=il1,il2
          IF ( sigma_t(i,l) < EPSILON(1.0_wp)) THEN
             losigma_t(i,l) = .FALSE.
          ENDIF
       ENDDO
       !
       !  end of level loop.
       !
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_wavnum

  SUBROUTINE hines_wind (v_alpha,vel_u,vel_v,  &
    &                    naz,il1,il2,lev1,lev2,nlons,nlevs,nazmth)
    !
    !  this routine calculates the azimuthal horizontal background wind components
    !  on a longitude by altitude grid for the case of 4 or 8 azimuths for
    !  the hines' doppler spread gwd parameterization scheme.
    !
    !  aug. 7/95 - c. mclandress
    !
    !  output arguement:
    !
    !     * v_alpha   = background wind component at each azimuth (m/s).
    !     *             (note: first azimuth is in eastward direction
    !     *              and rotate in counterclockwise direction.)
    !
    !  input arguements:
    !
    !     * vel_u     = background zonal wind component (m/s).
    !     * vel_v     = background meridional wind component (m/s).
    !     * naz       = actual number of horizontal azimuths used (must be 4 or 8).
    !     * il1       = first longitudinal index to use (il1 >= 1).
    !     * il2       = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1      = first altitude level to use (lev1 >=1).
    !     * lev2      = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons     = number of longitudes.
    !     * nlevs     = number of vertical levels.
    !     * nazmth    = azimuthal array dimension (nazmth >= naz).
    !
    !  constants in data statements.
    !
    !     * umin  = minimum allowable value for zonal or meridional
    !     *         wind component (m/s).
    !
    !  subroutine arguements.

    INTEGER  :: naz, il1, il2, lev1, lev2
    INTEGER  :: nlons, nlevs, nazmth
    REAL(wp) :: v_alpha(nlons,nlevs,nazmth)
    REAL(wp) ::  vel_u(nlons,nlevs), vel_v(nlons,nlevs)
    !
    !  internal variables.
    !
    INTEGER  :: i, l
    REAL(wp) :: u, v, vmu, vpu, umin
    !-----------------------------------------------------------------------


    umin  = 0.001_wp

    SELECT CASE (naz)
       !
    CASE(4)  !  case with 4 azimuths.

       DO l = lev1,lev2
!CDIR NODEP
          DO i = il1,il2
             !
             u = vel_u(i,l)
             IF (ABS(u) < umin)  u = umin
             v = vel_v(i,l)
             IF (ABS(v) < umin)  v = umin
             !
             v_alpha(i,l,1) = u     ! east
             v_alpha(i,l,2) = v     ! north
             !
             v_alpha(i,l,3) = - u   ! west
             v_alpha(i,l,4) = - v   ! south
             !
          END DO
       END DO
       !
    CASE (8)   !  case with 8 azimuths.
       !
       DO l = lev1,lev2
!CDIR NODEP
          DO i = il1,il2
             !
             u = vel_u(i,l)
             IF (ABS(u) < umin)  u = SIGN(umin,u)
             v = vel_v(i,l)
             IF (ABS(v) < umin)  v = SIGN(umin,v)
             vpu = v + u
             IF (ABS(vpu) < umin)  vpu = SIGN(umin,vpu)
             vmu = v - u
             IF (ABS(vmu) < umin)  vmu = SIGN(umin,vmu)
             !
             v_alpha(i,l,1) = u                  ! east
             v_alpha(i,l,2) = cos45 * vpu        ! north east
             v_alpha(i,l,3) = v                  ! north
             v_alpha(i,l,4) = cos45 * vmu        ! north west
             !
             v_alpha(i,l,5) = - v_alpha(i,l,1)   ! west
             v_alpha(i,l,6) = - v_alpha(i,l,2)   ! south west
             v_alpha(i,l,7) = - v_alpha(i,l,3)   ! south
             v_alpha(i,l,8) = - v_alpha(i,l,4)   ! south east
             !
          END DO
       END DO
       !
    END SELECT
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_wind

  SUBROUTINE hines_flux ( flux_u, flux_v, flux, drag_u, drag_v,        &
    &                     densb, mair,                                 &
    &                     m_alpha, ak_alpha, k_alpha,                  &
    &                     m_min, naz,                           &
    &                     il1, il2, lev1, lev2, nlons, nlevs, nazmth,  &
    &                     lorms )
    !
    !  calculate zonal and meridional components of the vertical flux
    !  of horizontal momentum and corresponding wave drag (force per unit mass)
    !  on a longitude by altitude grid for the hines' doppler spread
    !  gwd parameterization scheme.
    !  note: only 4 or 8 azimuths can be used.
    !
    !  aug. 6/95 - c. mclandress
    !       2001 - m. charron
    !
    !  output arguements:
    !
    !     * flux_u = zonal component of vertical momentum flux (pascals)
    !     * flux_v = meridional component of vertical momentum flux (pascals)
    !     * drag_u = zonal component of drag (m/s^2).
    !     * drag_v = meridional component of drag (m/s^2).
    !
    !  input arguements:
    !
    !     * alt       = altitudes (m).
    !     * densb     = background density at bottom level (kg/m^3).
    !     * mair      = background airmass (kg/m^2).
    !     * m_alpha   = cutoff vertical wavenumber (1/m).
    !     * ak_alpha  = spectral amplitude factor (i.e., {ajkj} in m^4/s^2).
    !     * k_alpha   = horizontal wavenumber (1/m).
    !     * slope     = slope of incident vertical wavenumber spectrum.
    !     * m_min     = minimum allowable cutoff wavenumber (1/m)
    !     *             for spectral slope of one.
    !     * naz       = actual number of horizontal azimuths used (must be 4 or 8).
    !     * il1       = first longitudinal index to use (il1 >= 1).
    !     * il2       = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1      = first altitude level to use (lev1 >=1).
    !     * lev2      = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons     = number of longitudes.
    !     * nlevs     = number of vertical levels.
    !     * nazmth    = azimuthal array dimension (nazmth >= naz).
    !     * lorms     = .true. for drag computation (column selector)
    !
    !
    !  subroutine arguements.
    !

    INTEGER  :: naz, il1, il2, lev1, lev2, lev2p
    INTEGER  :: nlons, nlevs, nazmth
    REAL(wp) ::  m_min
    REAL(wp) ::  flux_u(nlons,nlevs), flux_v(nlons,nlevs)
    REAL(wp) ::  flux(nlons,nlevs,nazmth)
    REAL(wp) ::  drag_u(nlons,nlevs), drag_v(nlons,nlevs)
    REAL(wp) ::  densb(nlons), mair(nlons,nlevs)
    REAL(wp) ::  m_alpha(nlons,nlevs,nazmth)
    REAL(wp) ::  ak_alpha(nlons,nazmth), k_alpha(nlons,nazmth)

    LOGICAL  :: lorms(nlons)
    !
    !  internal variables.
    !
    REAL(wp) ::  ak_k_alpha(nlons,nazmth)
    INTEGER  :: i, l, lev1p, k

    REAL(wp) ::  inv_slope
    REAL(wp) ::  densb_slope(nlons)
    !-----------------------------------------------------------------------
    !
    lev1p = lev1 + 1
    lev2p = lev2 + 1
    DO k = 1, nazmth
      DO i = il1,il2
        ak_k_alpha(i,k) = ak_alpha(i,k)*k_alpha(i,k)
      END DO
    END DO
    !
    !  for case where slope = 1.:
    !
    IF ( ABS(slope-1.0_wp) < EPSILON(1.0_wp) )  THEN
       !
       !  case with 4 azimuths.
       !
       IF (naz==4)  THEN
          DO l = lev1,lev2
             DO i = il1,il2
                flux(i,l,:) = ak_k_alpha(i,:)*(m_alpha(i,l,:)-m_min)
                flux_u(i,l) = flux(i,l,1) - flux(i,l,3)
                flux_v(i,l) = flux(i,l,2) - flux(i,l,4)
             END DO
          END DO
       END IF
       !
       !  case with 8 azimuths.
       !
       IF (naz==8)  THEN
          DO l = lev1,lev2
             DO k = 1, nazmth
                DO i = il1,il2
                  flux(i,l,k) = ak_k_alpha(i,k) * (m_alpha(i,l,k)-m_min)
                END DO
             END DO
             DO i = il1,il2
                flux_u(i,l) = flux(i,l,1) - flux(i,l,5) &
                  &         + cos45*( flux(i,l,2) - flux(i,l,4) - flux(i,l,6) + flux(i,l,8) )
                flux_v(i,l) = flux(i,l,3) - flux(i,l,7) &
                  &         + cos45*( flux(i,l,2) + flux(i,l,4) - flux(i,l,6) - flux(i,l,8) )
             END DO
          END DO
       END IF

    !
    !  for case where slope /= 1.:
    !
    ELSE
       !
       !  case with 4 azimuths.
       !
       IF (naz==4)  THEN
          DO l = lev1,lev2
             DO i = il1,il2
                flux(i,l,:) = ak_k_alpha(i,:)*m_alpha(i,l,:)**slope
                flux_u(i,l) = flux(i,l,1) - flux(i,l,3)
                flux_v(i,l) = flux(i,l,2) - flux(i,l,4)
             END DO
          END DO
       END IF
       !
       !  case with 8 azimuths.
       !
       IF (naz==8)  THEN
          DO l = lev1,lev2
             DO k = 1, nazmth
                DO i = il1,il2
                  flux(i,l,k) = ak_k_alpha(i,k)*m_alpha(i,l,k)**slope
                END DO
             END DO
             DO i = il1,il2
                flux_u(i,l) = flux(i,l,1) - flux(i,l,5) &
                  &         + cos45*( flux(i,l,2) - flux(i,l,4) - flux(i,l,6) + flux(i,l,8) )
                flux_v(i,l) = flux(i,l,3) - flux(i,l,7) &
                  &         + cos45*( flux(i,l,2) + flux(i,l,4) - flux(i,l,6) - flux(i,l,8) )
             END DO
          END DO
       END IF

    END IF
    !
    !  calculate flux from sum.
    !
    !
    !  for case where slope = 1.:
    !
    IF ( ABS(slope-1.0_wp) < EPSILON(1.0_wp) )  THEN
       DO l = lev1,lev2
          DO i = il1,il2
             flux_u(i,l) = flux_u(i,l) * densb(i)
             flux_v(i,l) = flux_v(i,l) * densb(i)
          END DO
          DO k = 1, nazmth
             DO i = il1,il2
                flux(i,l,k) = flux(i,l,k) * densb(i)
             END DO
          END DO
       END DO

    !
    !  for case where slope /= 1.:
    !
    ELSE
       inv_slope = 1._wp / slope
       densb_slope(il1:il2) = densb(il1:il2) * inv_slope
       DO l = lev1,lev2
          DO i = il1,il2
             flux_u(i,l) = flux_u(i,l) * densb_slope(i)
             flux_v(i,l) = flux_v(i,l) * densb_slope(i)
          END DO
          DO k = 1, nazmth
             DO i = il1,il2
                flux(i,l,k) = flux(i,l,k) * densb_slope(i)
             END DO
          END DO
       END DO
    END IF
    !
    !  calculate drag at intermediate levels
    !
    DO l = lev1p,lev2
!IBM* NOVECTOR
       DO i = il1,il2
          IF (lorms(i)) THEN
             drag_u(i,l) = - ( flux_u(i,l-1) - flux_u(i,l) ) / mair(i,l)
             drag_v(i,l) = - ( flux_v(i,l-1) - flux_v(i,l) ) / mair(i,l)
          ENDIF
       END DO
    END DO
    !
    !  drag at first and last levels using one-side differences.
    !
!IBM* NOVECTOR
    DO i = il1,il2
       IF (lorms(i)) THEN
          drag_u(i,lev1) =  flux_u(i,lev1)  / mair(i,lev1)
          drag_v(i,lev1) =  flux_v(i,lev1)  / mair(i,lev1)
       ENDIF
    END DO
    IF (nlevs > lev2) THEN
!IBM* NOVECTOR
       DO i = il1,il2
          IF (lorms(i)) THEN
             drag_u(i,lev2p) = - flux_u(i,lev2)  / mair(i,lev2p)
             drag_v(i,lev2p) = - flux_v(i,lev2)  / mair(i,lev2p)
          ENDIF
       END DO
    ENDIF
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_flux

  SUBROUTINE hines_heat ( heat, diffco,                                 &
    &                     bvfreq, mair,                                 &
    &                     sigma_t, sigma_alpha,                         &
    &                     flux, visc_mol, kstar, f1, f2, f3, f5, f6,    &
    &                     naz, il1, il2, lev1, lev2, nlons, nlevs,      &
    &                     nazmth, losigma_t )
    !
    !  this routine calculates the gravity wave induced heating and
    !  diffusion coefficient on a longitude by altitude grid for
    !  the hines' doppler spread gravity wave drag parameterization scheme.
    !
    !  This routine can be used for nonzero minimum cutoff wavenumber (m_min)
    !  only in the case of spectral slope=1, in which case m_min is not needed
    !  since its vertical derivative is zero.
    !
    !  aug. 6/95 - c. mclandress
    !  2001      - m. charron
    !
    !  output arguements:
    !
    !     * heat   = gravity wave heating (k/sec).
    !     * diffco = diffusion coefficient (m^2/sec)
    !
    !  input arguements:
    !
    !
    !     * bvfreq      = background brunt vassala frequency (rad/sec).
    !     * mair        = background airmass (kg/m^2).
    !
    !     * sigma_t     = total rms horizontal wind (m/s).
    !
    !     * visc_mol    = molecular viscosity (m^2/s).
    !     * kstar       = typical gravity wave horizontal wavenumber (1/m).
    !     * slope       = slope of incident vertical wavenumber spectrum.
    !     * f1,f2,f3,f5,f6 = hines's fudge factors.
    !
    !     * il1         = first longitudinal index to use (il1 >= 1).
    !     * il2         = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1        = first altitude level to use (lev1 >=1).
    !     * lev2        = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons       = number of longitudes.
    !     * nlevs       = number of vertical levels.
    !     * nazmth      = azimuthal array dimension (nazmth >= naz).
    !     * losigma_t   = .true. for total sigma not zero
    !

    IMPLICIT NONE

    INTEGER  ::  naz, il1, il2, lev1, lev2, nlons, nlevs, nazmth
    REAL(wp) ::  kstar, f1, f2, f3, f5, f6
    REAL(wp) ::  heat(nlons,nlevs), diffco(nlons,nlevs)
    REAL(wp) ::  bvfreq(nlons,nlevs), mair(nlons,nlevs)
    REAL(wp) ::  sigma_t(nlons,nlevs),  sigma_alpha(nlons,nlevs,nazmth)
    REAL(wp) ::  flux(nlons,nlevs,nazmth), visc_mol(nlons,nlevs)
    LOGICAL  ::  losigma_t(nlons,nlevs)
    !
    ! internal variables.
    !
    INTEGER  :: i, l, n, lev1p
    REAL(wp) :: m_sub_m_turb, m_sub_m_mol, m_sub_m, heatng
    REAL(wp) :: visc, visc_min

    REAL(wp) :: dfdz(nlons,nlevs,nazmth)
    !-----------------------------------------------------------------------

    visc_min = 1.e-10_wp

    lev1p = lev1 + 1

    DO l = lev1p,lev2
       DO i = il1,il2
          IF (losigma_t(i,l)) THEN
             visc    = MAX ( visc_mol(i,l), visc_min )
             m_sub_m_turb = bvfreq(i,l) / ( f2 * sigma_t(i,l) )
             m_sub_m_mol  = (bvfreq(i,l)*kstar/visc)**0.33333333_wp/f3
             m_sub_m      = MIN ( m_sub_m_turb, m_sub_m_mol )
!CDIR UNROLL=8
             dfdz(i,l,:) = ( flux(i,l-1,:) - flux(i,l,:) ) / mair(i,l) &
               &         * ( f1*sigma_alpha(i,l,:) + bvfreq(i,l)/m_sub_m )
          ENDIF
       END DO
    END DO

    DO i = il1,il2
       IF (losigma_t(i,lev1)) THEN
          visc    = MAX ( visc_mol(i,lev1), visc_min )
          m_sub_m_turb = bvfreq(i,lev1) / ( f2 * sigma_t(i,lev1) )
          m_sub_m_mol  = (bvfreq(i,lev1)*kstar/visc)**0.33333333_wp/f3
          m_sub_m      = MIN ( m_sub_m_turb, m_sub_m_mol )
!CDIR UNROLL=8
          dfdz(i,lev1,:) = -flux(i,lev1,:) / mair(i,l) &
            &             * ( f1*sigma_alpha(i,lev1,:) + bvfreq(i,lev1)/m_sub_m )
       ENDIF
    END DO

    !
    !  heating and diffusion.

    !
    !  maximum permissible value of cutoff wavenumber is the smaller
    !  of the instability-induced wavenumber (m_sub_m_turb) and
    !  that imposed by molecular viscosity (m_sub_m_mol).
    !
    !
    DO l = lev1,lev2
       DO i = il1,il2
          IF (losigma_t(i,l)) THEN
             visc    = MAX ( visc_mol(i,l), visc_min )
             m_sub_m_turb = bvfreq(i,l) / ( f2 * sigma_t(i,l) )
             m_sub_m_mol  = (bvfreq(i,l)*kstar/visc)**0.33333333_wp/f3
             m_sub_m      = MIN ( m_sub_m_turb, m_sub_m_mol )
             heatng = 0.0_wp
             DO n=1,naz
                heatng = heatng - f5 * dfdz(i,l,n)
             ENDDO
             diffco(i,l) = f6 * heatng**0.33333333_wp / m_sub_m**1.33333333_wp
             heat(i,l)   = heatng
          ENDIF
       END DO
    END DO

    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_heat

  SUBROUTINE hines_sigma (sigma_t,sigma_alpha,sigsqh_alpha,  &
    &                     naz,lev,il1,il2,nlons,nlevs,nazmth)
    !
    !  this routine calculates the total rms and azimuthal rms horizontal
    !  velocities at a given level on a longitude by altitude grid for
    !  the hines' doppler spread gwd parameterization scheme.
    !  note: only four or eight azimuths can be used.
    !
    !  aug. 7/95 - c. mclandress
    !
    !  output arguements:
    !
    !     * sigma_t      = total rms horizontal wind (m/s).
    !     * sigma_alpha  = total rms wind in each azimuth (m/s).
    !
    !  input arguements:
    !
    !     * sigsqh_alpha = portion of wind variance from waves having wave
    !     *                normals in the alpha azimuth (m/s).
    !     * naz       = actual number of horizontal azimuths used (must be 4 or 8).
    !     * lev       = altitude level to process.
    !     * il1       = first longitudinal index to use (il1 >= 1).
    !     * il2       = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * nlons     = number of longitudes.
    !     * nlevs     = number of vertical levels.
    !     * nazmth    = azimuthal array dimension (nazmth >= naz).
    !
    !  subroutine arguements.
    !

    INTEGER  :: lev, naz, il1, il2
    INTEGER  :: nlons, nlevs, nazmth
    REAL(wp) ::  sigma_t(nlons,nlevs)
    REAL(wp) :: sigma_alpha(nlons,nlevs,nazmth)
    REAL(wp) :: sigsqh_alpha(nlons,nlevs,nazmth)
    !
    !  internal variables.
    !
    INTEGER  :: i, n
    REAL(wp) :: sum_even, sum_odd
    !-----------------------------------------------------------------------
    !
    !  calculate azimuthal rms velocity for the 4 azimuth case.
    !
    IF (naz==4)  THEN
!CDIR NODEP
       DO i = il1,il2
          sigma_alpha(i,lev,1) = sigsqh_alpha(i,lev,1)+sigsqh_alpha(i,lev,3)
          sigma_alpha(i,lev,2) = sigsqh_alpha(i,lev,2)+sigsqh_alpha(i,lev,4)
       END DO

       sigma_alpha(il1:il2,lev,1) =  SQRT(sigma_alpha(il1:il2,lev,1))
       sigma_alpha(il1:il2,lev,2) =  SQRT(sigma_alpha(il1:il2,lev,2))

       DO i = il1,il2
          sigma_alpha(i,lev,3) = sigma_alpha(i,lev,1)
          sigma_alpha(i,lev,4) = sigma_alpha(i,lev,2)
       END DO
    END IF
    !
    !  calculate azimuthal rms velocity for the 8 azimuth case.
    !
    IF (naz==8)  THEN
!CDIR NODEP
       DO i = il1,il2
          sum_odd  = ( sigsqh_alpha(i,lev,1) + sigsqh_alpha(i,lev,3)            &
            &        + sigsqh_alpha(i,lev,5) + sigsqh_alpha(i,lev,7) ) * 0.5_wp
          sum_even = ( sigsqh_alpha(i,lev,2) + sigsqh_alpha(i,lev,4)            &
            &        + sigsqh_alpha(i,lev,6) + sigsqh_alpha(i,lev,8) ) * 0.5_wp
          sigma_alpha(i,lev,1) = sigsqh_alpha(i,lev,1) + sigsqh_alpha(i,lev,5) + sum_even
          sigma_alpha(i,lev,2) = sigsqh_alpha(i,lev,2) + sigsqh_alpha(i,lev,6) + sum_odd
          sigma_alpha(i,lev,3) = sigsqh_alpha(i,lev,3) + sigsqh_alpha(i,lev,7) + sum_even
          sigma_alpha(i,lev,4) = sigsqh_alpha(i,lev,4) + sigsqh_alpha(i,lev,8) + sum_odd
       END DO

       sigma_alpha(il1:il2,lev,1) = SQRT(sigma_alpha(il1:il2,lev,1))
       sigma_alpha(il1:il2,lev,2) = SQRT(sigma_alpha(il1:il2,lev,2))
       sigma_alpha(il1:il2,lev,3) = SQRT(sigma_alpha(il1:il2,lev,3))
       sigma_alpha(il1:il2,lev,4) = SQRT(sigma_alpha(il1:il2,lev,4))

       DO i = il1,il2
          sigma_alpha(i,lev,5) = sigma_alpha(i,lev,1)
          sigma_alpha(i,lev,6) = sigma_alpha(i,lev,2)
          sigma_alpha(i,lev,7) = sigma_alpha(i,lev,3)
          sigma_alpha(i,lev,8) = sigma_alpha(i,lev,4)
       END DO
    END IF
    !
    !  calculate total rms velocity.
    !
    DO i = il1,il2
       sigma_t(i,lev) = 0.0_wp
    END DO
    DO n = 1,naz
       DO i = il1,il2
          sigma_t(i,lev) = sigma_t(i,lev) + sigsqh_alpha(i,lev,n)
       END DO
    END DO
    sigma_t(il1:il2,lev) = SQRT(sigma_t(il1:il2,lev))
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_sigma

  SUBROUTINE hines_intgrl (i_alpha,                                     &
    &                      v_alpha, m_alpha, bvfb, m_min, naz,   &
    &                      lev, il1, il2, nlons, nlevs, nazmth,         &
    &                      lorms, do_alpha)
    !
    !  this routine calculates the vertical wavenumber integral
    !  for a single vertical level at each azimuth on a longitude grid
    !  for the hines' doppler spread gwd parameterization scheme.
    !  note: (1) only spectral slopes of 1, 1.5 or 2 are permitted.
    !        (2) the integral is written in terms of the product qm
    !            which by construction is always less than 1. series
    !            solutions are used for small |qm| and analytical solutions
    !            for remaining values.
    !
    !  aug. 8/95 - c. mclandress
    !  2001      - m. charron
    !  2003      - l. kornblueh
    !
    !  output arguement:
    !
    !     * i_alpha = hines' integral.
    !
    !  input arguements:
    !
    !     * v_alpha = azimuthal wind component (m/s).
    !     * m_alpha = azimuthal cutoff vertical wavenumber (1/m).
    !     * bvfb    = background brunt vassala frequency at model bottom.
    !     * m_min   = minimum allowable cutoff vertical wavenumber (1/m)
    !     *           for spectral slope of one.
    !     * slope   = slope of initial vertical wavenumber spectrum
    !     *           (must use slope = 1., 1.5 or 2.)
    !     * naz     = actual number of horizontal azimuths used.
    !     * lev     = altitude level to process.
    !     * il1     = first longitudinal index to use (il1 >= 1).
    !     * il2     = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * nlons   = number of longitudes.
    !     * nlevs   = number of vertical levels.
    !     * nazmth  = azimuthal array dimension (nazmth >= naz).
    !     * lorms     = .true. for drag computation (column selector)

    !
    !  constants in data statements:
    !
    !     * qmin = minimum value of q_alpha (avoids indeterminant form of integral)
    !     * qm_min = minimum value of q_alpha * m_alpha (used to avoid numerical
    !     *          problems).
    !

    IMPLICIT NONE

    INTEGER  :: lev, naz, il1, il2, nlons, nlevs, nazmth
    REAL(wp) :: i_alpha(nlons,nazmth)
    REAL(wp) :: v_alpha(nlons,nlevs,nazmth)
    REAL(wp) :: m_alpha(nlons,nlevs,nazmth)
    REAL(wp) :: bvfb(nlons), rbvfb(nlons), m_min

    LOGICAL  :: lorms(nlons), do_alpha(nlons,nazmth)
    !
    !  internal variables.
    !
    INTEGER  :: i, n, j
    REAL(wp) :: q_alpha, qm, qmm, sqrtqm, q_min, qm_min, ztmp

    !  variables for sparse vector optimization
    INTEGER  :: ixi(nlons*nazmth), ixnaz(nlons*nazmth)
    INTEGER  :: itx(nlons*nazmth)
    INTEGER  :: ilorms(nlons),icond(nlons)

    REAL(wp) :: vinp (nlons*nazmth), vout (nlons*nazmth)
    REAL(wp) :: winp (nlons*nazmth), wout (nlons*nazmth)
    REAL(wp) :: xinp (nlons*nazmth), xout (nlons*nazmth)
    INTEGER  :: vlen,nlorms,ic,ix,nerror

    CHARACTER(len=*), PARAMETER :: routine = 'mo_gw_hines:hines_hines_intgrl'

    !-----------------------------------------------------------------------
    !
    !  initialize local scalar and arrays

    q_min = 1.0_wp
    qm_min = 0.01_wp

!IBM *NOVECTOR
    DO i = il1,il2
      rbvfb(i)=1.0_wp/bvfb(i)
    ENDDO
    nlorms = 0
    DO i = il1,il2
       IF (lorms(i)) THEN
          nlorms = nlorms + 1
          ilorms(nlorms) = i
       END IF
    END DO
    !
    !  for integer value slope = 1.
    !
    IF ( ABS(slope-1.0_wp) < EPSILON(1.0_wp) )  THEN
       DO n = 1,naz
!IBM* ASSERT(NODEPS)
          DO j = 1,nlorms
             i = ilorms(j)
             icond(j) = INT(FSEL(m_min - m_alpha(i,lev,n),0._wp,1._wp))
          END DO

          vlen = 0
          DO j = 1,nlorms
             i = ilorms(j)
             IF (icond(j)>0) THEN
                vlen = vlen + 1
                ixi(vlen) = i
             ELSE
                i_alpha(i,n) = 0.0_wp
                do_alpha(i,n) = .FALSE.
             END IF
          END DO

!IBM* ASSERT(NODEPS)
          DO j = 1,vlen
             i = ixi(j)

             q_alpha = v_alpha(i,lev,n) * rbvfb(i)
             vinp(j) = q_alpha**2
             qm      = q_alpha * m_alpha(i,lev,n)
             winp(j) = 1.0_wp-qm
             qmm     = q_alpha * m_min
             xinp(j) = 1.0_wp-qmm
!             IF (ABS(q_alpha) < q_min .OR. ABS(qm)< qm_min)) ic = ic + 1
             ztmp     = FSEL(ABS(qm)-qm_min  ,0._wp,1._wp)
             ztmp     = FSEL(vinp(j)-q_min**2,ztmp ,1._wp)
             icond(j) = INT(ztmp)
          END DO

          ! if |qm| is small (rare event) we must do a Taylor expansion

          ic = 0
          IF (vlen>0) THEN
             DO j = 1,vlen
                IF (icond(j)<=0) CYCLE
                ic = ic + 1
                itx(ic) = ixi(j)
             END DO
          END IF

          vout(1:vlen) = 1._wp/vinp(1:vlen)
          wout(1:vlen) = LOG(winp(1:vlen))
          if (vlen > 0) xout(1:vlen) = LOG(xinp(1:vlen))   ! don't fuse vectorized LOGs

          DO j = 1,vlen
             i = ixi(j)

             i_alpha(i,n) = - (wout(j) - winp(j) - xout(j) + xinp(j) ) * vout(j)
             !
             !  If i_alpha negative due to round off error, set it to zero
             !
             i_alpha(i,n) = MAX( i_alpha(i,n) , 0.0_wp )
          END DO

          ! do Taylor expansion for small |qm|
!IBM* ASSERT(NODEPS,ITERCNT(1))
          DO j = 1,ic
             i = itx(j)
             q_alpha = v_alpha(i,lev,n) * rbvfb(i)
             qm      = q_alpha * m_alpha(i,lev,n)
             qmm     = q_alpha * m_min
             IF ( ABS(q_alpha) < EPSILON(1.0_wp) )  THEN
                i_alpha(i,n) = ( m_alpha(i,lev,n)**2  - m_min**2 ) * 0.5_wp
             ELSE
!                 i_alpha(i,n) = ( qm**2 * 0.50_wp + qm**3 * one_third   &
!                   &            + qm**4 * 0.25_wp + qm**5 * 0.2_wp   &
!                   &            -qmm**2 * 0.50_wp -qmm**3  * one_third    &
!                   &            -qmm**4 * 0.25_wp -qmm**5 * 0.2_wp ) &
!                   &           / q_alpha**2
! LL Aply Horner scheme
                i_alpha(i,n) = ( qm**2 * ( 0.50_wp + qm * ( one_third +  &
                  &     qm * ( 0.25_wp + qm * 0.2_wp)))   &
                  &  - (qmm**2 * (0.50_wp  + qmm * ( one_third +  &
                  &     qmm * ( 0.25_wp  + qmm * 0.2_wp ))))) &
                  &           / q_alpha**2
             END IF
             i_alpha(i,n) = MAX( i_alpha(i,n) , 0.0_wp )
          END DO

       END DO ! n
    END IF
    !
    !  for integer value slope = 2.
    !
    IF ( ABS(slope-2.0_wp) < EPSILON(1.0_wp) )  THEN
       DO n = 1,naz
!IBM* ASSERT(NODEPS)
          DO j = 1,nlorms
             i = ilorms(j)
             q_alpha = v_alpha(i,lev,n) * rbvfb(i)
             vinp(j) = q_alpha**3
             qm = q_alpha * m_alpha(i,lev,n)
             winp(j) = 1.0_wp-qm
!            IF (ABS(q_alpha) < q_min .OR. ABS(qm) < qm_min) ic = ic + 1
             ztmp = FSEL(ABS(qm)-qm_min,0._wp,1._wp)
             ztmp = FSEL(ABS(q_alpha)-q_min,ztmp,1._wp)
             icond(j) = INT(ztmp)
          END DO

          ! if |qm| is small (rare event) we must do a Taylor expansion

          ic = 0
          IF (nlorms>0) THEN
             DO j = 1,nlorms
                IF (icond(j)<=0) CYCLE
                ic = ic + 1
                itx(ic) = ilorms(j)
             END DO
          END IF

          vout(1:nlorms) = 1._wp/vinp(1:nlorms)
          wout(1:nlorms) = LOG(winp(1:nlorms))

!IBM* ASSERT(NODEPS)
          DO j = 1,nlorms
             i = ilorms(j)
             i_alpha(i,n) = - ( wout(j) + 1.0_wp - winp(j) + winp(j)**2 * 0.5_wp  ) * vout(j)
          END DO

          ! do Taylor expansion for small |qm|
!IBM* ASSERT(NODEPS,ITERCNT(1))
          DO j = 1,ic
             i = itx(j)
             q_alpha = v_alpha(i,lev,n) * rbvfb(i)
             qm = q_alpha * m_alpha(i,lev,n)
             IF ( ABS(q_alpha) < EPSILON(1.0_wp) )  THEN
                i_alpha(i,n) = m_alpha(i,lev,n)**3 * one_third
             ELSE
                i_alpha(i,n) = ( qm**3/3._wp + qm**4/4._wp + qm**5/5._wp + qm**6/6._wp ) &
                  &           / q_alpha**3
             END IF
          END DO
       END DO
    END IF
    !
    !  for real value slope = 1.5
    !
    IF ( ABS(slope-1.5_wp) < EPSILON(1.0_wp) )  THEN
       ic = 0
       DO n = 1,naz
          DO j = 1,nlorms
             i = ilorms(j)
             !
             q_alpha = v_alpha(i,lev,n) * rbvfb(i)
             qm = q_alpha * m_alpha(i,lev,n)
             !
             !  if |qm| is small then use first 4 terms series of taylor
             !  series expansion of integral in order to avoid
             !  indeterminate form of integral,
             !  otherwise use analytical form of integral.
             !
             IF (ABS(q_alpha) < q_min .OR. ABS(qm) < qm_min)  THEN
                ! taylor series expansion is a very rare event.
                ! do sparse processing separately
                ic = ic+1
                ixi(ic) = i
                ixnaz(ic) = n
             ELSE
                qm     = ABS(qm)
                sqrtqm = SQRT(qm)
                IF (q_alpha >= 0.0_wp)  THEN
                   i_alpha(i,n) = ( LOG( (1.0_wp+sqrtqm)/(1.0_wp-sqrtqm) )  &
                     &             -2.0_wp*sqrtqm*(1.0_wp+qm/3.0_wp) )      &
                     &           / q_alpha**2.5_wp
                ELSE
                   i_alpha(i,n) = 2.0_wp * ( ATAN(sqrtqm) + sqrtqm*(qm/3.0_wp-1.0_wp) ) &
                     &           / ABS(q_alpha)**2.5_wp
                ENDIF
             ENDIF
             !
          END DO
       END DO
       ! taylor series expansion is a very rare event.
       ! do sparse processing here separately
!CDIR NODEP
       DO ix = 1, ic
          n = ixnaz(ix)
          i = ixi(ix)
          q_alpha = v_alpha(i,lev,n) * rbvfb(i)
          qm = q_alpha * m_alpha(i,lev,n)
          IF ( ABS(q_alpha) < EPSILON(1.0_wp) )  THEN
             i_alpha(i,n) = m_alpha(i,lev,n)**2.5_wp / 2.5_wp
          ELSE
             i_alpha(i,n) = ( qm/2.5_wp + qm**2/3.5_wp + qm**3/4.5_wp + qm**4/5.5_wp )   &
               &            * m_alpha(i,lev,n)**1.5_wp / q_alpha
          END IF
       ENDDO
    END IF
    !
    !  if integral is negative (which in principle should not happen) then
    !  print a message and some info since execution will abort when calculating
    !  the variances.
    !
    DO n = 1,naz
       DO i = il1, il2
          icond(i-il1+1) = INT(FSEL(i_alpha(i,n),0._wp,1._wp))
       END DO
       IF (il2 >= il1) THEN  ! mpuetz: don't fuse with previous loop
          nerror = 0
          DO i = il1, il2
             nerror = nerror + icond(i-il1+1)
          END DO
       END IF

       IF (nerror > 0) THEN

          ! find the first error and dump

          DO i = il1, il2
             IF (icond(i-il1+1) > 0) EXIT
          END DO

          WRITE (message_text,*) ' '
          CALL message(TRIM(routine),message_text)
          WRITE (message_text,*) '******************************'
          CALL message(TRIM(routine),message_text)
          WRITE (message_text,*) 'hines integral i_alpha < 0 '
          CALL message(TRIM(routine),message_text)
          WRITE (message_text,*) '  column i =',i
          CALL message(TRIM(routine),message_text)
          WRITE (message_text,*) '  azimuth n=',n
          CALL message(TRIM(routine),message_text)
          WRITE (message_text,*) '  levellev =',lev
          CALL message(TRIM(routine),message_text)
          WRITE (message_text,*) '  i_alpha  =',i_alpha(i,n)
          CALL message(TRIM(routine),message_text)
          WRITE (message_text,*) '  v_alpha  =',v_alpha(i,lev,n)
          CALL message(TRIM(routine),message_text)
          WRITE (message_text,*) '  m_alpha  =',m_alpha(i,lev,n)
          CALL message(TRIM(routine),message_text)
          WRITE (message_text,*) '  q_alpha  =',v_alpha(i,lev,n)*rbvfb(i)
          CALL message(TRIM(routine),message_text)
          WRITE (message_text,*) '  qm       =',v_alpha(i,lev,n)*rbvfb(i)*m_alpha(i,lev,n)
          CALL message(TRIM(routine),message_text)
          WRITE (message_text,*) '******************************'

          CALL finish(TRIM(routine),'Hines i_alpha integral is negative')

       END IF

    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_intgrl

!!$  SUBROUTINE hines_print (flux_u, flux_v, drag_u, drag_v, alt, sigma_t,    &
!!$    &                     sigma_alpha, v_alpha, m_alpha,                   &
!!$    &                     iu_print, iv_print,                              &
!!$    &                     ilprt1, ilprt2, levprt1, levprt2, naz,           &
!!$    &                     nlons, nlevs, nazmth)
!!$    !
!!$    !  print out altitude profiles of various quantities from
!!$    !  hines' doppler spread gravity wave drag parameterization scheme.
!!$    !  (note: only for naz = 4 or 8).
!!$    !
!!$    !  aug. 8/95 - c. mclandress
!!$    !
!!$    !  input arguements:
!!$    !
!!$    !     * iu_print = 1 to print out values in east-west direction.
!!$    !     * iv_print = 1 to print out values in north-south direction.
!!$
!!$    !     * ilprt1   = first longitudinal index to print.
!!$    !     * ilprt2   = last longitudinal index to print.
!!$    !     * levprt1  = first altitude level to print.
!!$    !     * levprt2  = last altitude level to print.
!!$    !
!!$
!!$    INTEGER  :: naz, ilprt1, ilprt2, levprt1, levprt2
!!$    INTEGER  :: nlons, nlevs, nazmth
!!$    INTEGER  :: iu_print, iv_print
!!$    REAL(wp) :: flux_u(nlons,nlevs), flux_v(nlons,nlevs)
!!$    REAL(wp) :: drag_u(nlons,nlevs), drag_v(nlons,nlevs)
!!$    REAL(wp) :: alt(nlons,nlevs), sigma_t(nlons,nlevs)
!!$    REAL(wp) :: sigma_alpha(nlons,nlevs,nazmth)
!!$    REAL(wp) :: v_alpha(nlons,nlevs,nazmth), m_alpha(nlons,nlevs,nazmth)
!!$    !
!!$    !  internal variables.
!!$    !
!!$    INTEGER  :: n_east, n_west, n_north, n_south
!!$    INTEGER  ::  i, l
!!$    !-----------------------------------------------------------------------
!!$    !
!!$    !  azimuthal indices of cardinal directions.
!!$    !
!!$    n_east = 1
!!$    IF (naz==4)  THEN
!!$       n_west  = 3
!!$       n_north = 2
!!$       n_south = 4
!!$    ELSE IF (naz==8)  THEN
!!$       n_west  = 5
!!$       n_north = 3
!!$       n_south = 7
!!$    END IF
!!$    !
!!$    !  print out values for range of longitudes.
!!$    !
!!$    DO i = ilprt1,ilprt2
!!$       !
!!$       !  print east-west wind, sigmas, cutoff wavenumbers, flux and drag.
!!$       !
!!$       IF (iu_print==1)  THEN
!!$          WRITE (nout,*)
!!$          WRITE (nout,'(a,i3)') 'hines gw (east-west) at longitude i =',i
!!$          WRITE (nout,6005)
!!$6005      FORMAT (15x,' u ',2x,'sig_e',2x,'sig_t',3x,'m_e', 4x,'m_w',4x,'fluxu',5x,'gwdu')
!!$          DO l = levprt1,levprt2
!!$             WRITE (nout,6701) alt(i,l)/1.e3_wp, v_alpha(i,l,n_east),   &
!!$               &               sigma_alpha(i,l,n_east), sigma_t(i,l),   &
!!$               &               m_alpha(i,l,n_east)*1.e3_wp,             &
!!$               &               m_alpha(i,l,n_west)*1.e3_wp,             &
!!$               &               flux_u(i,l)*1.e5_wp, drag_u(i,l)*24.0_wp*3600.0_wp
!!$          END DO
!!$6701      FORMAT (' z=',f7.2,1x,3f7.1,2f7.3,f9.4,f9.3)
!!$       END IF
!!$       !
!!$       !  print north-south winds, sigmas, cutoff wavenumbers, flux and drag.
!!$       !
!!$       IF (iv_print==1)  THEN
!!$          WRITE(nout,*)
!!$          WRITE(nout,'(a,i3)') 'hines gw (north-south) at longitude i =',i
!!$          WRITE(nout,6006)
!!$6006      FORMAT (15x,' v ',2x,'sig_n',2x,'sig_t',3x,'m_n',4x,'m_s',4x,'fluxv',5x,'gwdv')
!!$          DO l = levprt1,levprt2
!!$             WRITE (nout,6701) alt(i,l)/1.e3_wp, v_alpha(i,l,n_north),    &
!!$               &               sigma_alpha(i,l,n_north), sigma_t(i,l),    &
!!$               &               m_alpha(i,l,n_north)*1.e3_wp,              &
!!$               &               m_alpha(i,l,n_south)*1.e3_wp,              &
!!$               &               flux_v(i,l)*1.e5_wp, drag_v(i,l)*24.0_wp*3600.0_wp
!!$          END DO
!!$       END IF
!!$       !
!!$    END DO
!!$    !
!!$    !-----------------------------------------------------------------------
!!$  END SUBROUTINE hines_print

  SUBROUTINE hines_exp (darr, data_zmax, alt, alt_exp,     &
    &                   il1, il2, lev1, lev2, nlons, nlevs)
    !
    !  this routine exponentially damps a longitude by altitude array
    !  of darr above a specified altitude.
    !
    !  aug. 13/95 - c. mclandress
    !
    !  output arguements:
    !
    !     * darr = modified data array.
    !
    !  input arguements:
    !
    !     * darr    = original data array.
    !     * alt     = altitudes.
    !     * alt_exp = altitude above which exponential decay applied.

    !     * il1     = first longitudinal index to use (il1 >= 1).
    !     * il2     = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1    = first altitude level to use (lev1 >=1).
    !     * lev2    = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons   = number of longitudes.
    !     * nlevs   = number of vertical
    !
    !  input work arrays:
    !
    !     * data_zmax = data values just above altitude alt_exp.
    !

    INTEGER  :: il1, il2, lev1, lev2, nlons, nlevs
    REAL(wp) :: alt_exp
    REAL(wp) :: darr(nlons,nlevs), data_zmax(nlons), alt(nlons,nlevs)
    !
    ! internal variables.
    !
    REAL(wp) :: vtmp(nlons)
    INTEGER  :: ialt(nlons)
    INTEGER  :: levbot, levtop, lincr, i, l, nalt, j
    REAL(wp) :: hscale

    CHARACTER(len=*), PARAMETER :: routine = 'mo_gw_hines:hines_exp'

    !-----------------------------------------------------------------------

    hscale = 5.e3_wp

    !  index of lowest altitude level (bottom of drag calculation).
    !
    levbot = lev2
    levtop = lev1
    lincr  = 1
    IF (levbot > levtop)  THEN
       levbot = lev1
       levtop = lev2
       lincr  = -1
    ELSE
       CALL finish(TRIM(routine),'level index not increasing downward')
    END IF
    !
    !  data values at first level above alt_exp.
    !
    DO i = il1,il2
       DO l = levtop,levbot,lincr
          IF (alt(i,l) >= alt_exp)  THEN
             data_zmax(i) = darr(i,l)
          END IF
       END DO
    END DO
    !
    !  exponentially damp field above alt_exp to model top at l=1.
    !
    DO l = 1,lev2
       nalt = 0
       DO i = il1,il2
          IF (alt(i,l) >= alt_exp)  THEN
             nalt = nalt + 1
             ialt(nalt) = i
          END IF
       END DO
       DO j = 1,nalt
          i = ialt(j)
          vtmp(j) = (alt_exp-alt(i,l))/hscale
       END DO

       vtmp(1:nalt) = EXP(vtmp(1:nalt))

       DO j = 1,nalt
          i = ialt(j)
          darr(i,l) = data_zmax(i) * vtmp(j)
       END DO
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_exp

  SUBROUTINE vert_smooth (darr, work, coeff, nsmooth,         &
    &                     il1, il2, lev1, lev2, nlons, nlevs)
    !
    !  smooth a longitude by altitude array in the vertical over a
    !  specified number of levels using a three point smoother.
    !
    !  note: input array darr is modified on output!
    !
    !  aug. 3/95 - c. mclandress
    !
    !  output arguement:
    !
    !     * darr    = smoothed array (on output).
    !
    !  input arguements:
    !
    !     * darr    = unsmoothed array of data (on input).
    !     * work    = work array of same dimension as darr.
    !     * coeff   = smoothing coefficient for a 1:coeff:1 stencil.
    !     *           (e.g., coeff = 2 will result in a smoother which
    !     *           weights the level l gridpoint by two and the two
    !     *           adjecent levels (l+1 and l-1) by one).
    !     * nsmooth = number of times to smooth in vertical.
    !     *           (e.g., nsmooth=1 means smoothed only once,
    !     *           nsmooth=2 means smoothing repeated twice, etc.)
    !     * il1     = first longitudinal index to use (il1 >= 1).
    !     * il2     = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1    = first altitude level to use (lev1 >=1).
    !     * lev2    = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons   = number of longitudes.
    !     * nlevs   = number of vertical levels.
    !
    !  subroutine arguments.
    !
    IMPLICIT NONE

    INTEGER  :: nsmooth, il1, il2, lev1, lev2, nlons, nlevs
    REAL(wp) :: coeff
    REAL(wp) :: darr(nlons,nlevs), work(nlons,nlevs)
    !
    !  internal variables.
    !
    INTEGER  :: i, l, ns, lev1p, lev2m
    REAL(wp) :: sum_wts
    !-----------------------------------------------------------------------
    !
    !  calculate sum of weights.
    !
    sum_wts = coeff + 2.0_wp
    !
    lev1p = lev1 + 1
    lev2m = lev2 - 1
    !
    !  smooth nsmooth times
    !
    DO ns = 1,nsmooth
       !
       !  copy darr into work array.
       !
       DO l = lev1,lev2
          DO i = il1,il2
             work(i,l) = darr(i,l)
          END DO
       END DO
       !
       !  smooth array work in vertical direction and put into darr.
       !
       DO l = lev1p,lev2m
          DO i = il1,il2
             darr(i,l) = (work(i,l+1)+coeff*work(i,l)+work(i,l-1) ) / sum_wts
          END DO
       END DO
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vert_smooth

!!$  SUBROUTINE  gw_fronts(krow, nc, nlev, nazmth, rmswind, ani)
!!$    !
!!$    !
!!$    !  may 22/2000 - m. charron
!!$    !
!!$    !  output arguements:
!!$    !
!!$    !     * rmswind    = root mean square gravity wave wind at lowest level (m/s).
!!$    !     * ani        = anisotropy factor (sum over azimuths is one)
!!$    !
!!$    !  input arguments:
!!$    !
!!$    !     * klon     = number of longitudes
!!$    !     * nazmth   = azimuthal array dimension (nazmth >= naz).
!!$    !
!!$    !  subroutine arguments.
!!$    !
!!$    INTEGER,  INTENT(in)  :: krow, nc, nlev, nazmth
!!$    REAL(wp), INTENT(out) :: rmswind(nc)
!!$    REAL(wp), INTENT(out) :: ani(nc,nazmth)
!!$    !
!!$    !  internal variables.
!!$    !
!!$    INTEGER  :: jl, jdir
!!$    INTEGER  :: opp(nazmth)
!!$    REAL(wp) :: angle(nc), gen(nc)
!!$
!!$    ! Shortcuts to components of echam_gwd_config
!!$    !
!!$    REAL(wp), POINTER :: rmscon        !< [m/s] root mean square gravity wave wind at emission level
!!$    REAL(wp), POINTER :: rms_front     !< [m/s] rms frontal gw wind at source level
!!$    REAL(wp), POINTER :: front_thres   !< [(K/m)^2/hr] minimum value of the frontogenesis function,
!!$                                       !< for which gravity waves are emitted from fronts
!!$    !
!!$    rmscon      => echam_gwd_config(jg)% rmscon
!!$    rms_front   => echam_gwd_config(jg)% rms_front
!!$    front_thres => echam_gwd_config(jg)% front_thres
!!$
!!$    IF ( naz == 8 ) THEN
!!$       opp(1)=5
!!$       opp(2)=6
!!$       opp(3)=7
!!$       opp(4)=8
!!$       opp(5)=1
!!$       opp(6)=2
!!$       opp(7)=3
!!$       opp(8)=4
!!$    ELSE IF ( naz == 4 ) THEN
!!$       opp(1)=3
!!$       opp(2)=4
!!$       opp(3)=1
!!$       opp(4)=2
!!$    END IF
!!$
!!$    CALL calculate_gen ( krow, nc, nlev, gen, angle )
!!$
!!$    DO jl=1,nc
!!$       IF ( gen(jl) >= front_thres*2.7777778E-14_wp ) THEN
!!$          rmswind(jl) = rms_front
!!$          jdir=INT( angle(jl)/360.0_wp*REAL(naz,wp) ) + 1
!!$          ani(jl,    :    ) = 0.0_wp
!!$          ani(jl,    jdir ) = 0.5_wp
!!$          ani(jl,opp(jdir)) = 0.5_wp
!!$       ELSE
!!$          rmswind(jl) = rmscon
!!$          ani(jl,:)   = 1.0_wp/REAL(naz,wp)
!!$       END IF
!!$    END DO
!!$    !
!!$    !-----------------------------------------------------------------------
!!$  END SUBROUTINE gw_fronts
!!$
!!$  SUBROUTINE calculate_gen(krow, nc, nlev, gen, angle)
!!$    !
!!$    ! determine the value of the frontogenesis function
!!$    !
!!$    !  may 22/2000   - m. charron - first version
!!$    !  march 29/2001 - m. charron - reorganized as a module
!!$    !
!!$    !  output arguement:
!!$    !
!!$    !     * gen    = frontogenesis function in (K/m)^2/hour
!!$    !     * angle  = orientation (in degrees) of the gradiant of temp
!!$    !
!!$    !  input arguements:
!!$    !
!!$    !     * nc     = last longitudinal index to use (1 <= nc <= nlons).
!!$    !
!!$    !  subroutine arguements.
!!$    !
!!$    INTEGER,                 INTENT(in)  :: krow, nc, nlev
!!$    REAL(wp), DIMENSION(nc), INTENT(out) :: gen, angle
!!$    !
!!$    !  internal variables.
!!$    !
!!$    INTEGER                     :: jl, jglat, jrow, iplev
!!$    REAL(wp), DIMENSION(nc)     :: dalpsdl, dalpsdm, dudl, dvdl, ps  , vort, div
!!$    REAL(wp), DIMENSION(nc)     :: dtdl   , dtdm   , uu  , vv  , uup1, uum1, vvp1
!!$    REAL(wp), DIMENSION(nc)     :: vvm1   , tt     , ttp1, ttm1
!!$    REAL(wp), DIMENSION(nc)     :: mu, cstmu, zrcst
!!$    REAL(wp)                    :: za, zb, zap1, zbp1
!!$    REAL(wp)                    :: zaplus, zaminus, zbplus, zbminus
!!$    REAL(wp)                    :: term_1, term_2, term_3l, term_3m, term_4l, term_4m
!!$    REAL(wp)                    :: term_5, term_6, term_7 , term_8
!!$    REAL(wp), PARAMETER         :: pr=100000.0_wp
!!$    REAL(wp), PARAMETER         :: kappa=2.0_wp/7.0_wp
!!$
!!$    ! Shortcuts to components of echam_gwd_config
!!$    !
!!$    INTEGER , POINTER :: emiss_lev     !< number of levels above the ground at which gw are emitted
!!$    !
!!$    emiss_lev => echam_gwd_config(jg)% emiss_lev
!!$
!!$!------------------------------------------------------------------------------
!!$
!!$    jrow=krow
!!$    iplev=nlev-emiss_lev
!!$
!!$    ps     (:)=EXP(alpsm1(:,jrow))
!!$    za        =vct_a   (iplev)
!!$    zap1      =vct_a   (iplev+1)
!!$    zb        =vct_b   (iplev)
!!$    zbp1      =vct_b   (iplev+1)
!!$    zaplus    =zap1+za
!!$    zaminus   =zap1-za
!!$    zbplus    =zbp1+zb
!!$    zbminus   =zbp1-zb
!!$    DO jl=1,nc
!!$      jglat     =ilat(jl,jrow)
!!$      mu(jl)    =0.5_wp*gl_twomu(jglat)
!!$      zrcst(jl) =1.0_wp/gl_sqcst(jglat)
!!$    ENDDO
!!$    cstmu  (:)=1.0_wp-mu(:)*mu(:)
!!$    dtdl   (:)=re*cstmu(:)*dtlm1 (:,iplev,jrow)
!!$    dtdm   (:)=re*dtmm1 (:,iplev,jrow)
!!$    dalpsdl(:)=re*cstmu(:)*dalpslm1(:,jrow)
!!$    dalpsdm(:)=re*dalpsmm1(:,jrow)
!!$
!!$    uu  (:) = um1 (:,iplev  ,jrow)
!!$    vv  (:) = vm1 (:,iplev  ,jrow)
!!$    uup1(:) = um1 (:,iplev+1,jrow)
!!$    uum1(:) = um1 (:,iplev-1,jrow)
!!$    vvp1(:) = vm1 (:,iplev+1,jrow)
!!$    vvm1(:) = vm1 (:,iplev-1,jrow)
!!$    tt  (:) = tm1 (:,iplev  ,jrow)
!!$    ttp1(:) = tm1 (:,iplev+1,jrow)
!!$    ttm1(:) = tm1 (:,iplev-1,jrow)
!!$    vort(:) = vom1(:,iplev  ,jrow)
!!$    div (:) = dm1 (:,iplev  ,jrow)
!!$
!!$    dudl(:) = dudlm1(:,iplev  ,jrow)*zrcst(:)
!!$    dvdl(:) = dvdlm1(:,iplev  ,jrow)*zrcst(:)
!!$
!!$    DO jl=1,nc
!!$      term_1  = (2.0_wp*pr/(zaplus+ps(jl)*zbplus))**kappa
!!$      term_2  = 1.0_wp/(zaminus+ps(jl)*zbminus)
!!$      term_3l = ps(jl)*dalpsdl(jl)*zbplus*term_2
!!$      term_3m = ps(jl)*dalpsdm(jl)*zbplus*term_2
!!$      term_4l = ( dtdl(jl) - 0.25_wp*(ttp1(jl)-ttm1(jl)) * term_3l ) * term_1
!!$      term_4m = ( dtdm(jl) - 0.25_wp*(ttp1(jl)-ttm1(jl)) * term_3m ) * term_1
!!$      term_5  = 0.25_wp * (uup1(jl) - uum1(jl))
!!$      term_6  = 0.25_wp * (vvp1(jl) - vvm1(jl))
!!$      term_7  = 1.0_wp/(re*re*cstmu(jl))
!!$      term_8  = re * term_7
!!$      gen(jl) = -term_7*term_4l**2                                                  &
!!$        &       *(term_8*(dudl(jl)-term_5*term_3l)-mu(jl)*vv(jl)*term_8)            &
!!$        &       -cstmu(jl)/(re*re) * term_4m**2                                     &
!!$        &       *(1.0_wp/re*(re*div(jl)-dudl(jl)/cstmu(jl)-term_6*term_3m)          &
!!$        &       +mu(jl)*vv(jl)*term_8)-term_7* term_4l*term_4m                      &
!!$        &       *( 1.0_wp/re*(dvdl(jl)-term_6*term_3l)                              &
!!$        &       +cstmu(jl)/re*(dvdl(jl)/cstmu(jl)-re*vort(jl)-term_5*term_3m)       &
!!$        &       +2.0_wp*mu(jl)*uu(jl)/re)
!!$      angle(jl) = ATAN2(cstmu(jl)*dtdm(jl),dtdl(jl))*180.0_wp/pi+180.0_wp/REAL(naz,wp)
!!$      IF ( angle(jl) < 0.0_wp ) angle(jl) = angle(jl) + 360.0_wp
!!$    END DO
!!$    !
!!$    !-----------------------------------------------------------------------
!!$  END SUBROUTINE calculate_gen

END MODULE mo_gw_hines
