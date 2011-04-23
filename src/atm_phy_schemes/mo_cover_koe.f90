!>
!! Computation of cloud cover and grid mean cloud liquid water and cloud ice
!!
!! This routine takes information from turbulence, convection and grid-scale
!! to produce cloud properties used in radiation (and microphysics).
!!
!! Possible future options
!! - simple diagnostic (from turbulence, convection and grid scale)
!! - prognostic total water variance AND prognostic ice
!!
!!
!! @author Martin Koehler, DWD
!!
!!
!! @par Revision History
!! Developed by Martin Koehler  (starting 2010-08-23)
!! Modification by Martin Koehler, DWD (2010-11-20)
!! - options 0,1,2,3 ready
!!
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_cover_koe

  USE mo_kind,               ONLY: wp, i4

  USE mo_physical_constants, ONLY: re
  USE mo_run_nml,            ONLY: nproma, nlev
  USE mo_math_utilities,     ONLY: gamma_fct
  USE mo_math_constants,     ONLY: dbl_eps, pi

  USE mo_physical_constants, ONLY: rv     , & !! gas constant for water vapour
                                   rd     , & !! gas constant for dry air
                                   vtmpc1 , & !! r_v/r_d - 1
                                   o_m_rdv, & !! 1 - r_d/r_v
                                   rdv    , & !! r_d / r_v
                                   alv    , & !! latent heat of vapourization
                                   als    , & !! latent heat of sublimation
                                   alf    , & !! latent heat of fusion
                                   cpd    , & !! specific heat of dry air at constant press
                                   rcpd   , & !! (spec. heat of dry air at constant press)^-1
                                   tmelt  , & !! melting temperature of ice/snow
                                   rhoh2o , & !! density of liquid water (kg/m^3)
                                   grav   , & !! acceleration due to gravity
                                   tmelt      !! melting temperature of ice/snow

  USE mo_convect_tables,     ONLY: c1es   , & !! constants for computing the sat. vapour
                                   c3les  , & !! pressure over water (l) and ice (i)
                                   c3ies  , & !!               -- " --
                                   c4les  , & !!               -- " --
                                   c4ies  , & !!               -- " --
                                   c5les      !!               -- " --

  USE mo_cufunctions,        ONLY: foealfa    !! liquid fraction as in Tiedtke/Bechtold convection

  USE mo_exception,          ONLY: message, finish, message_text

  USE mo_cloud_diag,         ONLY: cloud_diag

  USE mo_cover_cosmo,        ONLY: cover_cosmo

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: cover_koe

!-------------------------------------------------------------------------

CONTAINS


!-------------------------------------------------------------------------
!
!  Cloud cover and cloud water/ice calculation.
!
!
!  @par Revision History
!  Initial version by Martin Koehler, DWD (2010-08-23)
!
!  Options:
!  (0) no clouds
!  (1) grid-scale cloud cover [1 or 0]
!  (2) clouds as in turbulence
!  (3) clouds as in COSMO
!  (4) diagnostic cloud cover          (in progress)
!  (5) prognostic total water variance (not started)
!  
!-------------------------------------------------------------------------

SUBROUTINE cover_koe( &
  & kidia, kfdia, klon, kstart, klev, & ! in:    dimensions (turn off physics above kstart)
  & icldscheme                      , & ! in:    cloud cover framework
  & tt                              , & ! in:    temperature (main levels)
  & pp                              , & ! in:    pressure (")
  & ps                              , & ! in:    surface pressure
  & pgeo                            , & ! in:    geopotential
  & rho                             , & ! in:    density
  & rcld                            , & ! inout: standard deviation of saturation deficit
  & ldland                          , & ! in:    land/sea mask
  & ldcum, kcbot, kctop             , & ! in:    convection: on/off, bottom, top
  & pmfude_rate                     , & ! in:    convection: updraft detrainment rate
  & plu                             , & ! in:    convection: updraft condensate 
  & qv, qc, qi                      , & ! inout: prognostic cloud variables
  & cc_tot, qv_tot, qc_tot, qi_tot )    ! out:   cloud output diagnostic


!! Subroutine arguments:
!! --------------------

INTEGER(KIND=i4), INTENT(IN) ::  &
  & kidia            , & ! horizontal start index
  & kfdia            , & ! horizontal end   index
  & klon             , & ! horizontal dimension
  & kstart           , & ! vertical start index (turn off physics above)
  & klev                 ! vertical dimension

INTEGER(KIND=i4), INTENT(IN) ::  &
  & icldscheme           ! cloud cover framework: see option above

REAL(KIND=wp), DIMENSION(klon,klev), INTENT(IN) ::   &
  & pp               , & ! full pressure                                 (  Pa )
  & pgeo             , & ! geopotential                                  (m2/s2)
  & rho                  ! density                                       (kg/m3)

REAL(KIND=wp), DIMENSION(klon,klev), INTENT(IN) ::   &
  & tt               , & ! temperature                                   (  K  )
  & qv               , & ! specific water vapor content                  (kg/kg)
  & qc               , & ! specific cloud water content                  (kg/kg)
  & qi                   ! specific cloud ice   content                  (kg/kg)

REAL(KIND=wp), DIMENSION(klon), INTENT(IN) ::   &
  & ps                   ! standard deviation of saturation deficit

REAL(KIND=wp), DIMENSION(klon,klev+1), INTENT(INOUT) ::   &
  & rcld                 ! standard deviation of saturation deficit

LOGICAL, DIMENSION(klon), INTENT(IN) ::  &
  & ldland           , & ! true for land points
  & ldcum                ! true for convection points

INTEGER(KIND=i4), DIMENSION(klon), INTENT(IN) ::  &
  & kcbot            , & ! convective cloud base level (klev: bottom level, -1: no conv)
  & kctop                ! convective cloud top level

REAL(KIND=wp), DIMENSION(klon,klev), INTENT(IN) ::   &
  & pmfude_rate      , & ! convective updraft detrainment rate           (kg/(m3*s))
  & plu                  ! updraft condensate                            (kg/kg)

REAL(KIND=wp), DIMENSION(klon,klev), INTENT(OUT) ::   &
  & cc_tot           , & ! cloud cover diagnostic
  & qv_tot           , & ! specific water vapor content diagnostic       (kg/kg)
  & qc_tot           , & ! specific cloud water content diagnostic       (kg/kg)
  & qi_tot               ! specific cloud ice   content diagnostic       (kg/kg)


!! Local variables:
!! ----------------

LOGICAL ::  &
  & lprog_qi

INTEGER (KIND=i4) :: &
  & jl, jk,          &
  & itype_wcld, icldm_rad

REAL(KIND=wp), DIMENSION(klon,klev)  :: &
  & cc_turb, qc_turb, qi_turb, &
  & cc_conv, qc_conv, qi_conv, &
  & cc_grid, qc_grid, qi_grid, &
  & p0

REAL(KIND=wp) :: &
  & t_g(klon)

!! Local parameters:
!! -----------------

REAL(wp), PARAMETER  :: &
  & zcldlim = 1.0e-6_wp  ! threshold of cloud water/ice for cloud cover  (kg/kg)

REAL(KIND=wp) :: &
  & taudecay

!-----------------------------------------------------------------------

IF (kstart > 1) THEN
  ! Set cloud fields for stratospheric levels to zero
  DO jk = 1,kstart-1
    DO jl = kidia,kfdia
      qc_tot(jl,jk) = 0.0_wp
      qi_tot(jl,jk) = 0.0_wp
      cc_tot(jl,jk) = 0.0_wp
    ENDDO
  ENDDO
ENDIF

! Select desired cloud cover framework
SELECT CASE( icldscheme )

!-----------------------------------------------------------------------

! no clouds
CASE( 0 )

  DO jk = kstart,klev
    DO jl = kidia,kfdia
      qc_tot(jl,jk) = 0.0_wp
      qi_tot(jl,jk) = 0.0_wp
      cc_tot(jl,jk) = 0.0_wp
    ENDDO
  ENDDO

!-----------------------------------------------------------------------

! grid-scale cloud cover [1 or 0]
CASE( 1 )

  DO jk = kstart,klev
    DO jl = kidia,kfdia
      qc_tot(jl,jk) = qc(jl,jk)
      qi_tot(jl,jk) = qi(jl,jk)
      IF ( qc(jl,jk) + qi(jl,jk) > zcldlim ) THEN
        cc_tot(jl,jk) = 1.0_wp
      ELSE
        cc_tot(jl,jk) = 0.0_wp
      ENDIF
    ENDDO
  ENDDO

!-----------------------------------------------------------------------

! cloud cover as in turbulence
CASE( 2 )

  itype_wcld = 2     ! 2: Gaussian calculation; 1: Sundqvist type
  p0         = 0.0_wp! base state presssure set to zero, pp is full pressure

  CALL cloud_diag ( cc_tot, qc_tot,                   &
                    kidia, kfdia, 1, 1, kstart, klev, &
                    klon , 1 , klev,                  &
                    tt, qv, qc, pp, p0, rcld, ps,     &
                    itype_wcld )

  qi_tot     = 0.0_wp

!-----------------------------------------------------------------------

! clouds as in COSMO
CASE( 3 )

  lprog_qi   = .true.       ! .true.: running with cloud ice
  icldm_rad  = 3            ! 3:      standard COSMO/GME
  itype_wcld = 2            ! cloud calculation in cloud_diag (only active if used)
  p0         = 0.0_wp       ! base state presssure set to zero, pp is full pressure
  DO jl = kidia,kfdia
    t_g(jl)  = tt(jl,klev)  ! should be surface temperature!
  ENDDO

  CALL cover_cosmo ( &
    klon       ,    & ! number of grid points in zonal direction
    1          ,    & ! number of grid points in meridional direction
    klev       ,    & ! number of grid points in vertical direction
    kidia      ,    & ! start- and end-indices for computing the radiation (i-start)
    1          ,    & !    -"- (j-start)
    kfdia      ,    & !    -"- (i-end)
    1          ,    & !    -"- (j-end)
    kstart     ,    & ! vertical start index
    tt         ,    & ! temperature                                   (  k  )
    ps         ,    & ! surface pressure                              ( pa  )
    p0         ,    & ! reference pressure at full levels             ( pa  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )
    pgeo       ,    & ! geopotential                                  (m2/s2)
    qv         ,    & ! specific water vapor content                  (kg/kg)
    qc         ,    & ! specific cloud water content                  (kg/kg)
    qi         ,    & ! specific cloud ice content                    (kg/kg)
    rcld       ,    & ! standard deviation of saturation deficit
    t_g        ,    & ! weighted surface temperature                  (  k  )
    cc_turb    ,    & ! OUT: subgrid-scale stratiform cloud cover     (  1  )
    cc_conv    ,    & ! OUT: cloud cover due to convection            (  1  )
    icldm_rad  ,    & ! mode of cloud representation in radiation parametr.
    itype_wcld ,    & ! type of water cloud diagnosis
    lprog_qi   ,    & ! if .TRUE., running with cloud ice
    ldland     ,    & ! landpoint mask
    ldcum      ,    & ! in:  convection on/off
    kcbot      ,    & ! in:  convective cloud base
    kctop      ,    & ! in:  convective cloud top
    qv_tot     ,    & ! OUT: water vapour mixing ratio
    cc_tot     ,    & ! OUT: cloud cover in each layer
    qc_tot     ,    & ! OUT: liquid water mixing ratio
    qi_tot     )      ! OUT: ice mixing ratio

 !cc_tot = cc_turb   ! CC=turb for testing
 !cc_tot = cc_conv   ! CC=conv for testing

!-----------------------------------------------------------------------

! diagnostic cloud cover
CASE( 4 )

  DO jk = kstart,klev
    DO jl = kidia,kfdia

     !turbulence
      cc_turb(jl,jk) = 0.0_wp
      qc_turb(jl,jk) = 0.0_wp
      qi_turb(jl,jk) = 0.0_wp

     !convection
      taudecay = 1800_wp
      cc_conv(jl,jk) = pmfude_rate(jl,jk) / rho(jl,jk) * taudecay             ! cc = detrainment / rho * tau,decay
      qc_conv(jl,jk) = cc_conv(jl,jk) * plu(jl,jk)*      foealfa(tt(jl,jk))   ! ql up  foealfa = liquid/(liquid+ice)
      qi_conv(jl,jk) = cc_conv(jl,jk) * plu(jl,jk)*( 1 - foealfa(tt(jl,jk)) ) ! qi up
      cc_conv(jl,jk) = min(max(0.0_wp,cc_conv(jl,jk)),1.0_wp)
      qc_conv(jl,jk) = min(max(0.0_wp,qc_conv(jl,jk)),0.1_wp*qv(jl,jk))       ! qc limit to 10%qv
      qi_conv(jl,jk) = min(max(0.0_wp,qi_conv(jl,jk)),0.1_wp*qv(jl,jk))       ! qi limit to 10%qv
     !cc_conv(jl,jk) = 0.0_wp
     !qc_conv(jl,jk) = 0.0_wp
     !qi_conv(jl,jk) = 0.0_wp

     !grid-scale
      IF ( qc(jl,jk) + qi(jl,jk) > zcldlim ) THEN
        cc_grid(jl,jk) = 1.0_wp
      ELSE
        cc_grid(jl,jk) = 0.0_wp
      ENDIF
      qc_grid(jl,jk) = qc(jl,jk)
      qi_grid(jl,jk) = qi(jl,jk)

     !combination
      cc_tot(jl,jk)  = max(cc_turb(jl,jk), cc_conv(jl,jk), cc_grid(jl,jk))
      qc_tot(jl,jk)  = max(qc_turb(jl,jk), qc_conv(jl,jk), qc_grid(jl,jk))
      qi_tot(jl,jk)  = max(qi_turb(jl,jk), qi_conv(jl,jk), qi_grid(jl,jk))

    ENDDO
  ENDDO

!-----------------------------------------------------------------------

! prognostic total water variance
CASE( 5 )

  cc_tot = 0.0_wp
  qc_tot = 0.0_wp
  qi_tot = 0.0_wp

!-----------------------------------------------------------------------

END SELECT

! total water vapor by conservation of grid-scale total water

DO jk = 1,klev
  DO jl = kidia,kfdia
    qv_tot(jl,jk) = qv(jl,jk) + qc(jl,jk) + qi(jl,jk) - qc_tot(jl,jk) - qi_tot(jl,jk)
    qv_tot(jl,jk) = MAX(qv_tot(jl,jk), 0.0_wp)
  ENDDO
ENDDO

END SUBROUTINE cover_koe


END MODULE mo_cover_koe
