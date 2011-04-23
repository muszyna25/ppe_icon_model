MODULE mo_cover_cosmo

USE mo_kind, ONLY: &
    ireals   =>wp       , &
    iintegers=>i4

USE mo_physical_constants, ONLY: &
    o_m_rdv             , & !! 1 - r_d/r_v
    rdv                 , & !! r_d / r_v
    cpdr     => rcpd    , & !! 1 / cp_d
    rdocp    => rd_o_cpd, & !! r_d / cp_d
    lh_v     => alv     , & !! latent heat of vapourization
    b3       => tmelt   , & !! melting temperature of ice/snow
    t0_melt  => tmelt   , & !! melting temperature of ice/snow
    grav                , & !! acceleration due to gravity
    uc1                 , & !! variables for computing the rate of cloud cover in
    ucl                     !! the unsaturated case

USE mo_math_constants, ONLY: &
    uc2 => sqrt3

USE mo_convect_tables, ONLY: &
    b1       => c1es    , & !! constants for computing the sat. vapour
    b2w      => c3les   , & !! pressure over water (l) and ice (i)
    b2i      => c3ies   , & !!               -- " --
    b4w      => c4les   , & !!               -- " --
    b4i      => c4ies       !!               -- " --

USE mo_cloud_diag,         ONLY: cloud_diag


!--------------------------------------------------------------------------------
! module variables from subroutine organized_radiation now used as arguments:

!USE data_modelconfig, ONLY :   &
!    ie,           & ! number of grid points in zonal direction
!    je,           & ! number of grid points in meridional direction
!    ke,           & ! number of grid points in vertical direction

!USE data_radiation, ONLY : &
!    jstartrad,       & !    fesft are computed on all grid points, to compute an
!    jendparrad,      & ! end-index just for preparations
!    istartrad,       & ! start- and end-indices for computing the radiation
!    iendrad,         & !   (when running on a coarser grid, the input values for

!USE data_fields     , ONLY :   &
!    t          ,    & ! temperature                                   (  k  )
!    ps         ,    & ! surface pressure                              ( pa  )
!    p0         ,    & ! reference pressure at full levels             ( pa  )
!    pp         ,    & ! deviation from the reference pressure         ( pa  )
!    qc         ,    & ! specific cloud water content                  (kg/kg)
!    qi         ,    & ! specific cloud ice content                    (kg/kg)
!    t_g        ,    & ! weighted surface temperature                  (  k  )
!    clc_sgs    ,    & ! subgrid-scale stratiform cloud cover          ( --- )
!    clc_con    )      ! cloud cover due to convection                 ( --- )

!USE data_runcontrol , ONLY :   &
!    icldm_rad,    & ! mode of cloud representation in radiation  parametr.
!    lprog_qi,     & ! if .TRUE., running with cloud ice
!--------------------------------------------------------------------------------


IMPLICIT NONE

PRIVATE

PUBLIC :: cover_cosmo

!--------------------------------------------------------------------------------

CONTAINS

  !------------------------------------------------------------------------------
  !
  !  Cloud cover and cloud water/ice calculation from COSMO model.
  !  (extracted from subroutine organize_radiation
  !   in file cosmo/src_radiation.f90, section 3.1 & 3.2)
  !  Changes: 
  !  * 500hPa threshold for thin upper-level ice clouds has been taken out.
  !  * fix to SGS/grid-scale cloud water/ice interaction:
  !    if small grid-scale qc or qi then 
  !      qc,total=qc,sgs+qi,sgs
  !      qi,total=qc,sgs+qi,sgs
  !    This overestimates the SGS cloud water.
  !
  !  @par Revision History
  !  Initial version by Martin Koehler, DWD (2010-08-23)
  !
  !
  !  LITERATURE
  !  ....
  !
  !------------------------------------------------------------------------------

SUBROUTINE cover_cosmo ( &
    ie         ,    & ! number of grid points in zonal direction
    je         ,    & ! number of grid points in meridional direction
    ke         ,    & ! number of grid points in vertical direction
    istartrad  ,    & ! start- and end-indices for computing the radiation (i-start)
    jstartrad  ,    & !    -"- (j-start)
    iendparrad ,    & !    -"- (i-end)
    jendparrad ,    & !    -"- (j-end)
    kstart     ,    & ! vertical start index of region where moist physics is active
    t          ,    & ! temperature                                   (  k  )
    ps         ,    & ! surface pressure                              ( pa  )
    p0         ,    & ! reference pressure at full levels             ( pa  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )
    pgeo       ,    & ! geopotential                                  (m2/s2 )
    qv         ,    & ! specific water vapor content                  (kg/kg)
    qc         ,    & ! specific cloud water content                  (kg/kg)
    qi         ,    & ! specific cloud ice content                    (kg/kg)
    rcld       ,    & ! standard deviation of saturation deficit
    t_g        ,    & ! weighted surface temperature                  (  k  )
    clc_sgs    ,    & ! subgrid-scale stratiform cloud cover          ( --- )
    clc_con    ,    & ! cloud cover due to convection                 ( --- )
    icldm_rad  ,    & ! mode of cloud representation in radiation parametr.
    itype_wcld ,    & ! type of water cloud diagnosis
    lprog_qi   ,    & ! if .TRUE., running with cloud ice
    llandmask  ,    & ! landpoint mask
    ldcum      ,    & ! in:  convection on/off
    kcbot      ,    & ! in:  convective cloud base
    kctop      ,    & ! in:  convective cloud top
    zwv        ,    & ! OUT: water vapour mixing ratio
    zclc       ,    & ! OUT: cloud cover in each layer
    zclwc      ,    & ! OUT: liquid water mixing ratio
    zciwc      )      ! OUT: ice mixing ratio

!--------------------------------------------------------------------------------

!! Subroutine arguments:
!! --------------------

INTEGER(KIND=iintegers), INTENT(IN) ::  &
  & ie               , & ! number of grid points in zonal direction
  & je               , & ! number of grid points in meridional direction
  & ke               , & ! number of grid points in vertical direction
  & istartrad        , & ! start- and end-indices for computing the radiation
  & jstartrad        , & !  -"-
  & iendparrad       , & ! end-index just for preparations
  & jendparrad       , & !  -"-
  & kstart

INTEGER(KIND=iintegers), PARAMETER ::  &
  & nnow=1           , & ! corresponds to ntstep
  & nnew=1           , & ! corresponds to ntstep + 1
  & nzx=1                ! time level of prognostic variables

REAL(KIND=ireals), DIMENSION(ie,je,ke,nnew), INTENT(IN) ::   &
  & t                , & ! temperature                                   (  k  )
  & pp               , & ! deviation from the reference pressure         ( pa  )
  & qv               , & ! specific water vapor content                  (kg/kg)
  & qc               , & ! specific cloud water content                  (kg/kg)
  & qi                   ! specific cloud ice content                    (kg/kg)

REAL(KIND=ireals), DIMENSION(ie,je,ke+1), INTENT(IN) ::   &
  & rcld                 ! standard deviation of saturation deficit

REAL(KIND=ireals), DIMENSION(ie,je,nnew), INTENT(IN) ::   &
  & ps               , & ! surface pressure                              ( pa  )
  & t_g                  ! weighted surface temperature                  (  k  )

REAL(KIND=ireals), DIMENSION(ie,je,ke), INTENT(IN) ::   &
  & p0               , & ! reference pressure at full levels             ( pa  )
  & pgeo                 ! geopotential                                  (m2/s2 )

INTEGER(KIND=iintegers), INTENT(IN) ::  &
  & icldm_rad        , & ! mode of cloud representation in radiation  parametr.
  & itype_wcld           ! type of water cloud diagnosis

LOGICAL, INTENT(IN) :: &
  & lprog_qi         , & ! if .TRUE., running with cloud ice
  & llandmask(ie,je) , & ! landpoint mask
  & ldcum(ie,je)         ! true for convection points

INTEGER(KIND=iintegers), DIMENSION(ie,je), INTENT(IN) ::  &
  & kcbot            , & ! convective cloud base level (klev: bottom level, -1: no conv)
  & kctop                ! convective cloud top level

REAL(KIND=ireals), DIMENSION(ie,je,ke), INTENT(OUT) ::   &
  & clc_sgs          , & ! subgrid-scale stratiform cloud cover          ( --- )
  & clc_con          , & ! cloud cover due to convection                 ( --- )
  & zwv              , & ! OUT: Water vapour mixing ratio
  & zclc             , & ! OUT: Cloud cover in each layer
  & zclwc            , & ! OUT: liquid water mixing ratio
  & zciwc                ! OUT: ice mixing ratio

! Local (automatic) arrays:
! ------------------------

INTEGER(KIND=iintegers) ::    &
  i, j, k
! klv500,                     & ! k index of the LM-mainlevel, on 500 hPa

REAL(KIND=ireals), DIMENSION(ie,je,ke) :: &
  zsw                       , & ! Saturation water vapour mixing ratio over water
  zse                           ! Saturation water vapour mixing ratio over ice

REAL(KIND=ireals) ::          &
  fgew   , fgee   , fgqv   ,  & ! name of statement functions
  zt_ice1, zt_ice2,           &
  zph    , zsigma , zdthdz ,  &
  zpio   , zpiu   , zpim   ,  &
  zthvo  , zthvu  , zuc    ,  &
  zclwcs ,                    &
  zclics , zclws  ,           &
  zclick , zclwck , zclwk  ,  &
  zclwfk ,                    &
  zcs    , zck    ,           &
  zqdw   , zsex   , zf_ice ,  &
  ztt    , zzpv   , zzpa   ,  &
  zdtbot, zdttop

REAL(KIND=ireals), PARAMETER :: &
  zepclc = 1.0E-8_ireals,     & ! avoids cloud cover =1.0 and = 0.0
  zeph2o = 1.0E-9_ireals,     & ! minimum value for specific humidity
  zclwcm = 1.0E-9_ireals        ! avoids cloud water content = 0.0


!--------------------------------------------------------------------------------

! statement function to calculate saturation vapour pressure over water
  fgew(ztt)       = b1 * EXP( b2w*(ztt - b3)/(ztt - b4w) ) ! ztt: temperature

! statement function to calculate saturation vapour pressure over ice
  fgee(ztt)       = b1 * EXP( b2i*(ztt - b3)/(ztt - b4i) ) ! ztt: temperature

! statement function to calculate specific humitdity
  fgqv(zzpv,zzpa) = rdv*zzpv/(zzpa - o_m_rdv*zzpv)   ! zzpv: vapour pressure


! maximum (in-)cloud water content: 0.5% of specific humidity at saturation
  zclwfk=0.010_ireals !1.0% of specific humidity at saturation in convective clouds

!--------------------------------------------------------------------------------
! Section 3:  Set cloudiness and humidity on input for fesft;
!             Store cloud cover on corresponding global arrays
!--------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  ! Section 3.0: Calculate convective cloud cover by a simple empirical
  !              relation (following B.Ritter, FE14). Anvils are assumed for
  !              a temperature increase at top level
  !              (taken from cosmo/src_convection.f90)
  !------------------------------------------------------------------------------

  clc_con(:,:,:) = 0.0_ireals
  DO k = kstart, ke
    DO j = jstartrad, jendparrad       ! jstartpar, jendpar
      DO  i = istartrad, iendparrad    ! istartpar, iendpar
        IF( ldcum(i,j) .AND. kctop(i,j) > 0 .AND. k <= kcbot(i,j) .AND. k >= kctop(i,j) ) THEN
 
          clc_con(i,j,k) = 0.7_ireals/10000.0_ireals &
                       & * ( pgeo(i,j,kctop(i,j)) - pgeo(i,j,kcbot(i,j)) ) / grav

          IF ( k ==  kctop(i,j) ) THEN
            zdtbot = t(i,j,k+1,nzx) - t(i,j,k  ,nzx)
            zdttop = t(i,j,k  ,nzx) - t(i,j,k-1,nzx)
            IF ( zdtbot > 0.0_ireals .AND. zdttop <= 0.0_ireals ) THEN
              clc_con(i,j,k) = 2.0_ireals * clc_con(i,j,k)
            ENDIF
          ENDIF

          clc_con(i,j,k) = MIN ( 1.0_ireals, MAX(0.05_ireals, clc_con(i,j,k)) )

        ENDIF
      ENDDO
    ENDDO
  ENDDO


  !------------------------------------------------------------------------------
  ! Section 3.1: Calculate water vapour saturation mixing ratios of
  !              over water and over ice
  !------------------------------------------------------------------------------

  zt_ice1= t0_melt -  5.0_ireals
  zt_ice2= t0_melt - 25.0_ireals

  DO k = kstart, ke
    DO j = jstartrad, jendparrad       ! jstartpar, jendpar
      DO  i = istartrad, iendparrad    ! istartpar, iendpar

        ! specific humidity (zwv) specific total water content (zqdw),
        ! specific humidity at saturation
        ! over water (zsw ) and ice (zse)
        zph       = p0(i,j,k) + pp(i,j,k,nzx)
        zse (i,j,k) = fgqv ( fgee(t(i,j,k,nzx)), zph)
        zsw (i,j,k) = fgqv ( fgew(t(i,j,k,nzx)), zph)
      ENDDO
    ENDDO
  ENDDO

  !------------------------------------------------------------------------------
  ! Section 3.2: Calculate stratiform cloud cover (non-convective)
  !------------------------------------------------------------------------------

  IF     ( icldm_rad == 0 ) THEN

    ! a) No interpretation of clouds at all for radiative calculations
    !-----------------------------------------------------------------

    DO  k = kstart, ke
      DO j = jstartrad, jendparrad       ! jstartpar, jendpar
        DO  i = istartrad, iendparrad    ! istartpar, iendpar
          zclwc(i,j,k) = 0.0_ireals
          zciwc(i,j,k) = 0.0_ireals
          zclc (i,j,k) = 0.0_ireals
        ENDDO
      ENDDO
    ENDDO

  ELSEIF ( icldm_rad == 1 ) THEN

    ! b) Only grid-sale water clouds are passed to the radiation routine
    !-------------------------------------------------------------------

    DO  k = kstart, ke
      DO j = jstartrad, jendparrad       ! jstartpar, jendpar
        DO  i = istartrad, iendparrad    ! istartpar, iendpar
          zclwc(i,j,k) = qc(i,j,k,nzx)
          IF (lprog_qi) THEN
            IF ( qc(i,j,k,nzx)+qi(i,j,k,nnow) > 0.0_ireals ) THEN
              zclc(i,j,k) = 1.0_ireals
            ELSE
              zclc(i,j,k) = 0.0_ireals
            END IF
            zciwc(i,j,k) = qi(i,j,k,nnow)
          ELSE
            IF ( qc(i,j,k,nzx) > 0.0_ireals ) THEN
              zclc(i,j,k) = 1.0_ireals
            ELSE
              zclc(i,j,k) = 0.0_ireals
            END IF
            zciwc(i,j,k) = 0.0_ireals
          ENDIF
          clc_sgs(i,j,k) = zclc(i,j,k)
        ENDDO
      ENDDO
    ENDDO

  ELSEIF (icldm_rad == 2) THEN

    ! c) Cloud cover and water content from statistical diagnosis
    !------------------------------------------------------------

    CALL cloud_diag(zclc, zclwc,                                     &
         istartrad, iendparrad, jstartrad, jendparrad, kstart, ke,   &
         ie, je, ke,                                                 &
         t(:,:,:,nzx), qv(:,:,:,nzx), qc(:,:,:,nzx), pp(:,:,:,nzx),  &
         p0, rcld, ps(:,:,nzx),                                      &
         itype_wcld )

    DO  k = kstart, ke
      DO j = jstartrad, jendparrad       ! jstartpar, jendpar
        DO  i = istartrad, iendparrad    ! istartpar, iendpar
          ! convective (in-)cloud water content
          ! as a function of specific humidity at saturation
          IF ( t(i,j,k,nzx) >= t0_melt )  THEN
             zclwck = zsw(i,j,k)*zclwfk  !
          ELSE
             zclwck = zse(i,j,k)*zclwfk
          ENDIF

          ! cloud cover of the non convective part of the grid box and cloud ice
          zcs = zclc(i,j,k)
          zciwc(i,j,k) = 0.0_ireals

          IF (lprog_qi) THEN
            ! if there is a grid scale cloud with cloud ice,
            ! even there might has been diagnosed subgrid scale water clouds,
            ! their water is thought to be distributed over the
            ! whole grid volume:
            IF ( qi(i,j,k,nnow) > 0.0_ireals ) THEN
              zcs = 1.0_ireals
            ENDIF
            zciwc(i,j,k) = qi(i,j,k,nnow)
          ENDIF
          clc_sgs(i,j,k) = zcs

          ! convective cloud cover
          zck = clc_con(i,j,k)

          ! grid scale cloud cover and water content
          zclc (i,j,k) = zcs + zck*(1.0_ireals-zcs)
          zclwc(i,j,k) = zclwc(i,j,k)*(1.0_ireals-zck) + zclwck*zck
        ENDDO
      ENDDO
    ENDDO

  ELSEIF ( icldm_rad == 4 .OR. icldm_rad == 3 ) THEN

    ! a) Standard diagnosis
    ! ---------------------

    DO  k = kstart, ke
      DO j = jstartrad, jendparrad       ! jstartpar, jendpar
        DO  i = istartrad, iendparrad    ! istartpar, iendpar
          ! Critical relative humidity as function of thermal stability
          zph      = p0(i,j,k) + pp(i,j,k,nzx)
          zsigma   = zph / ps(i,j,nzx)
          zdthdz   = 0.0_ireals

          zsex      = zsw(i,j,k)
          zqdw      = qv(i,j,k,nzx) + qc(i,j,k,nzx)
          IF (lprog_qi) THEN
            zf_ice      = 1.0_ireals - MIN( 1.0_ireals, MAX( 0.0_ireals, &
                          (t(i,j,k,nzx)-zt_ice2)/(zt_ice1-zt_ice2) ) )
            zqdw        = zqdw      + qi(i,j,k,nnow)
            zsex        = zsw(i,j,k) * (1.0_ireals - zf_ice) + zse(i,j,k)*zf_ice
          ENDIF

          IF(k == ke) THEN
            zpio    = ( 1.0E-5_ireals *( p0(i,j,k)+pp(i,j,k,nzx) ) )**rdocp
            zpiu    = ( 1.0E-5_ireals * ps(i,j,nzx)  )**rdocp
            zpim    = 0.5_ireals*(zpio+zpiu)
            zthvo   = t  (i,j,k  ,nzx)/zpio
            zthvu   = t_g(i,j,    nzx)/zpiu
            zdthdz  = zthvo - zthvu
          ELSE IF(zsigma.GT.0.95_ireals) THEN
            zpio    = ( 1.0E-5_ireals *( p0(i,j,k  )+pp(i,j,k  ,nzx) ) )**rdocp
            zpiu    = ( 1.0E-5_ireals *( p0(i,j,k+1)+pp(i,j,k+1,nzx) ) )**rdocp
            zpim    = 0.5_ireals*(zpio+zpiu)
            zthvo   = t(i,j,k  ,nzx)/zpio
            zthvu   = t(i,j,k+1,nzx)/zpiu
            zdthdz  = zthvo - zthvu + (lh_v*cpdr/zpim)*(qv(i,j,k,nzx)-qv(i,j,k+1,nzx))
          ENDIF

          ! grid scale cloud cover as function of relative humidity
          zuc     = 0.95_ireals - uc1*zsigma*(1._ireals-zsigma)*(1._ireals+uc2*(zsigma-0.5_ireals))
          zcs     = MAX ( 0.0_ireals,                &
                    MIN ( 1.0_ireals, (zqdw/zsex-zuc)/(ucl-zuc) ) )**2

          ! Corrections and limitations
          IF ( (zsigma > 0.95_ireals) .AND. (zdthdz < 0.0_ireals) ) THEN
            zcs = 0.0_ireals  ! no cloud cover in unstable stratification
          ENDIF
          IF ( qc(i,j,k,nzx) > 0.0_ireals ) THEN  ! grid scale clouds
            IF ( llandmask(i,j) .AND. k < ke ) zcs = 1.0_ireals
          ENDIF
          IF (lprog_qi) THEN
            IF (qi(i,j,k,nnow) > 1.0E-7_ireals) THEN
              zcs = 1.0_ireals ! grid scale clouds with cloud ice
            ENDIF
          ENDIF

          ! store grid-scale cloud cover on global array
          clc_sgs(i,j,k) = zcs

          ! Maximum in-cloud water content:  1.0% of specific humidity at saturation
          !                                  except for convective clouds (fixed)
          ! Standard diagnosis

          IF (lprog_qi) THEN
            zclws  = 0.005_ireals*zsex
            zclwcs = zclws*(1.0_ireals-zf_ice)
            zclics = zclws*zf_ice

! OLD: Check for grid-scale water or ice-clouds
!            IF ( qc(i,j,k,nzx) > 0.0 .OR. qi(i,j,k,nnow) > 1.0E-7_ireals ) THEN
!              zclics = 0.0_ireals
!              zclwcs = 0.0_ireals
!              IF ( qc(i,j,k,nzx) > 0.0 ) THEN  ! grid scale cloud water
!                 zclwcs = MAX( zclws, 0.5_ireals*qc(i,j,k,nzx) )
!              ENDIF
!              IF ( qi(i,j,k,nnow) > 1.0E-7_ireals ) THEN  ! grid scale cloud ice
!                 zclics = MAX( zclws, 0.5_ireals*qi(i,j,k,nnow) )
!              ENDIF
!            ENDIF

            ! Check for grid-scale water or ice-clouds
            IF ( qc(i,j,k,nzx) > 0.0_ireals ) THEN      ! grid scale cloud water
              zclwcs = MAX( zclwcs, 0.5_ireals*qc(i,j,k,nzx) )
            ENDIF
            IF ( qi(i,j,k,nnow) > 1.0E-7_ireals ) THEN  ! grid scale cloud ice
              zclics = MAX( zclics, 0.5_ireals*qi(i,j,k,nnow) )
            ENDIF

            ! Convective cloud water / ice content
            zclwk  = MAX( 2.0_ireals*zclws, 0.0002_ireals )
            zclwck = zclwk*(1.0_ireals-zf_ice)
            zclick = zclwk*zf_ice

            ! Reduce the cloud cover of ice clouds in the upper troposphere
            ! for the diagnosis of clch and clct
!           IF ((k <= klv500) .AND. (zclwcs <= 1.0E-10_ireals) .AND. &
            IF (                    (zclwcs <= 1.0E-10_ireals) .AND. &
                                    (zclics  > 0.0_ireals) ) THEN
              clc_sgs(i,j,k) = clc_sgs(i,j,k)*MIN( 1._ireals, MAX(0.2_ireals, &
                                 ( LOG(zclics)       - LOG(1.E-7_ireals) )/ &
                                 ( LOG(5.E-5_ireals) - LOG(1.E-7_ireals) )) )
            ENDIF

            ! set area-average cloud water/ice content
            zclwc(i,j,k) = zclwck*clc_con(i,j,k) + &
                           zclwcs*clc_sgs(i,j,k)*(1.0_ireals-clc_con(i,j,k))
            zciwc(i,j,k) = zclick*clc_con(i,j,k) + &
                           zclics*clc_sgs(i,j,k)*(1.0_ireals-clc_con(i,j,k))

          ELSE

            zclwcs = 0.005_ireals*zsw(i,j,k)
            zclwck = MAX( zclwcs, 0.0002_ireals )
            IF ( qc(i,j,k,nzx) > 0.0_ireals ) THEN  ! grid scale clouds
               zclwcs = MAX( zclwcs, 0.5_ireals*qc(i,j,k,nzx) )
            ENDIF

            ! set area-average cloud water/ice content
            zclwc(i,j,k) = zclwck*clc_con(i,j,k) + &
                           zclwcs*clc_sgs(i,j,k)*(1.0_ireals-clc_con(i,j,k))

            ! set area average cloud ice content (in-cloud)
            zciwc(i,j,k) = 0.0_ireals

          ENDIF

          ! calculate combined cloud cover
          zclc (i,j,k) = clc_sgs(i,j,k) + &
                         clc_con(i,j,k)*( 1.0_ireals - clc_sgs(i,j,k) )

        ENDDO
      ENDDO
    ENDDO

  ENDIF  ! icldm_rad

  ! Restrictions for radiative calculations
  ! ---------------------------------------

  DO  k = kstart, ke
    DO j = jstartrad, jendparrad       ! jstartpar, jendpar
      DO  i = istartrad, iendparrad    ! istartpar, iendpar
        zwv  (i,j,k) = MIN( MAX(zeph2o,qv(i,j,k,nzx)), zsw(i,j,k) )
        zclc (i,j,k) = MAX( zepclc, MIN(1.0_ireals-zepclc,zclc(i,j,k)) )
        zclwc(i,j,k) = MAX( zclwcm, zclwc(i,j,k) )
        zciwc(i,j,k) = MAX( zclwcm, zciwc(i,j,k) )
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE cover_cosmo

!================================================================================

END MODULE mo_cover_cosmo
