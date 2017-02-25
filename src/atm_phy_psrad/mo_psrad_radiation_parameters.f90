!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!++mgs: new module 14.03.2010
!++mgs: added decl_sun_cur (for MOZ photolysis) 02.06.2010
!! @brief Module to provide parameters to radiation routines and avoid circular dependencies.
!!
!! @remarks
!!   This module contains the public parameters provided by the radiation module
!!   mo_radiation.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2009-09-19):
!!
!!         Martin Schultz, FZJ, Juelich (2010-04-13):
!!              extracted parameters from mo_radiation
!!         Dagmar Popke, MPI-M, Hamburg (2013-11-15):
!!              Implementation of RCE
!!
!! $ID: n/a$
!!
!! @par Origin
!!   see mo_radiation.f90
!!
!
MODULE mo_psrad_radiation_parameters

  USE mo_kind,            ONLY: wp
  USE mo_model_domain,    ONLY: t_patch
  USE mo_parallel_config, ONLY: nproma
  USE mo_math_constants,  ONLY: pi, pi2, pi_2, pi_4 ! pi, pi*2, pi/2, pi/4

IMPLICIT NONE

  PRIVATE

  PUBLIC :: irad_aero_forcing
  PUBLIC :: lw_spec_samp, sw_spec_samp, lw_gpts_ts, sw_gpts_ts, rad_perm
  PUBLIC :: i_overlap, l_do_sep_clear_sky
  PUBLIC :: cemiss, rad_undef
  PUBLIC :: psctm, ssi_factor
  
  ! 1.0 NAMELIST global variables and parameters
  ! --------------------------------
  INTEGER :: irad_aero_forcing = 0 !< reference aerosols for instantaneous 
                                   !< radiative forcing
  
  ! Spectral sampling
  INTEGER :: lw_spec_samp = 1, sw_spec_samp = 1 ! 1 is broadband, 2 is MCSI, 3 and up are teams
  INTEGER :: lw_gpts_ts   = 1, sw_gpts_ts   = 1 ! Number of g-points per time step using MCSI
  INTEGER :: rad_perm = 0                       ! Integer for perturbing random number seeds

  ! Radiation driver
  LOGICAL :: l_do_sep_clear_sky = .TRUE. ! Compute clear-sky fluxes by removing clouds
  INTEGER :: i_overlap = 1               ! 1 = max-ran, 2 = max, 3 = ran
  !
  ! --- radiative transfer parameters
  !
  REAL(wp), PARAMETER :: cemiss = 0.996_wp  !< LW Emissivity Factor
  REAL(wp), PARAMETER :: rad_undef = -999._wp
  !
  !++hs
  REAL(wp) :: psctm                         !< orbit and time dependent solar constant for radiation time step
  REAL(wp) :: ssi_factor(14)                !< fraction of TSI in the 14 RRTM SW bands
  !--hs

  public solar_parameters

contains

  !---------------------------------------------------------------------------
  !>
  !! @brief Scans a block and fills with solar parameters
  !! 
  !! @remarks: This routine calculates the solar zenith angle for each
  !! point in a block of data.  For simulations with no diurnal cycle 
  !! the cosine of the zenith angle is set to its average value (assuming 
  !! negatives to be zero and for a day divided into nds intervals).  
  !! Additionally a field is set indicating the fraction of the day over 
  !! which the solar zenith angle is greater than zero.  Otherwise the field 
  !! is set to 1 or 0 depending on whether the zenith angle is greater or 
  !! less than 1. 
  !
  SUBROUTINE solar_parameters(decl_sun,    time_of_day,     &
       &                      icosmu0,     dt_ext,          &
       &                      ldiur,       l_sph_symm_irr,  &
       &                      p_patch,                      &
       &                      cos_mu0,     daylight_frc     )

    REAL(wp), INTENT(in)  :: &
         decl_sun,           & !< delination of the sun
         time_of_day,        & !< time_of_day (in radians)
         dt_ext                !< time interval overfor which the insolated area is extended
    INTEGER , INTENT(in)  :: &
         icosmu0               !< defines if and how cosmu0 is defined for extended sunlit areas
    LOGICAL               :: &
         ldiur,              & !< diurnal cycle ON (ldiur=.TRUE.) or OFF (ldiur=.FALSE.)
         l_sph_symm_irr        !< spherical symmetric irradiation ON (l_sph_symm_irr=.TRUE.)
                               !< or OFF (l_sph_symm_irr=.FALSE.)
    TYPE(t_patch), INTENT(in) ::      p_patch
    REAL(wp), INTENT(out) :: &
         cos_mu0(:,:),       & !< cos_mu_0, cosine of the solar zenith angle
         daylight_frc(:,:)     !< daylight fraction (0 or 1) with diurnal cycle

    INTEGER     :: i, j
    REAL(wp)    :: zen1, zen2, zen3, xx

    INTEGER, PARAMETER :: nds = 128 !< number of diurnal samples

    LOGICAL  , SAVE :: initialized_sincos = .FALSE.
    REAL (wp), SAVE :: cosrad(nds), sinrad(nds)
    REAL (wp)       :: xsmpl(nds), xsmpl_day(nds), xnmbr_day(nds)
    REAL (wp), ALLOCATABLE :: sinlon(:,:), sinlat(:,:), coslon(:,:), coslat(:,:)
    REAL (wp), ALLOCATABLE :: mu0(:,:)

    REAL (wp), PARAMETER   :: eps=1.e-9_wp ! to converge in ca. 30 iterations
    LOGICAL  , SAVE        :: initialized_mu0s = .FALSE.
    REAL (wp), SAVE        :: mu0s, sin_mu0s
    REAL (wp)              :: dmu0, dcos_mu0, mu0min, mu0max
    INTEGER  , SAVE        :: nmu0

    INTEGER :: nprom, npromz, nblks

    nprom=nproma
    npromz=p_patch%npromz_c
    nblks=p_patch%nblks_c

    IF (.NOT.ldiur.AND..NOT.initialized_sincos) THEN
       !
       ! sin and cos arrays for zonal mean radiation
       DO i = 1, nds
          xx = pi2*(i-1.0_wp)/nds
          sinrad(i) = SIN(xx)
          cosrad(i) = COS(xx)
       END DO
       !
       initialized_sincos = .TRUE.
    END IF
    !
    IF (ldiur.AND.(icosmu0==4).AND.(.NOT.initialized_mu0s).AND.(dt_ext/=0.0_wp)) THEN
       !
       ! find mu0s for radiation with diurnal cycle
       dmu0   = dt_ext/86400._wp*pi
       mu0min = 0.0_wp
       mu0max = pi_2
       mu0s   = mu0max
       nmu0   = 0
       DO
          xx = (pi_2+dmu0-mu0s)*SIN(mu0s) - COS(mu0s)
          IF      (xx >  eps ) THEN
             mu0max = mu0s
          ELSE IF (xx < -eps ) THEN
             mu0min = mu0s
          ELSE
             EXIT  ! ABS(xx)<=eps, mu0s is found
          END IF
          mu0s = 0.5_wp*(mu0min+mu0max)
          nmu0 = nmu0+1
          !
       END DO
       !
       ! slope of the linear section of cos_mu0, see below, where used
       sin_mu0s = SIN(mu0s)
       !
       initialized_mu0s = .TRUE.
    END IF
    !
    zen1 = SIN(decl_sun)
    zen2 = COS(decl_sun)*COS(time_of_day)
    zen3 = COS(decl_sun)*SIN(time_of_day)

    IF (.NOT.l_sph_symm_irr) THEN       ! - spherically variable irradiation
       !
       ALLOCATE(sinlon(nprom,nblks))
       ALLOCATE(sinlat(nprom,nblks))
       ALLOCATE(coslon(nprom,nblks))
       ALLOCATE(coslat(nprom,nblks))
       ALLOCATE(mu0(nprom,nblks))
       sinlon(:,:)=0._wp
       sinlat(:,:)=0._wp
       coslon(:,:)=0._wp
       coslat(:,:)=0._wp
       sinlon(1:nprom,1:nblks-1)=SIN(p_patch%cells%center(1:nprom,1:nblks-1)%lon)
       sinlat(1:nprom,1:nblks-1)=SIN(p_patch%cells%center(1:nprom,1:nblks-1)%lat)
       coslon(1:nprom,1:nblks-1)=COS(p_patch%cells%center(1:nprom,1:nblks-1)%lon)
       coslat(1:nprom,1:nblks-1)=COS(p_patch%cells%center(1:nprom,1:nblks-1)%lat)
       sinlon(1:npromz,nblks)=SIN(p_patch%cells%center(1:npromz,nblks)%lon)
       sinlat(1:npromz,nblks)=SIN(p_patch%cells%center(1:npromz,nblks)%lat)
       coslon(1:npromz,nblks)=COS(p_patch%cells%center(1:npromz,nblks)%lon)
       coslat(1:npromz,nblks)=COS(p_patch%cells%center(1:npromz,nblks)%lat)
       !
       IF (ldiur) THEN                  ! - with local diurnal cycle
          !
          ! cos(zenith angle), positive for sunlit hemisphere
          cos_mu0(:,:)     =  zen1*sinlat(:,:)                &
               &             -zen2*coslat(:,:)*coslon(:,:)    &
               &             +zen3*coslat(:,:)*sinlon(:,:)
          !
          ! zenith angle
          mu0(:,:) = ACOS(cos_mu0(:,:))
          !
          ! increment of mu0 to include a rim of width dmu0= (dt_ext/2)*2pi/1day around
          ! the sunlit hemisphere at the radiation time so that the extended area includes
          ! the sunlit areas of all time step within the radiation interval of length dt_ext
          ! centered at the radiation time
          dmu0     = dt_ext/86400._wp*pi
          dcos_mu0 = SIN(dmu0)
          !
          ! add increment dcos_mu0 for the definition of the extended daylight area
          ! set day/night indicator to 1/0
          WHERE (cos_mu0(:,:)+dcos_mu0 < 0.0_wp)
             daylight_frc(:,:) = 0.0_wp
          ELSEWHERE
             daylight_frc(:,:) = 1.0_wp
          END WHERE
          !
          IF (dt_ext/=0.0_wp) THEN
             !
             
             SELECT CASE (icosmu0)
             CASE (0)
                !
                ! no modification -> nothing to do
                !
             CASE (1)
                !
                ! minimum value
                !
                WHERE (daylight_frc(:,:) == 1.0_wp)
                   cos_mu0(:,:) = MAX(0.1_wp,cos_mu0(:,:))
                END WHERE
                !
             CASE (2)
                !
                ! shift and rescale
                !
                WHERE (daylight_frc(:,:) == 1.0_wp)
                   cos_mu0(:,:) = (cos_mu0(:,:)+dcos_mu0)/(1._wp+dcos_mu0)
                END WHERE
                !
             CASE (3)
                !
                ! slope in band [pi/2-dmu0,pi/2+dmu0]
                !
                ! redefine cos_mu0 in [pi/2-dmu0,pi/2+dmu0] using a linear function of mu0 such that:
                !   inner edge                       : mu0=pi/2-dmu0 --> cos_mu0 = cos(pi/2-dmu0)
                !   center     = original terminator : mu0=pi/2      --> cos_mu0 = cos(pi/2-dmu0)/2
                !   outer edge = extended terminator : mu0=pi/2+dmu0 --> cos_mu0 = 0
                !
                WHERE (ABS(mu0(:,:)-pi_2)<dmu0)
                   cos_mu0(:,:) = 0.5_wp*SIN(dmu0)*(1._wp-(mu0(:,:)-pi_2)/dmu0)
                END WHERE
                !
             CASE (4)
                !
                ! tangent slope in [mu0s,pi/2+dmu0]
                !
                ! redefine cos_mu0 in [mu0s,pi/2+dmu0] using a linear function of mu0 such that:
                !   inner edge                       : mu0=mu0s      --> cos_mu0 = cos(mu0s)
                !   in between                       : mu0           --> cos_mu0 = sin(mu0s)*(pi/2+dmu0-mu0)
                !   outer edge = extended terminator : mu0=pi/2+dmu0 --> cos_mu0 = 0
                !
                ! mu0s is the solution of : cos(mu0s) = sin(mu0s)*(pi/2+dmu0-mu0s)
                ! so that cos_mu0 is a C1 function in [0,pi/2+dmu0]
                !
                WHERE ((mu0s<mu0(:,:)).AND.(mu0(:,:)<(pi_2+dmu0)))
                   cos_mu0(:,:) = sin_mu0s*(pi_2+dmu0-mu0(:,:))
                END WHERE
                !
             END SELECT
          END IF
          !
       ELSE                             ! - with zonally symmetric irradiation
          !
          ! For each cell (i,j), compute first cos(mu0) for nds regularly spaced longitudes on the
          ! latitude circle of the cell. Then compute the zonal mean of cos(mu0) for the latitude
          ! of this cell as the average over all longitudes, where cos(mu0)>=epsilon.
          !
          DO j = 1, SIZE(cos_mu0,2)
             DO i = 1, SIZE(cos_mu0,1)
                !
                ! cos(zenith angle) for nds longitudes on the latitude circle of cell (i,j)
                xsmpl(:) =  zen1*sinlat(i,j)              &
                     &     -zen2*coslat(i,j)*cosrad(:)    &
                     &     +zen3*coslat(i,j)*sinrad(:)
                !
                ! set day/night indicator of sampled longitudes to 1/0
                ! and mask out cosmu0 values in day-only and night-only cells
                !
                WHERE (xsmpl(:) < EPSILON(1.0_wp))     ! night side
                   xsmpl_day(:) = 0.0_wp
                   xnmbr_day(:) = 0.0_wp
                ELSEWHERE                              ! day side
                   xsmpl_day(:) = xsmpl(:)
                   xnmbr_day(:) = 1.0_wp
                END WHERE
                !
                cos_mu0(i,j)      = SUM(xsmpl_day(:))
                daylight_frc(i,j) = SUM(xnmbr_day(:))
                !
                IF (daylight_frc(i,j) > EPSILON(1.0_wp)) THEN         ! at least one point on the latitude
                   cos_mu0(i,j)      = cos_mu0(i,j)/daylight_frc(i,j) ! circle of the cell has daylight
                   daylight_frc(i,j) = daylight_frc(i,j)/REAL(nds,wp) ! and a zonal mean can be defined
                ELSE
                   cos_mu0(i,j)      = 0.0_wp                         ! set cos(mu0) = 0 on polar night
                   daylight_frc(i,j) = 0.0_wp                         ! latitude circle
                END IF
                !
             END DO
          END DO
          !
       END IF
       !
       DEALLOCATE(sinlon)
       DEALLOCATE(sinlat)
       DEALLOCATE(coslon)
       DEALLOCATE(coslat)
       DEALLOCATE(mu0)
       !
    ELSE                                ! - spherically symmetric irradiation (for RCE),
       !                                    all grid points have the same solar incoming
       !                                    flux at TOA
       !
       IF (ldiur) THEN                  ! - with diurnal cycle of (0degE,0degN), i.e.
          !                                 local noon is at 12:00 UTC in all points
          !
          cos_mu0(:,:) = -zen2          !   = cos_mu0 at (0degE,0degN)
          IF (-zen2 < 0.0_wp) THEN
             daylight_frc(:,:) = 0.0_wp
          ELSE
             daylight_frc(:,:) = 1.0_wp
          END IF
          !
       ELSE                             ! - without diurnal cycle
          !                                 all grid points have the same constant cos_mu0
          !
          cos_mu0(:,:) = pi_4           !  = pi/4 (why this choice?)
          daylight_frc(:,:) = 1.0_wp
          !
       END IF
       !
       !
    END IF

  END SUBROUTINE solar_parameters

END MODULE mo_psrad_radiation_parameters
