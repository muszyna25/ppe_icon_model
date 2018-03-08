!>
!! Configuration of the parameterization for atmospheric gravity wave effects,
!! that is used in the ECHAM physics package.
!!
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!     First version by Marco Giorgetta, MPI-M (2017-11)
!!
!! Based on earlier codes of:
!!     Elisa Manzini, MPI-M, 1997
!!
!! References:
!!     Hines, C. O. (1997a). Doppler spread parameterization of gravity wave momentum deposition
!!                           in the middle atmosphere. Part I: Basic formulation.
!!                           J. Atmos. Solar Terr. Phys., 59, 371– 386.
!!     Hines, C. O. (1997b). Doppler spread parameterization of gravity wave momentum deposition
!!                           in the middle atmosphere. Part II: Broad and quasi monochromatic
!!                           spectra and implemen- tation. J. Atmos. Solar Terr. Phys., 59, 387–400.
!!     Manzini, E., McFarlane, N. A. and McLandress, C. (1997). Impact of the Doppler spread
!!                           parameterization on the simulation of the middle atmosphere circulation using
!!                           the MA/ECHAM4 general circulation model. J. Geophys. Res., 102, 25751–25762.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_gwd_config

  USE mo_exception            ,ONLY: message, print_value, finish
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom
  USE mo_grid_config          ,ONLY: n_dom
  USE mo_run_config           ,ONLY: nlev

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         echam_gwd_config   !< user specified configuration parameters
  PUBLIC ::    init_echam_gwd_config   !< allocate and initialize echam_gwd_config
  PUBLIC ::    eval_echam_gwd_config   !< evaluate echam_gwd_config
  PUBLIC ::   print_echam_gwd_config   !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'echam_gwd'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the ECHAM atmospheric gravity wave drag
  !!
  TYPE t_echam_gwd_config
     !
     ! configuration parameters
     ! ------------------------
     !
     LOGICAL  :: lheatcal      !< true : compute momentum flux deposition, heating and diffusion coefficient
     !                            false: compute momentum flux deposition only
     !
     INTEGER  :: emiss_lev     !< number of levels above the ground at which gw are emitted
     REAL(wp) :: rmscon        !< [m/s] root mean square gravity wave wind at emission level
     REAL(wp) :: kstar         !< [1/m] typical gravity wave horizontal wavenumber
     REAL(wp) :: m_min         !< [1/m] minimum bound in vertical wavenumber
     !
!!$     LOGICAL  :: lfront        !< true: compute gw sources emerging from fronts and background
!!$                               !< (Charron and Manzini, 2002)
!!$     REAL(wp) :: rms_front     !< [m/s] rms frontal gw wind at source level
!!$     REAL(wp) :: front_thres   !< [(K/m)^2/hour] minimum value of the frontogenesis function,
!!$     !                            for whichgravity waves are emitted from fronts
!!$     !
!!$     LOGICAL  :: lozpr         !< true: for background enhancement associated with precipitation
!!$                               !< (Manzini et al., 1997)
!!$     REAL(wp) :: pcrit         !< [mm/d] critical precipitation value, above which 
!!$                               !< gravity wave rms wind enhancement is applied
!!$     REAL(wp) :: pcons         !< [] adimensional factor for background enhancement 
!!$     !                            associated with precipitation
!!$     !
!!$     LOGICAL  :: lrmscon_lat   !< true:  use latitude dependent rmscon
!!$     !                            - |latitude| >= lat_rmscon:
!!$     !                                 use rmscon
!!$     !                            - |latitude| <= lat_rmscon_eq:
!!$     !                                 use rmscon_eq
!!$     !                            - lat_rmscon_eq < |latitude| < lat_rmscon: 
!!$     !                                 use linear interpolation between rmscon_eq and rmscon
!!$     !                            false: use rmscon for all latitudes
!!$     !                            attention: may be overwritten if lfront or lozpr is true
!!$     REAL(wp) :: lat_rmscon_eq !< [degN] rmscon_eq is used equatorward of this latitude
!!$     REAL(wp) :: lat_rmscon    !< [degN] rmscon is used poleward of this latitude
!!$     REAL(wp) :: rmscon_eq     !< [m/s]  rms constant used equatorward of lat_rmscon_eq
     !
  END TYPE t_echam_gwd_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_echam_gwd_config), TARGET :: echam_gwd_config(max_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_echam_gwd_config
    !
    ! ECHAM atmospheric gravity wave drag configuration
    ! -------------------------------------------------
    !
    echam_gwd_config(:)% lheatcal      = .FALSE.    ! normally do not use heating
    !
    ! For global uniform gravity wave sources:
    echam_gwd_config(:)% emiss_lev     = 10         ! is correct for L31 and L47
    echam_gwd_config(:)% rmscon        = 0.87_wp    ! default value used in ECHAM5
    echam_gwd_config(:)% kstar         = 5.0e-5_wp  ! = 2*pi/(126000 m)
    echam_gwd_config(:)% m_min         = 0.0_wp
    !
!!$    ! For latitude dependent rms wind:
!!$    echam_gwd_config(:)% lrmscon_lat   = .FALSE.    ! ICON default: no latitude dependence
!!$    echam_gwd_config(:)% lat_rmscon_eq =  5.0_wp    ! as in ECHAM6 T63
!!$    echam_gwd_config(:)% lat_rmscon    = 10.0_wp    ! as in ECHAM6 T63
!!$    echam_gwd_config(:)% rmscon_eq     = 1.2_wp     ! as in ECHAM6 T63
    !
  END SUBROUTINE init_echam_gwd_config

  !----

  !>
  !! Evaluate additional derived parameters
  !!
  SUBROUTINE eval_echam_gwd_config
    !
    CHARACTER(LEN=*), PARAMETER :: routine = 'eval_echam_gwd_config'
    !
    !------------------------------------------------------------------
    ! Sanity Check
    !------------------------------------------------------------------
    IF ( ANY(echam_gwd_config(:)%emiss_lev < 0     ) ) CALL finish(routine,'emiss_lev < 0 is not allowed')
    IF ( ANY(echam_gwd_config(:)%emiss_lev > nlev  ) ) CALL finish(routine,'emiss_lev > nlev is not allowed')
    IF ( ANY(echam_gwd_config(:)%rmscon    < 0.0_wp) ) CALL finish(routine,'rmscon < 0. is not allowed')
    IF ( ANY(echam_gwd_config(:)%kstar     < 0.0_wp) ) CALL finish(routine,'kstar  < 0. is not allowed')
    IF ( ANY(echam_gwd_config(:)%m_min     < 0.0_wp) ) CALL finish(routine,'m_min  < 0. is not allowed')
    !
  END SUBROUTINE eval_echam_gwd_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_echam_gwd_config
    !
    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','ECHAM atmospheric gravity wave drag configuration')
    CALL message    ('','=================================================')
    CALL message    ('','')
    CALL message    ('','')
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    echam_gwd_config('//TRIM(cg)//')% emiss_lev  ',echam_gwd_config(jg)% emiss_lev )
       CALL print_value('    echam_gwd_config('//TRIM(cg)//')% rmscon     ',echam_gwd_config(jg)% rmscon    )
       CALL print_value('    echam_gwd_config('//TRIM(cg)//')% kstar      ',echam_gwd_config(jg)% kstar     )
       CALL print_value('    echam_gwd_config('//TRIM(cg)//')% m_min      ',echam_gwd_config(jg)% m_min     )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_echam_gwd_config

  !----

END MODULE mo_echam_gwd_config
