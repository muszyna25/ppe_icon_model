!>
!! Configuration of the parameterization for sub-grid scale orographic effects,
!! that is used in the ECHAM physics package.
!!
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!     First version by Marco Giorgetta, MPI-M (2017-04)
!!
!! Based on earlier codes of:
!!     Martin Miller, ECMWF, Jan 1990
!!     Francois Lott, LMD,   Jul 1999  
!!     Elisa Manzini, MPI-M, Aug 2000
!!
!! References: 
!!     Lott, 1999: Alleviation of stationary biases in a GCM through...
!!                 Monthly Weather Review, 127, pp 788-801.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_sso_config

  USE mo_exception            ,ONLY: message, print_value
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom
  USE mo_grid_config          ,ONLY: n_dom
  USE mo_run_config           ,ONLY: nlev
  USE mo_vertical_coord_table ,ONLY: vct

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         echam_sso_config   !< user specified configuration parameters
  PUBLIC ::    init_echam_sso_config   !< allocate and initialize echam_sso_config
  PUBLIC ::    eval_echam_sso_config   !< evaluate echam_sso_config
  PUBLIC ::   print_echam_sso_config   !< print out

  ! parameters
  PUBLIC :: gfrcrit, grcrit, grahilo
  PUBLIC :: gsigcr , gssec , gtsec , gvsec

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'echam_sso'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the ECHAM subgrid scale orographic drag
  !!
  TYPE t_echam_sso_config
     !
     ! configuration parameters
     ! ------------------------
     !
     ! threshold values for defining the mask of active points
     REAL(wp) :: gpicmea   ! minimum difference "SSO peak height - SSO mean height" [m]
     REAL(wp) :: gstd      ! minimum standard deviation of SSO height [m]
     !
     ! parameters controling the strength of the effects
     REAL(wp) :: gkdrag    ! Gravity wave drag coefficient                  (G  in (3), LOTT 1999)
     REAL(wp) :: gkwake    ! Bluff-body drag coefficient for low level wake (Cd in (2), LOTT 1999)
     REAL(wp) :: gklift    ! Mountain Lift coefficient                      (Cl in (4), LOTT 1999)
     !
     ! parameters related to the vertical grid
     INTEGER  :: nktopg    ! Security value for blocked flow level
     INTEGER  :: ntop      ! An estimate to qualify the upper levels of the
     !                       model where one wants to impose stress profiles
     !
     ! scaling with sftlf, the cell area fraction of land incl. lakes
     LOGICAL  :: lsftlf    ! true: *lsftlf, false: *1

  END TYPE t_echam_sso_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_echam_sso_config), TARGET :: echam_sso_config(max_dom)
  
  ! "tunable parameters" of the various SSO schemes, same for all domains
  !
  REAL(wp), PARAMETER :: gfrcrit = 0.50_wp ! Critical Non-dimensional mountain Height (HNC in (1), LOTT 1999)
  REAL(wp), PARAMETER :: grcrit  = 0.25_wp ! Critical Richardson Number (Ric, end of first column p791, LOTT 1999)
  REAL(wp), PARAMETER :: grahilo = 1.00_wp ! Set-up the trapped waves fraction (Beta , end of first column, LOTT 1999)

  ! numerical security parameters, same for all domains
  !
  REAL(wp), PARAMETER :: gsigcr = 0.80_wp      ! Security value for blocked flow depth
  REAL(wp), PARAMETER :: gssec  = 0.0001_wp    ! Security min value for low-level B-V frequency
  REAL(wp), PARAMETER :: gtsec  = 0.00001_wp   ! Security min value for anisotropy and GW stress.
  REAL(wp), PARAMETER :: gvsec  = 0.10_wp      ! Security min value for ulow

CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_echam_sso_config
    !
    ! ECHAM subgrid scale orographic drag configuration
    ! -------------------------------------------------
    !
    ! Define the mask for the SSO parameterization:
    echam_sso_config(:)% gpicmea = 40.0_wp ! only where  (peak - mean height) is typically > 1st layer depth
    echam_sso_config(:)% gstd    = 10.0_wp ! only where SSO slope, asymmetry and orientation are defined by EXTPAR
    !
    ! Define the tuning parameters for SSO drag. These values depend on:
    ! (1) the resolution of the topography data used to compute the SSO parameters, and
    ! (2) the model resolution.
    ! A 0-value switches the relevant effect off.
    echam_sso_config(:)% gkdrag  = 0.10_wp
    echam_sso_config(:)% gkwake  = 0.01_wp
    echam_sso_config(:)% gklift  = 0.00_wp
    !
    ! parameters related to the vertical grid
    echam_sso_config(:)% ntop    = 1
    echam_sso_config(:)% nktopg  = 0 ! needs to be derived
    !
    ! cell area fraction scaling for SSO effects, .FALSE.: *1, .TRUE.: *sftlf
    echam_sso_config(:)% lsftlf  = .TRUE.
    !
  END SUBROUTINE init_echam_sso_config

  !----

  !>
  !! Evaluate additional derived parameters
  !!
  SUBROUTINE eval_echam_sso_config
    !
    INTEGER          :: jk
    REAL(wp)         :: zsigt, zpm1r, zpr
    !
    ! height sigma grid
    !
    zpr   =  1950._wp ! m (800 hPa in International Standard Atmosphere)
    zsigt =  2460._wp ! m (750 hPa in International Standard Atmosphere)
    !
    DO jk=nlev,1,-1
       !
       ! full level height zf(jk) = (zh(jk)+zh(jk+1)/2 for z_sfc=zpr
       zpm1r = 0.5_wp*(vct(jk)+vct(jk+1)+zpr*(vct(nlev+1+jk)+vct(nlev+1+jk+1)))
       !
       ! Find highest full level with zf(jk) <= zsigt 
       IF (zpm1r <= zsigt) THEN
          echam_sso_config(:)% nktopg = jk
       END IF
       !
    END DO
    !
  END SUBROUTINE eval_echam_sso_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_echam_sso_config
    !
    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','ECHAM subgrid scale orographic drag configuration')
    CALL message    ('','=================================================')
    CALL message    ('','')
    !
    CALL print_value('    Critical Froude     number gfrcrit',gfrcrit )
    CALL print_value('    Critical Richardson number grcrit ',grcrit  )
    CALL message    ('','')
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL message    ('','Use SSO where peak-mean height > gpicmea and orostdh>gstd:')
       CALL print_value('    echam_sso_config('//TRIM(cg)//')% gpicmea  ',echam_sso_config(jg)% gpicmea  )
       CALL print_value('    echam_sso_config('//TRIM(cg)//')% gstd     ',echam_sso_config(jg)% gstd     )
       CALL message    ('','Coefficients for gravity wave drag, low level blocking and lift:')
       CALL print_value('    echam_sso_config('//TRIM(cg)//')% gkdrag   ',echam_sso_config(jg)% gkdrag   )
       CALL print_value('    echam_sso_config('//TRIM(cg)//')% gkwake   ',echam_sso_config(jg)% gkwake   )
       CALL print_value('    echam_sso_config('//TRIM(cg)//')% gklift   ',echam_sso_config(jg)% gklift   )
       CALL message    ('','')
       !
       CALL print_value('    echam_sso_config('//TRIM(cg)//')% lsftlf   ',echam_sso_config(jg)% lsftlf   )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_echam_sso_config

  !----

END MODULE mo_echam_sso_config
