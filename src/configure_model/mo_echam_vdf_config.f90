!>
!! Configuration of the parameterization for vertical diffusion,
!! that is used in the ECHAM physics package.
!!
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!! First version by Marco Giorgetta, MPI-M (2017-04)
!!
!! Based on earlier codes of:
!!     ...
!!
!! References: 
!!     Angevine, W. M., Jiang, H., & Mauritsen T. (2010).
!!           Performance of an eddy diffusivity- mass flux scheme for shallow cumulus boundary layers.
!!           Monthly Weather Review, 138(7), 2895–2912. https://doi.org/10.1175/2010MWR3142.1
!!     Mauritsen, T., & Svensson, G. (2007).
!!           Observations of stably stratified shear-driven atmospheric turbulence at low and high Richardson numbers.
!!           Journal of the Atmospheric Sciences, 64(2), 645–655. https://doi.org/10.1175/JAS3856.1
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_vdf_config

  USE mo_exception            ,ONLY: finish, message, print_value
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom
  USE mo_grid_config          ,ONLY: n_dom

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         echam_vdf_config   !< user specified configuration parameters
  PUBLIC ::    init_echam_vdf_config   !< allocate and initialize echam_vdf_config
  PUBLIC ::    eval_echam_vdf_config   !< evaluate echam_vdf_config
  PUBLIC ::   print_echam_vdf_config   !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'echam_vdf'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the ECHAM physics package
  !!
  TYPE t_echam_vdf_config
     !
     ! configuration parameters
     ! ------------------------
     !
     LOGICAL  :: lsfc_mom_flux   !< switch on/off surface momentum flux
     LOGICAL  :: lsfc_heat_flux  !< switch on/off surface heat flux (sensible AND latent)
     !
     REAL(wp) :: pr0             !< neutral limit Prandtl number, can be varied from about 0.6 to 1.0
     REAL(wp) :: f_tau0          !< neutral non-dimensional stress factor
     REAL(wp) :: f_theta0        !< neutral non-dimensional heat flux factor 
     REAL(wp) :: c_f             !< mixing length: coriolis term tuning parameter
     REAL(wp) :: c_n             !< mixing length: stability term tuning parameter
     REAL(wp) :: c_e             !< dissipation coefficient (=f_tau0^(3/2))
     REAL(wp) :: wmc             !< ratio of typical horizontal velocity to wstar at free convection
     REAL(wp) :: fsl             !< fraction of first-level height at which surface fluxes
     !                              are nominally evaluated, tuning param for sfc stress
     REAL(wp) :: fbl             !< 1/fbl: fraction of BL height at which lmix hat its max
     !
     INTEGER :: turb             !< 1: TTE scheme, 2: 3D Smagorinsky
     REAL(wp) :: smag_constant
     REAL(wp) :: turb_prandtl
     REAL(wp) :: rturb_prandtl     !inverse turbulent prandtl number
     REAL(wp) :: km_min        !min mass weighted turbulent viscosity
     REAL(wp) :: max_turb_scale !max turbulence length scale
     REAL(wp) :: min_sfc_wind  !min sfc wind in free convection limit
     !
  END TYPE t_echam_vdf_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_echam_vdf_config), TARGET :: echam_vdf_config(max_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_echam_vdf_config
    !
    ! ECHAM VDF configuration
    ! -----------------------
    !
    echam_vdf_config(:)% lsfc_mom_flux  = .TRUE.
    echam_vdf_config(:)% lsfc_heat_flux = .TRUE.
    !
    echam_vdf_config(:)% pr0      =  1.0_wp
    echam_vdf_config(:)% f_tau0   =  0.17_wp
    echam_vdf_config(:)% f_theta0 = -SQRT(echam_vdf_config(:)%f_tau0**2/2.0_wp/echam_vdf_config(:)%pr0)
    echam_vdf_config(:)% c_f      =  0.185_wp
    echam_vdf_config(:)% c_n      =  2.0_wp
    echam_vdf_config(:)% c_e      =  SQRT(echam_vdf_config(:)%f_tau0**3)
    echam_vdf_config(:)% wmc      =  0.5_wp
    echam_vdf_config(:)% fsl      =  0.4_wp
    echam_vdf_config(:)% fbl      =  3._wp
    echam_vdf_config(:)%turb      =  1
    echam_vdf_config(:)%smag_constant = 0.23_wp
    echam_vdf_config(:)%max_turb_scale= 300._wp
    echam_vdf_config(:)%turb_prandtl  = 0.33333333333_wp
    echam_vdf_config(:)%rturb_prandtl = 1/echam_vdf_config(:)%turb_prandtl
    echam_vdf_config(:)%km_min        =  0.001_wp
    echam_vdf_config(:)%min_sfc_wind  =  1._wp
    !
  END SUBROUTINE init_echam_vdf_config

  !----

  !>
  !! Evaluate additional derived parameters
  !!
  SUBROUTINE eval_echam_vdf_config
    !
    CHARACTER(LEN=*), PARAMETER :: routine = 'eval_echam_vdf_config'

#ifdef _OPENACC
    IF ( ANY (echam_vdf_config(:)%turb == 2) ) THEN ! turb=2 is Smagorinsky
      CALL finish (routine, 'Smagorinsky not yet implemented with OpenACC')
    END IF
#endif
    !
    ! check range of pr0
    IF ( ANY(echam_vdf_config(:)% pr0 <0.6_wp ) ) CALL message(routine//' WARNING:',' echam_vdf_config(:)% pr0 < 0.6 ')
    IF ( ANY(echam_vdf_config(:)% pr0 >1.0_wp ) ) CALL message(routine//' WARNING:',' echam_vdf_config(:)% pr0 > 1.0 ')
    !
    ! re-compute dependent parameters 
    echam_vdf_config(:)% f_theta0 = -SQRT(echam_vdf_config(:)%f_tau0**2/2.0_wp/echam_vdf_config(:)%pr0)
    echam_vdf_config(:)% c_e      =  SQRT(echam_vdf_config(:)%f_tau0**3)
    !
  END SUBROUTINE eval_echam_vdf_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_echam_vdf_config
    !
    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','ECHAM vertical diffusion configuration')
    CALL message    ('','======================================')
    CALL message    ('','')
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% lsfc_mom_flux  ',echam_vdf_config(jg)% lsfc_mom_flux  )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% lsfc_heat_flux ',echam_vdf_config(jg)% lsfc_heat_flux )
       CALL message    ('','')
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% pr0            ',echam_vdf_config(jg)% pr0            )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% f_tau0         ',echam_vdf_config(jg)% f_tau0         )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% f_theta0       ',echam_vdf_config(jg)% f_theta0       )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% c_f            ',echam_vdf_config(jg)% c_f            )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% c_n            ',echam_vdf_config(jg)% c_n            )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% c_e            ',echam_vdf_config(jg)% c_e            )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% wmc            ',echam_vdf_config(jg)% wmc            )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% fsl            ',echam_vdf_config(jg)% fsl            )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% fbl            ',echam_vdf_config(jg)% fbl            )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% turb           ',echam_vdf_config(jg)% turb           )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% smag_constant  ',echam_vdf_config(jg)% smag_constant  )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% turb_prandtl   ',echam_vdf_config(jg)% turb_prandtl   )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% rturb_prandtl  ',echam_vdf_config(jg)% rturb_prandtl  )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% km_min         ',echam_vdf_config(jg)% km_min         )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% max_turb_scale ',echam_vdf_config(jg)% max_turb_scale )
       CALL print_value('    echam_vdf_config('//TRIM(cg)//')% min_sfc_wind   ',echam_vdf_config(jg)% min_sfc_wind   )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_echam_vdf_config

  !----

END MODULE mo_echam_vdf_config
