!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_diffusion_config

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text
  USE mo_impl_constants,      ONLY: max_dom

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_diffusion_config, diffusion_config  !< derived type and variable
  PUBLIC :: configure_diffusion                   !< subroutine

  !--------------------------------------------------------------------------
  ! Basic configuration setup for diffusion
  !--------------------------------------------------------------------------
  TYPE t_diffusion_config

    ! variables from namelist

    INTEGER :: hdiff_order  ! order of horizontal diffusion
                            ! -1: no diffusion
                            ! 2: 2nd order linear diffusion on all vertical levels 
                            ! 3: Smagorinsky diffusion for hexagonal model
                            ! 4: 4th order linear diffusion on all vertical levels 
                            ! 5: Smagorinsky diffusion for triangular model
                            

    REAL(wp) :: hdiff_efdt_ratio      ! ratio of e-folding time to (2*)time step
    REAL(wp) :: hdiff_w_efdt_ratio    ! ratio of e-folding time to time step for w diffusion (NH only)
    REAL(wp) :: hdiff_min_efdt_ratio  ! minimum value of hdiff_efdt_ratio 
                                      ! (for upper sponge layer)
    REAL(wp) :: hdiff_tv_ratio        ! the ratio of diffusion coefficient: temp:mom
    REAL(wp) :: hdiff_smag_fac        ! scaling factor for Smagorinsky diffusion
    REAL(wp) :: hdiff_multfac         ! multiplication factor of normalized diffusion
                                      ! coefficient for nested domains
    INTEGER  :: itype_vn_diffu        ! options for discretizing the Smagorinsky momentum diffusion
    INTEGER  :: itype_t_diffu         ! options for discretizing the Smagorinsky temperature diffusion

    LOGICAL :: lhdiff_temp   ! if .TRUE., apply horizontal diffusion to temp.
    LOGICAL :: lhdiff_vn     ! if .TRUE., apply horizontal diffusion to momentum.
    LOGICAL :: lhdiff_w      ! if .TRUE., apply horizontal diffusion to vertical momentum.
    LOGICAL :: lsmag_3d      ! if .TRUE., compute 3D Smagorinsky diffusion coefficient.

    ! variables not from namelist

    REAL(wp) :: k6, k4, k2, k4w  ! numerical diffusion coefficients
                                 ! Values for these parameters are not directly
                                 ! specified by the user, but derived from the ratio 
                                 ! between the e-folding time and the model time step
                                 ! (hdiff_efdt_ratio above), and the horizontal 
                                 ! resolution of the model


  END TYPE t_diffusion_config
  !>
  !!
  TYPE(t_diffusion_config),SAVE :: diffusion_config(max_dom)

CONTAINS
  !>
  !!
  SUBROUTINE configure_diffusion( n_dom, dynamics_parent_grid_id )

    INTEGER, INTENT(IN) :: n_dom
    INTEGER, INTENT(IN) :: dynamics_parent_grid_id(max_dom)

    INTEGER  :: jg, jk, jgp
    REAL(wp) :: tmp

    CHARACTER(len=*), PARAMETER :: &
      routine = 'mo_diffusion_config:configure_diffusion'


    !-----------------------------------------------------------
    ! Compute diffusion coefficients
    !-----------------------------------------------------------

    IF (ANY(diffusion_config(1:n_dom)%hdiff_efdt_ratio<=0._wp)) THEN

      diffusion_config(:)%k2 = 0._wp
      diffusion_config(:)%k4 = 0._wp
      diffusion_config(:)%k6 = 0._wp

      CALL message(TRIM(routine),'Background linear diffusion is '//&
                                 'switched off in all domains')

    ELSE

      tmp = 1._wp/(diffusion_config(1)%hdiff_efdt_ratio*8._wp)
      diffusion_config(1)%k2 = tmp
      tmp = 1._wp/(diffusion_config(1)%hdiff_efdt_ratio*64._wp)
      diffusion_config(1)%k4 = tmp
      tmp = 1._wp/(diffusion_config(1)%hdiff_efdt_ratio*512._wp)
      diffusion_config(1)%k6 = tmp
  
      tmp = 1._wp/(diffusion_config(1)%hdiff_w_efdt_ratio*36._wp)
      diffusion_config(:)%k4w = tmp

      DO jg = 2, n_dom

         jgp = dynamics_parent_grid_id(jg)
         
         diffusion_config(jg)%k2 = diffusion_config(jgp)%k2 * diffusion_config(jg)%hdiff_multfac
         diffusion_config(jg)%k4 = diffusion_config(jgp)%k4 * diffusion_config(jg)%hdiff_multfac
         diffusion_config(jg)%k6 = diffusion_config(jgp)%k6 * diffusion_config(jg)%hdiff_multfac
      ENDDO

    ENDIF

  END SUBROUTINE configure_diffusion

END MODULE mo_diffusion_config
