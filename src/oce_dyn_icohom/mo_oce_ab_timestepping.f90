!>
!! Contains the implementation of the semi-implicit Adams-Bashforth timestepping
!! for the ICON ocean model.
!! 
!! 
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2010/04)
!!  Modified by Stephan Lorenz,     MPI-M (2010-06)
!!   - renaming and adjustment to ocean domain and patch
!!   - implementation of continuity equation for vertical velocities
!! 
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_oce_ab_timestepping
!-------------------------------------------------------------------------  
!
!    ProTeX FORTRAN source: Style 2  
!    modified for ICON project, DWD/MPI-M 2006
!  
!-------------------------------------------------------------------------  
!  
!   
! 
!USE mo_kind,                           ONLY: wp
!USE mo_mpi,                            ONLY: p_pe, p_io
USE mo_impl_constants,                 ONLY: max_char_length!sea_boundary, &
!  &                                          min_rlcell, min_rledge, min_rlcell, &
USE mo_ocean_nml,                      ONLY: idisc_scheme
USE mo_dynamics_nml,                   ONLY: nnew, nold
USE mo_oce_state,                      ONLY: t_hydro_ocean_state!, t_hydro_ocean_diag
USE mo_oce_forcing,                     ONLY: t_ho_sfc_flx
USE mo_interpolation,                  ONLY: t_int_state
USE mo_model_domain,                   ONLY: t_patch
!USE mo_exception,                      ONLY: message, finish!, message_text
!USE mo_loopindices,                    ONLY: get_indices_c, get_indices_e !, get_indices_v
!USE mo_oce_index,                      ONLY: c_i, c_b, c_k, ne_b, ne_i, nc_b, nc_i, form4ar, ldbg
USE mo_oce_ab_timestepping_mimetic,   ONLY: solve_free_sfc_ab_mimetic,&
                                           & calc_normal_velocity_ab_mimetic, &
                                           & calc_vert_velocity_mimetic
USE mo_oce_ab_timestepping_rbf,     ONLY: solve_free_sfc_ab_RBF,    &
                                           & calc_normal_velocity_ab_RBF,     &
                                           & calc_vert_velocity_RBF
USE mo_oce_physics,                ONLY: t_ho_params
IMPLICIT NONE

PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'
INTEGER, PARAMETER :: MIMETIC_TYPE = 1
INTEGER, PARAMETER :: RBF_TYPE     = 2
!
! PUBLIC INTERFACE
!
PUBLIC :: solve_free_surface_eq_ab
PUBLIC :: calc_normal_velocity_ab
PUBLIC :: calc_vert_velocity

CONTAINS
!-------------------------------------------------------------------------  
!
!  
!>
!! !  Solves the free surface equation.
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
SUBROUTINE solve_free_surface_eq_ab(p_patch, p_os, p_sfc_flx, p_phys_param, timestep, p_int)
!
TYPE(t_patch), TARGET, INTENT(in)             :: p_patch
TYPE(t_hydro_ocean_state), TARGET             :: p_os 
TYPE(t_ho_sfc_flx), INTENT(INOUT)             :: p_sfc_flx    
TYPE (t_ho_params)                            :: p_phys_param
INTEGER                                       :: timestep
TYPE(t_int_state),TARGET,INTENT(IN), OPTIONAL :: p_int
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_ab_timestepping:solve_free_surface_eq_2tl_ab')
!-------------------------------------------------------------------------------

IF(idisc_scheme==MIMETIC_TYPE)THEN

  CALL solve_free_sfc_ab_mimetic(p_patch, p_os, p_sfc_flx, p_phys_param, timestep, p_int)

ELSEIF(idisc_scheme==RBF_TYPE)THEN

  CALL solve_free_sfc_ab_RBF(p_patch, p_os, p_sfc_flx, p_phys_param, timestep, p_int)

ENDIF

!CALL message (TRIM(routine), 'end')
END SUBROUTINE solve_free_surface_eq_ab
!-------------------------------------------------------------------------  
!
!  
!>
!! Computation of new velocity in Adams-Bashforth timestepping.
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
SUBROUTINE calc_normal_velocity_ab(p_patch, p_os)
!
! Patch on which computation is performed
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
!
! Type containing ocean state
TYPE(t_hydro_ocean_state), TARGET :: p_os
!
!
!  local variables
!
! CHARACTER(len=max_char_length), PARAMETER ::     &
!   &      routine = ('mo_oce_ab_timestepping: calc_normal_velocity_2tl_ab')
!-----------------------------------------------------------------------  
IF(idisc_scheme==MIMETIC_TYPE)THEN

  CALL calc_normal_velocity_ab_mimetic(p_patch, p_os)

ELSEIF(idisc_scheme==RBF_TYPE)THEN

  CALL calc_normal_velocity_ab_RBF(p_patch, p_os)

ENDIF

END SUBROUTINE calc_normal_velocity_ab
!-------------------------------------------------------------------------  
!
!  
!>
!! Computation of new vertical velocity using continuity equation

!! Calculate diagnostic vertical velocity from horizontal velocity using the
!! incommpressibility condition in the continuity equation.
!! For the case of the semi-implicit-AB scheme the land-sea-mask may be applied
!! at least after collecting the whole explicit term.
!! 
!! @par Revision History
!! Developed  by  Peter Korn,   MPI-M (2006).
!! 
SUBROUTINE calc_vert_velocity( p_patch, p_os)
!
TYPE(t_patch), TARGET, INTENT(IN) :: p_patch       ! patch on which computation is performed
TYPE(t_hydro_ocean_state)         :: p_os
!
!
! Local variables
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_ab_timestepping:calc_vert_velocity')
!-----------------------------------------------------------------------  

!Store current vertical velocity before the new one is calculated
p_os%p_diag%w_old = p_os%p_diag%w

IF(idisc_scheme==MIMETIC_TYPE)THEN

  CALL calc_vert_velocity_mimetic( p_patch,            &
                             & p_os,                   &
                             & p_os%p_diag,            &
                             & p_os%p_prog(nnew(1))%h, &
                             & p_os%p_aux%bc_top_w,    &
                             & p_os%p_aux%bc_bot_w,    &
                             & p_os%p_diag%w )
!  CALL calc_vert_velocity_RBF( p_patch,&
!                              & p_os%p_prog(nnew(1))%vn,&
!                              & p_os%p_prog(nold(1))%vn,&
!                              & p_os%p_prog(nnew(1))%h, &
!                              & p_os%p_aux%bc_top_w,    &
!                              & p_os%p_aux%bc_bot_w,    &
!                              & p_os%p_diag%w )
ELSEIF(idisc_scheme==RBF_TYPE)THEN

  CALL calc_vert_velocity_RBF( p_patch,&
                             & p_os%p_prog(nnew(1))%vn,&
                             & p_os%p_prog(nold(1))%vn,&
                             & p_os%p_prog(nnew(1))%h, &
                             & p_os%p_diag%w(:,1,:),    &
                             & p_os%p_aux%bc_bot_w,    &
!                              & p_os%p_aux%bc_top_w,    &
!                              & p_os%p_aux%bc_bot_w,    &
                             & p_os%p_diag%w )
ENDIF


END SUBROUTINE calc_vert_velocity
!-------------------------------------------------------------------------  


END MODULE mo_oce_ab_timestepping
