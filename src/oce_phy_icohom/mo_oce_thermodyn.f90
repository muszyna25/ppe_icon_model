!>
!! Provide an implementation of the ocean thermodynamics.
!!
!! Provide an implementation of the parameters used for the thermodynamics
!! of the hydrostatic ocean model.
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
!!  Modified by Stephan Lorenz,     MPI-M, 2010-07
!!   - adapted to structures discussed in 2010-01.
!!  Modified by Stephan Lorenz,     MPI-M, 2010-10
!!   - structured input/output parameters
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_oce_thermodyn
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2007
!
!-------------------------------------------------------------------------
!
!
!
!
USE mo_kind,                ONLY: wp
USE mo_ocean_nml,           ONLY: n_zlev, t_ref, s_ref
USE mo_model_domain,        ONLY: t_patch
USE mo_impl_constants,      ONLY: min_rlcell,sea_boundary, sea_boundary, &
  &                               toplev
! &                               success, max_char_length,                &
!USE mo_exception,           ONLY: message, finish
USE mo_loopindices,         ONLY: get_indices_c!, get_indices_e, get_indices_v
USE mo_oce_physics,          ONLY: t_ho_params!params_oce
USE mo_physical_constants,  ONLY: grav
!USE mo_model_domain_import, ONLY: ??

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

CHARACTER(len=*), PARAMETER :: this_mod_name = 'mo_oce_thermodyn'

PUBLIC  :: calc_internal_press
PUBLIC  :: calc_density
PRIVATE :: calc_density_lin_EOS
PRIVATE :: calc_density_JM_EOS
CONTAINS
  !
  !-------------------------------------------------------------------------
  !
  !>
  !! Calculation the hydrostatic pressure
  !!
  !! Calculation the hydrostatic pressure by computing the weight of the
  !! fluid column above a certain level.
  !!
  !!
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2009)
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !! Modified by Stephan Lorenz,        MPI-M (2010-10-22)
  !!  - division by rho_ref included
  !!
  SUBROUTINE calc_internal_press(ppatch, p_phys_param, rho, h, press_hyd)
  !
  TYPE(t_patch), INTENT(IN)     :: ppatch
  TYPE (t_ho_params),INTENT(IN) :: p_phys_param
  REAL(wp),    INTENT(IN)       :: rho      (:,:,:)  !< density
  REAL(wp),    INTENT(IN)       :: h        (:,:)    !< surface elevation at cells
  REAL(wp),   INTENT(INOUT)     :: press_hyd(:,:,:)  !< hydrostatic pressure

  ! local variables:

  !CHARACTER(len=max_char_length), PARAMETER :: &
  !       & routine = (this_mod_name//':calc_internal_pressure')

  INTEGER :: slev, end_lev     ! vertical start and end level
  INTEGER :: jc, jk, jb
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

  REAL(wp) :: z_full, z_box!, max_den, min_den, max_pres, min_pres!, p_hyd(4)

  !-------------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')
  ! #slo# due to nag -nan compiler-option set intent(out) variables to zero
  !press_hyd(:,:,:) = 0.0_wp

  rl_start = 1
  rl_end = min_rlcell
  slev = 1

  ! values for the blocking
  i_startblk = ppatch%cells%start_blk(rl_start,1)
  i_endblk  = ppatch%cells%end_blk(rl_end,1)
  !
  !  loop through all patch cells
  !
  DO jb = i_startblk, i_endblk
    CALL get_indices_c( ppatch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

    DO jc = i_startidx, i_endidx
      !
      !  #slo# calculation of pressure due to elevation is done here
      !  including actual density of surface water
      !
      !z_full = grav * rho(jc,toplev,jb) * h(jc,jb)
      !  #slo# 2011-01-19 - elevation not considered:
      !   - in SWM ok, since density is constant
      !   - check to include h if tracers (T, S) are active
      z_full  = 0.0_wp
      end_lev = ppatch%patch_oce%dolic_c(jc,jb)

      DO jk = slev, end_lev
        IF(ppatch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN   
          z_box               = grav*(rho(jc,jk,jb)-p_phys_param%rho_ref)&         !pressure in single box at layer jk
                                    &*ppatch%patch_oce%del_zlev_m(jk) 
          press_hyd(jc,jk,jb) = ( z_full + 0.5_wp*z_box ) * p_phys_param%rho_inv   !hydrostatic press at level jk
                                                                                   ! =half of pressure at actual box+ sum of all boxes above
          z_full              = z_full + z_box
        ELSE
          press_hyd(jc,jk,jb) = 0.0_wp
       ENDIF
      END DO

!       !the code within the level loop below is identical to the level loop above.
!       !It is a bit more transparentr 
!       p_hyd(:) = 0.0_wp
!       p_hyd(1) = grav*(rho(jc,1,jb)-p_phys_param%rho_ref)*ppatch%patch_oce%del_zlev_m(1)* p_phys_param%rho_inv*0.5_wp
!       DO jk = slev+1, end_lev
!         IF(ppatch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN   
! 
!           z_box = (rho(jc,jk-1,jb)-p_phys_param%rho_ref)*ppatch%patch_oce%del_zlev_m(jk-1)&
!                &+ (rho(jc,jk,jb)  -p_phys_param%rho_ref)*ppatch%patch_oce%del_zlev_m(jk)
!           p_hyd(jk) = p_hyd(jk-1)&
!                    &+ p_phys_param%rho_inv*grav*0.5_wp*z_box
!         ELSE
!           press_hyd(jc,jk,jb) = 0.0_wp
!        ENDIF
!      END DO 
! DO jk = slev, end_lev
!  IF(ppatch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN 
! ! IF(jk==1)THEN
!  write(*,*)'pressure',jk,jc,jb,press_hyd(jc,jk,jb),p_hyd(jk)
! ! ENDIF
!  ENDIF
! END DO
! ! IF(jb==900.and. jk==3)THEN
! ! write(*,*)'pressure sample',jc,jk,jb, press_hyd(jc,jk,jb), rho(jc,jk,jb)
! ! ENDIF
    END DO
  END DO
!  DO jk = 1, n_zlev
!   max_den=0.0_wp
!   min_den=0.0_wp
!   max_pres=0.0_wp
!   min_pres=0.0_wp
!    DO jb = i_startblk, i_endblk
!      CALL get_indices_c( ppatch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
!       DO jc = i_startidx, i_endidx
!         IF(ppatch%patch_oce%lsm_oce_c(jc,jk,jb) == sea ) THEN   
!           IF(max_pres> press_hyd(jc,jk,jb))max_pres= press_hyd(jc,jk,jb)
!           IF(min_pres< press_hyd(jc,jk,jb))min_pres= press_hyd(jc,jk,jb)
! 
!           IF(max_den> rho(jc,jk,jb))max_den= rho(jc,jk,jb)
!           IF(min_den<= rho(jc,jk,jb))min_den= rho(jc,jk,jb)
!       ENDIF
!     END DO
!   END DO
! write(*,*)'max-min dens:', jk,max_pres, min_pres, max_den, min_den
! END DO
  !CALL message (TRIM(routine), 'end')
  END SUBROUTINE calc_internal_press

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Calculates the density via a call to the equation-of-state.
  !! Several options for EOS are provided.
  !!
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2009)
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !!
  SUBROUTINE calc_density(ppatch,tracer, rho)
  !
  !!
  TYPE(t_patch), INTENT(IN) :: ppatch
  TYPE (t_ho_params)        :: p_param
  REAL(wp),    INTENT(IN)   :: tracer(:,:,:,:)     !< input of S and T
  REAL(wp), INTENT(INOUT)   :: rho   (:,:,:)       !< density

  ! local variables:
   INTEGER :: EOS_TYPE = 1
  ! CHARACTER(len=max_char_length), PARAMETER :: &
  !      & routine = (this_mod_name//':calc_density')
  !---------------------------------------------------------------------

  ! CALL message (TRIM(routine), 'start')

   SELECT CASE (EOS_TYPE)
     CASE(1)
      CALL calc_density_lin_EOS(ppatch, p_param, tracer, rho)
     CASE(2)
       CALL calc_density_JM_EOS(ppatch, p_param, tracer, rho)!(ppatch, params_oce, tracer, rho)
  !   CASE(3)
  !     CALL calc_density_JMDWFG06_EOS(ppatch, pprog, pdiag)
     CASE DEFAULT

   END SELECT

  END SUBROUTINE calc_density
  !-------------------------------------------------------------------------
  !
  !
  !>
  !!Calculates the density via a linear equation-of-state.
  !!
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2009)
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !!
  SUBROUTINE  calc_density_lin_EOS(ppatch, params_oce, tracer, rho)
  !
  !!
  !!
  TYPE(t_patch), INTENT(IN)     :: ppatch
  TYPE (t_ho_params),INTENT(IN) :: params_oce
  REAL(wp),    INTENT(IN)       :: tracer(:,:,:,:)     !< input of S and T
  REAL(wp), INTENT(INOUT)       :: rho   (:,:,:)       !< density

  ! local variables:

  INTEGER :: jc, jk, jb
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
  !-------------------------------------------------------------------------
  rl_start   = 1
  rl_end     = min_rlcell
  i_startblk = ppatch%cells%start_blk(rl_start,1)
  i_endblk   = ppatch%cells%end_blk(rl_end,1)

  DO jb = i_startblk, i_endblk
    CALL get_indices_c( ppatch, jb, i_startblk, i_endblk,&
                       & i_startidx, i_endidx, rl_start, rl_end)
    !  tracer 1: potential temperature
    !  tracer 2: salinity
    DO jk=1, n_zlev
      DO jc = i_startidx, i_endidx
        IF(ppatch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN   
!         rho(jc,jk,jb) =  ( params_oce%rho_ref                              +  &
!           &                params_oce%a_T * (tracer(jc,jk,jb,1) - t_ref)   +  &
!           &                params_oce%b_S * (tracer(jc,jk,jb,2) - s_ref) )
          rho(jc,jk,jb) =    params_oce%rho_ref                    &
            &              - params_oce%a_T * tracer(jc,jk,jb,1)   &
            &              + params_oce%b_S * tracer(jc,jk,jb,2)
         ELSE
           rho(jc,jk,jb) = 0.0_wp 
        ENDIF
!          IF(jk==1)THEN
!         IF(ppatch%patch_oce%lsm_oce_c(jc,jk,jb) == sea ) THEN   
!          write(*,*)'density',jk,jc,jb,params_oce%rho_ref, tracer(jc,jk,jb,1),&
!           &rho(jc,jk,jb)
!         ENDIF
!         ENDIF
      END DO
    END DO
  END DO

  END SUBROUTINE  calc_density_lin_EOS
  !-------------------------------------------------------------------------
 !
  !
  !>
  !!  Calculates density as a function of potential temperature and salinity
  !! using the Jackett and McDougall equation of state
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2011)
  !! Code below is an adaption of Sergey Danilov's implementation in
  !! the AWI Finitre-Volume model. 
  !!
  SUBROUTINE calc_density_JM_EOS(ppatch,params_oce, tracer, rho)
  !
  !!
  !!
  TYPE(t_patch), INTENT(IN)     :: ppatch
  TYPE (t_ho_params),INTENT(IN) :: params_oce
  REAL(wp),    INTENT(IN)       :: tracer(:,:,:,:)  
  REAL(wp), INTENT(OUT)         :: rho(:,:,:) 

  ! local variables:              
  REAL(wp) :: z_t
  REAL(wp) :: z_s
  REAL(wp) :: z_rhopot, z_bulk, pz 
  !INTEGER  :: slev, end_lev
  INTEGER  :: jc, jk, jb
  INTEGER  :: rl_start, rl_end
  INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
  !---------------------------------------------------------------------------
  rl_start   = 1
  rl_end     = min_rlcell
  i_startblk = ppatch%cells%start_blk(rl_start,1)
  i_endblk   = ppatch%cells%end_blk(rl_end,1)

  DO jb = i_startblk, i_endblk
    CALL get_indices_c( ppatch, jb, i_startblk, i_endblk,&
                       & i_startidx, i_endidx, rl_start, rl_end)

    DO jk=1, n_zlev
      DO jc = i_startidx, i_endidx

      pz=ppatch%patch_oce%zlev_m(jk)
      z_t = tracer(jc,jk,jb,1)
      z_s = tracer(jc,jk,jb,2)

      !compute secant bulk modulus
      z_bulk = 19092.56_wp + z_t*(209.8925_wp             &
       - z_t*(3.041638_wp - z_t*(-1.852732e-3_wp          &
       - z_t*(1.361629e-5_wp))))                          &
       + z_s*(104.4077_wp - z_t*(6.500517_wp              &
       - z_t*(.1553190_wp - z_t*(-2.326469e-4_wp))))      &
       + sqrt(z_s**3)*(-5.587545_wp                       &
       + z_t*(0.7390729_wp - z_t*(1.909078e-2_wp)))       &
       - pz *(4.721788e-1_wp + z_t*(1.028859e-2_wp        &
       + z_t*(-2.512549e-4_wp - z_t*(5.939910e-7_wp))))   &
       - pz*z_s*(-1.571896e-2_wp                          &
       - z_t*(2.598241e-4_wp + z_t*(-7.267926e-6_wp)))    &
       - pz*sqrt(z_s**3)                                  &
       *2.042967e-3_wp + pz*pz*(1.045941e-5_wp            &
       - z_t*(5.782165e-10_wp - z_t*(1.296821e-7_wp)))    &
       + pz*pz*z_s                                        &
       *(-2.595994e-7_wp                                  &
       + z_t*(-1.248266e-9_wp + z_t*(-3.508914e-9_wp)))

      z_rhopot = ( 999.842594_wp                     &
       + z_t*( 6.793952e-2_wp                        &
       + z_t*(-9.095290e-3_wp                        &
       + z_t*( 1.001685e-4_wp                        &
       + z_t*(-1.120083e-6_wp                        &
       + z_t*( 6.536332e-9_wp)))))                   &
       + z_s*( 0.824493_wp                           &
       + z_t *(-4.08990e-3_wp                        &
       + z_t *( 7.64380e-5_wp                        &
       + z_t *(-8.24670e-7_wp                        &
       + z_t *( 5.38750e-9_wp)))))                   &
       + sqrt(z_s**3)*(-5.72466e-3_wp                &
       + z_t*( 1.02270e-4_wp                         &
       + z_t*(-1.65460e-6_wp)))                      &
       + 4.8314e-4_wp*z_s**2)

       rho(jc,jk,jb) = z_rhopot/(1.0_wp + 0.1_wp*pz/z_bulk)&
                     & - params_oce%rho_ref
      END DO
    END DO
  END DO

end subroutine calc_density_JM_EOS

END MODULE mo_oce_thermodyn

