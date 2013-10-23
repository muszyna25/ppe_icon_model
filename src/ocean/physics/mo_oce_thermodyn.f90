
!---------------------------------------------------------------------------------------
!>
!! Provide an implementation of the ocean thermodynamics
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
!!  mpi parallelized LL (no sync required)
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
  USE mo_kind,                ONLY: wp
  USE mo_ocean_nml,           ONLY: n_zlev, eos_type, no_tracer
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_impl_constants,      ONLY: sea_boundary, sea_boundary, min_dolic !, &
  !USE mo_exception,           ONLY: message, finish
  USE mo_loopindices,         ONLY: get_indices_c!, get_indices_e, get_indices_v
  USE mo_physical_constants,  ONLY: grav, rho_ref, sal_ref, rho_inv, a_t, b_s, &
    & sitodbar, sfc_press_bar
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_parallel_config,     ONLY: nproma
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_oce_thermodyn'
  
  PUBLIC :: calc_internal_press
  PUBLIC :: calc_internal_press_new
  PUBLIC :: calc_density,calc_potential_density
  !each specific EOS comes as a sbr and as a function. The sbr version is private as it is
  !only used in "calc_internal_press", whilethe function version is used in mo_oce_physics
  !(sbr "update_ho_params") to calculate the local Richardson number.
  ! PRIVATE :: calc_density_lin_eos
  PUBLIC :: calc_density_lin_eos_func
  ! PRIVATE :: calc_density_jmdwfg06_eos
  PUBLIC :: calc_density_jmdwfg06_eos_func
  PUBLIC :: calc_density_mpiom_func
  PUBLIC :: calc_density_mpiom_elemental
  ! PRIVATE :: calc_density_mpiom
  ! PRIVATE :: convert_insitu2pot_temp
  PUBLIC :: convert_insitu2pot_temp_func
  PUBLIC :: adisit
  REAL(wp), PARAMETER :: eosmdjwfnum(0:11) = (/                                 &
    & 9.99843699e+02_wp,  7.35212840e+00_wp, -5.45928211e-02_wp,                 &
    & 3.98476704e-04_wp,  2.96938239e+00_wp, -7.23268813e-03_wp,                 &
    & 2.12382341e-03_wp,  1.04004591e-02_wp,  1.03970529e-07_wp,                 &
    & 5.18761880e-06_wp, -3.24041825e-08_wp, -1.23869360e-11_wp /)
  
  REAL(wp), PARAMETER :: eosmdjwfden(0:12) = (/                                 &
    & 1.00000000e+00_wp,  7.28606739e-03_wp, -4.60835542e-05_wp,                  &
    & 3.68390573e-07_wp,  1.80809186e-10_wp,  2.14691708e-03_wp,                  &
    & -9.27062484e-06_wp, -1.78343643e-10_wp,  4.76534122e-06_wp,                  &
    & 1.63410736e-09_wp,  5.30848875e-06_wp, -3.03175128e-16_wp,                  &
    & -1.27934137e-17_wp /)
  
  ! ! constants used within the Jackett et al. (2006) nonlinear equation of
  !   state
  
  REAL(wp), PARAMETER :: eosjmdwfgnum(0:11) = (/                                &
    & 9.9984085444849347e+02_wp,  7.3471625860981584e+00_wp,                     &
    & -5.3211231792841769e-02_wp,  3.6492439109814549e-04_wp,                     &
    & 2.5880571023991390e+00_wp, -6.7168282786692355e-03_wp,                     &
    & 1.9203202055760151e-03_wp,  1.1798263740430364e-02_wp,                     &
    & 9.8920219266399117e-08_wp,  4.6996642771754730e-06_wp,                     &
    & -2.5862187075154352e-08_wp, -3.2921414007960662e-12_wp /)
  
  REAL(wp), PARAMETER :: eosjmdwfgden(0:12) = (/    1.0_wp,                     &
    & 7.2815210113327091e-03_wp, -4.4787265461983921e-05_wp,                     &
    & 3.3851002965802430e-07_wp,  1.3651202389758572e-10_wp,                     &
    & 1.7632126669040377e-03_wp, -8.8066583251206474e-06_wp,                     &
    & -1.8832689434804897e-10_wp,  5.7463776745432097e-06_wp,                     &
    & 1.4716275472242334e-09_wp,  6.7103246285651894e-06_wp,                     &
    & -2.4461698007024582e-17_wp, -9.1534417604289062e-18_wp /)
  
  REAL (wp), PARAMETER ::  dbl_eps   = EPSILON(1._wp)
  
  
  REAL(wp), PARAMETER :: &
    & a_a1=3.6504E-4_wp, a_a2=8.3198E-5_wp, a_a3=5.4065E-7_wp, &
    & a_a4=4.0274E-9_wp, &
    & a_b1=1.7439E-5_wp, a_b2=2.9778E-7_wp, &
    & a_c1=8.9309E-7_wp, a_c2=3.1628E-8_wp, a_c3=2.1987E-10_wp, &
    & a_d=4.1057E-9_wp, &
    & a_e1=1.6056E-10_wp, a_e2=5.0484E-12_wp
  
  REAL(wp), PARAMETER :: &
    & r_a0=999.842594_wp, r_a1=6.793952e-2_wp, r_a2=-9.095290e-3_wp, &
    & r_a3=1.001685e-4_wp, r_a4=-1.120083e-6_wp, r_a5=6.536332e-9_wp, &
    & r_b0=8.24493e-1_wp, r_b1=-4.0899e-3_wp, r_b2=7.6438e-5_wp, &
    & r_b3=-8.2467e-7_wp, r_b4=5.3875e-9_wp, &
    & r_c0=-5.72466e-3_wp, r_c1=1.0227e-4_wp, r_c2=-1.6546e-6_wp, &
    & r_d0=4.8314e-4_wp, &
    & r_e0=19652.21_wp, r_e1=148.4206_wp, r_e2=-2.327105_wp, &
    & r_e3=1.360477e-2_wp, r_e4=-5.155288e-5_wp, &
    & r_f0=54.6746_wp, r_f1=-0.603459_wp, r_f2=1.09987e-2_wp, &
    & r_f3=-6.1670e-5_wp, &
    & r_g0=7.944e-2_wp, r_g1=1.6483e-2_wp, r_g2=-5.3009e-4_wp, &
    & r_h0=3.239908_wp, r_h1=1.43713e-3_wp, r_h2=1.16092e-4_wp, &
    & r_h3=-5.77905e-7_wp, &
    & r_ai0=2.2838e-3_wp, r_ai1=-1.0981e-5_wp, r_ai2=-1.6078e-6_wp, &
    & r_aj0=1.91075e-4_wp, &
    & r_ak0=8.50935e-5_wp, r_ak1=-6.12293e-6_wp, r_ak2=5.2787e-8_wp, &
    & r_am0=-9.9348e-7_wp, r_am1=2.0816e-8_wp, r_am2=9.1697e-10_wp
  
CONTAINS

  !-------------------------------------------------------------------------
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
  SUBROUTINE calc_internal_press_new(p_patch_3d, trac_t, trac_s, h, calc_density_func, press_hyd)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in):: p_patch_3d
    REAL(wp),    INTENT(in)       :: trac_t   (:,:,:)  !temperature
    REAL(wp),    INTENT(in)       :: trac_s   (:,:,:)  !salinity
    REAL(wp),    INTENT(in)       :: h        (:,:)    !< surface elevation at cells
    REAL(wp),   INTENT(inout)     :: press_hyd(:,:,:)  !< hydrostatic pressure
    INTERFACE !This contains the function version of the actual EOS as chosen in namelist
      FUNCTION calc_density_func(tpot, sal, press) result(rho)
        USE mo_kind, ONLY: wp
        REAL(wp), INTENT(in) :: tpot
        REAL(wp), INTENT(in) :: sal
        REAL(wp), INTENT(in) :: press
        REAL(wp) :: rho
      ENDFUNCTION calc_density_func
    END INTERFACE
    ! local variables:
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !       & routine = (this_mod_name//':calc_internal_pressure')
    INTEGER :: slev, end_lev     ! vertical start and end level
    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    REAL(wp) :: z_full, z_box, z_press, z_rho_up, z_rho_down
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    ! #slo# due to nag -nan compiler-option set intent(inout) variables to zero
    !press_hyd(:,:,:) = 0.0_wp
    all_cells => p_patch%cells%ALL
    
    slev = 1
    press_hyd    = 0.0_wp
    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
      
      DO jc = i_startidx, i_endidx
        
        z_press      = (p_patch_3d%p_patch_1d(1)%zlev_i(1)+h(jc,jb))*rho_ref*sitodbar ! grav
        z_rho_up = calc_density_func(&
          & trac_t(jc,1,jb),&
          & trac_s(jc,1,jb),&
          & z_press)
        
        press_hyd(jc,slev,jb) = grav*z_rho_up*p_patch_3d%p_patch_1d(1)%del_zlev_m(1)*rho_inv!*0.5_wp
        
        ! write(*,*)'press',jc,jb,1,&
        !  &press_hyd(jc,1,jb), z_press
        
        end_lev = p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
        DO jk = slev+1, end_lev
          
          z_press = p_patch_3d%p_patch_1d(1)%zlev_i(jk)*rho_ref*sitodbar!grav
          !density of upper cell w.r.t.to pressure at intermediate level
          z_rho_up = calc_density_func(&
            & trac_t(jc,jk-1,jb),&
            & trac_s(jc,jk-1,jb),&
            & z_press)
          !density of lower cell w.r.t.to pressure at intermediate level
          z_rho_down = calc_density_func(&
            & trac_t(jc,jk,jb),&
            & trac_s(jc,jk,jb),&
            & z_press)
          
          z_box = ( z_rho_up*p_patch_3d%p_patch_1d(1)%del_zlev_m(jk-1)&
            & + z_rho_down*p_patch_3d%p_patch_1d(1)%del_zlev_m(jk))&
            & /(p_patch_3d%p_patch_1d(1)%del_zlev_m(jk)+p_patch_3d%p_patch_1d(1)%del_zlev_m(jk-1))
          press_hyd(jc,jk,jb) = press_hyd(jc,jk-1,jb) + rho_inv*grav*z_box
          !  write(*,*)'press',jc,jb,jk,&
          !  &press_hyd(jc,jk,jb), z_rho_up, z_rho_down,z_press
        END DO
        ! DO jk = slev, end_lev
        !  IF(v_base%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
        ! ! IF(jk==1)THEN
        !  write(*,*)'pressure',jk,jc,jb,press_hyd(jc,jk,jb),p_hyd(jk)
        ! ! ENDIF
        !  ENDIF
        ! END DO
      END DO
    END DO
  END SUBROUTINE calc_internal_press_new
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
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
  SUBROUTINE calc_internal_press(p_patch_3d, rho, prism_thick_c, h, press_hyd)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: p_patch_3d
    REAL(wp), INTENT(in)              :: rho          (1:nproma,1:n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)  !< density
    REAL(wp), INTENT(in), TARGET :: prism_thick_c(1:nproma,1:n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)              :: h            (1:nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)     !< surface elevation at cells
    REAL(wp), INTENT(inout)           :: press_hyd    (1:nproma,1:n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)  !< hydrostatic pressure
    
    ! local variables:
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !       & routine = (this_mod_name//':calc_internal_pressure')
    INTEGER :: slev, end_lev     ! vertical start and end level
    INTEGER :: jc, jk, jb
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    REAL(wp) :: z_full, z_box
    !   REAL(wp), POINTER :: del_zlev_m(:)
    REAL(wp),PARAMETER :: z_grav_rho_inv=rho_inv*grav
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    ! #slo# due to nag -nan compiler-option set intent(inout) variables to zero
    !press_hyd(:,:,:) = 0.0_wp
    all_cells => p_patch%cells%ALL
    press_hyd(:,:,:) = 0.0_wp
    
    slev = 1
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
      
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
        end_lev = p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
        
        IF(end_lev>=min_dolic)THEN
          DO jk = slev, end_lev
            IF(p_patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
              !             del_zlev_m => prism_thick_c(jc,:,jb)
              !             z_box      = del_zlev_m(jk)*rho(jc,jk,jb)      !-rho_ref!&!     pressure in single box at layer jk
              z_box = prism_thick_c(jc,jk,jb) * rho(jc,jk,jb)      !-rho_ref!&!     pressure in single box at layer jk
              
              press_hyd(jc,jk,jb) = ( z_full + 0.5_wp*z_box ) * z_grav_rho_inv
              ! rho_inv*grav  !hydrostatic press at level jk
              ! =half of pressure at actual box+ sum of all boxes above
              z_full              = z_full + z_box
              !           ELSE
              !             press_hyd(jc,jk,jb) = 0.0_wp
            ENDIF
          END DO
        ENDIF
      END DO
    END DO
    
  END SUBROUTINE calc_internal_press
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Calculates the density via a call to the equation-of-state.
  !! Several options for EOS are provided.
  !!
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2009)
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !!
  SUBROUTINE calc_density(p_patch_3d,tracer, rho)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: p_patch_3d
    REAL(wp),    INTENT(in), TARGET :: tracer(:,:,:,:)     !< input of S and T
    REAL(wp), INTENT(inout), TARGET :: rho   (:,:,:)       !< density
    
    ! local variables:
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !      & routine = (this_mod_name//':calc_density')
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3d%p_patch_2d(1)
    !---------------------------------------------------------------------
    ! CALL message (TRIM(routine), 'start')
    
    !For calc_density_lin_EOS and calc_density_MPIOM the conversion to in-situ temperature is done
    !internally.
    SELECT CASE (eos_type)
    CASE(1)
      CALL calc_density_lin_eos(p_patch_3d, tracer, rho)
    CASE(2)
      CALL calc_density_mpiom(p_patch_3d, tracer, rho)
    CASE(3)
      CALL calc_density_jmdwfg06_eos(p_patch_3d, tracer, rho)
      !CALL calc_density_JM_EOS(p_patch, tracer, rho)
    CASE default
      
    END SELECT
    
  END SUBROUTINE calc_density
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE calc_potential_density(p_patch_3d,tracer, rhopot)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: p_patch_3d
    REAL(wp),    INTENT(in), TARGET :: tracer(:,:,:,:)     !< input of S and T
    REAL(wp), INTENT(inout), TARGET :: rhopot(:,:,:)       !< density
    
    ! local variables:
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !      & routine = (this_mod_name//':calc_density')
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3d%p_patch_2d(1)
    !---------------------------------------------------------------------
    ! CALL message (TRIM(routine), 'start')
    
    !For calc_density_lin_EOS and calc_density_MPIOM the conversion to in-situ temperature is done
    !internally.
    !  SELECT CASE (EOS_TYPE)
    !    CASE(1)
    !      CALL calc_density_lin_EOS(p_patch_3D, tracer, rhopot)
    !    CASE(2)
    CALL calc_potential_density_mpiom(p_patch_3d, tracer, rhopot)
    !    CASE(3)
    !      CALL calc_density_JMDWFG06_EOS(p_patch_3D, tracer, rhopot)
    !      !CALL calc_density_JM_EOS(p_patch, tracer, rho)
    !    CASE DEFAULT
    
    !  END SELECT
    
  END SUBROUTINE calc_potential_density
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !!Calculates the density via a linear equation-of-state.
  !!
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2009)
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !!
  SUBROUTINE calc_density_lin_eos(p_patch_3d, tracer, rho)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: p_patch_3d
    REAL(wp),    INTENT(in)       :: tracer(:,:,:,:)     !< input of S and T
    REAL(wp), INTENT(inout)       :: rho   (:,:,:)       !< density
    
    ! local variables:
    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    all_cells => p_patch%cells%ALL
    
    IF(no_tracer==2)THEN
      
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
        !  tracer 1: potential temperature
        !  tracer 2: salinity
        DO jk=1, n_zlev
          DO jc = i_startidx, i_endidx
            IF(p_patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
              rho(jc,jk,jb) = rho_ref                      &
                & - a_t * tracer(jc,jk,jb,1)   &
                & + b_s * tracer(jc,jk,jb,2)
              !write(123,*)'density',jk,jc,jb,rho_ref, tracer(jc,jk,jb,1),&
              ! &tracer(jc,jk,jb,2),rho(jc,jk,jb), a_T, b_S
            ELSE
              ! rho(jc,jk,jb) = 0.0_wp
              rho(jc,jk,jb) = rho_ref   !  plotting purpose
            ENDIF
          END DO
        END DO
      END DO
      
    ELSEIF(no_tracer==1)THEN
      
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
        !  tracer 1: potential temperature
        DO jk=1, n_zlev
          DO jc = i_startidx, i_endidx
            IF(p_patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
              rho(jc,jk,jb) = rho_ref - a_t * tracer(jc,jk,jb,1) +b_s*sal_ref
              !write(123,*)'density',jk,jc,jb,rho(jc,jk,jb), tracer(jc,jk,jb,1),a_T
            ELSE
              rho(jc,jk,jb) = rho_ref   !  plotting purpose
            ENDIF
          END DO
        END DO
      END DO
      
    ENDIF
    
  END SUBROUTINE calc_density_lin_eos
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !!Calculates the density via a linear equation-of-state.
  !!
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2009)
  !!
  FUNCTION calc_density_lin_eos_func(t,s,p) result(rho)
    !
    REAL(wp),INTENT(in) :: t
    REAL(wp),INTENT(in) :: s
    REAL(wp),INTENT(in) :: p     !  pressure is unused
    REAL(wp)            :: rho   !< density
    
    rho = rho_ref - a_t * t  + b_s * s
    
  END FUNCTION calc_density_lin_eos_func
  !---------------------------------------------------------------------------
  
  !---------------------------------------------------------------------------
  !>
  ! !DESCRIPTION:
  !
  !  Calculation of the density as a function of salinity and temperature
  !  using the Jackett et al. (2006) equation of state. It uses exactly the
  !  same polynomial formulation as in McDougall et al. (2003), implemented
  !  in subroutine calculate_density_MDJWF03_EOS, but with a revised set of
  !  coefficients.
  !
  !  Check values are:
  !
  !    rho(theta=25 degC, S=35 PSU, p=2000 dbar) = 1031.65212332355 kg/m^3
  !    rho(theta=20 degC, S=20 PSU, p=1000 dbar) = 1017.84289041198 kg/m^3
  !
  !  Reference:
  !
  !    Jackett, D.R., T.J. McDougall, D.G. Wright, R. Feistel, and S.M. Griffies,
  !    2006: Algorithms for Density, Potential Temperature, Conservative
  !    Temperature, and the Freezing Temperature of Seawater. JAOT, 23, 1709-1728
  !
  !    McDougall, T.J., D.R. Jackett, D.G. Wright, and R. Feistel,  2003:
  !    Accurate and Computationally Efficient Algorithms for Potential
  !    Temperature and Density of Seawater. JAOT, 20, 730-741
  !
  ! !REVISION HISTORY:
  ! implemented by Peter Herrmann (2009)
  !
  SUBROUTINE calc_density_jmdwfg06_eos(p_patch_3d, tracer, rho)
    TYPE(t_patch_3d ),TARGET, INTENT(in):: p_patch_3d
    REAL(wp),    INTENT(in)                    :: tracer(:,:,:,:)
    REAL(wp), INTENT(inout)                    :: rho(:,:,:)       !< density
    
    ! !LOCAL VARIABLES:
    REAL(wp)::  z_p
    
    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------------------------------------
    !write(*,*)'inside EOS 06'
    all_cells => p_patch%cells%ALL
    
    !  tracer 1: potential temperature
    !  tracer 2: salinity
    IF(no_tracer==2)THEN
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
        DO jk=1, n_zlev
          DO jc = i_startidx, i_endidx
            ! operate on wet ocean points only
            IF(p_patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
              z_p=sfc_press_bar ! rho_ref*v_base%zlev_m(jk)*SItodBar
              rho(jc,jk,jb) = calc_density_jmdwfg06_eos_func(tracer(jc,jk,jb,1),&
                & tracer(jc,jk,jb,2),&
                & z_p)
              !           write(*,*)'rho',jc,jk,jb,rho(jc,jk,jb)
            END IF
          END DO
        END DO
      END DO
    ELSE IF(no_tracer==1)THEN
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
        DO jk=1, n_zlev
          DO jc = i_startidx, i_endidx
            ! operate on wet ocean points only
            IF(p_patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
              z_p=sfc_press_bar ! rho_ref*v_base%zlev_m(jk)*SItodBar
              rho(jc,jk,jb) = calc_density_jmdwfg06_eos_func(tracer(jc,jk,jb,1),&
                & sal_ref,      &
                & z_p)
              !           write(*,*)'rho',jc,jk,jb,rho(jc,jk,jb)
            END IF
          END DO
        END DO
      END DO
    ENDIF
    
  END SUBROUTINE calc_density_jmdwfg06_eos
  !----------------------------------------------------------------


  !----------------------------------------------------------------
  !>
  ! !DESCRIPTION:
  !
  !  Calculation of the density as a function of salinity and temperature
  !  using the Jackett et al. (2006) equation of state. It uses exactly the
  !  same polynomial formulation as in McDougall et al. (2003), implemented
  !  in subroutine calculate_density_MDJWF03_EOS, but with a revised set of
  !  coefficients.
  !
  !  Check values are:
  !
  !    rho(theta=25 degC, S=35 PSU, p=2000 dbar) = 1031.65212332355 kg/m^3 1031.6505605657569 1031.6505605657569
  !    rho(theta=20 degC, S=20 PSU, p=1000 dbar) = 1017.84289041198 kg/m^3
  !
  !  Reference:
  !
  !    Jackett, D.R., T.J. McDougall, D.G. Wright, R. Feistel, and S.M. Griffies,
  !    2006: Algorithms for Density, Potential Temperature, Conservative
  !    Temperature, and the Freezing Temperature of Seawater. JAOT, 23, 1709-1728
  !
  !    McDougall, T.J., D.R. Jackett, D.G. Wright, and R. Feistel,  2003:
  !    Accurate and Computationally Efficient Algorithms for Potential
  !    Temperature and Density of Seawater. JAOT, 20, 730-741
  !
  ! !REVISION HISTORY:
  ! implemented by Peter Herrmann (2009)
  !
  FUNCTION calc_density_jmdwfg06_eos_func(tracer_t, tracer_s, p) result(rho)
    REAL(wp), INTENT(in)       :: tracer_t
    REAL(wp), INTENT(in)       :: tracer_s
    REAL(wp), INTENT(in)       :: p
    REAL(wp)                   :: rho       !< density
    
    ! EOS variables, following the naming of the MITgcm implementation
    REAL (wp)  :: locpres, t1, t2, s1, p1, rhonum, sp5, p1t1, den, rhoden
    !-------------------------------------------------------------------------------------------------------
    !write(*,*)'inside EOS 06'
    
    ! tracer_t=25.0_wp
    ! tracer_s=35.0_wp
    ! p=2000.0_wp
    
    ! abbreviations
    t1 = tracer_t
    t2 = t1*t1
    s1 = tracer_s
    p1 = p
    
    rhonum = eosjmdwfgnum(0)                               &
      & + t1*(eosjmdwfgnum(1)                                  &
      & +     t1*(eosjmdwfgnum(2) + eosjmdwfgnum(3)*t1) )      &
      & + s1*(eosjmdwfgnum(4)                                  &
      & +     eosjmdwfgnum(5)*t1  + eosjmdwfgnum(6)*s1)        &
      & + p1*(eosjmdwfgnum(7) + eosjmdwfgnum(8)*t2             &
      & +     eosjmdwfgnum(9)*s1                               &
      & +     p1*(eosjmdwfgnum(10) + eosjmdwfgnum(11)*t2) )
    
    ! calculate the denominator of the Jackett et al.
    ! equation of state
    IF ( s1 .GT. 0.0_wp ) THEN
      sp5 = SQRT(s1)
    ELSE
      s1  = 0.0_wp
      sp5 = 0.0_wp
    END IF
    
    p1t1 = p1*t1
    den = eosjmdwfgden(0)                                          &
      & + t1*(eosjmdwfgden(1)                                     &
      & +     t1*(eosjmdwfgden(2)                                 &
      & +         t1*(eosjmdwfgden(3) + t1*eosjmdwfgden(4) ) ) )  &
      & + s1*(eosjmdwfgden(5)                                     &
      & +     t1*(eosjmdwfgden(6)                                 &
      & +         eosjmdwfgden(7)*t2)                             &
      & +     sp5*(eosjmdwfgden(8) + eosjmdwfgden(9)*t2) )        &
      & + p1*(eosjmdwfgden(10)                                    &
      & +     p1t1*(eosjmdwfgden(11)*t2 + eosjmdwfgden(12)*p1) )
    
    rhoden = 1.0_wp/(dbl_eps+den)
    
    !rhoLoc  = rhoNum*rhoDen - rho_ref
    rho     = rhonum*rhoden
    
    ! &rhoConst*9.80665_wp*dz*SItodBar,rhoConst*9.80616_wp*dz*SItodBar,dz,&
    ! &locPres*SItodBar
  END FUNCTION calc_density_jmdwfg06_eos_func
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!  Calculates density as a function of potential temperature and salinity
  !! using the Jackett and McDougall equation of state
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2011)
  !! Code below is an adaption of Sergey Danilov's implementation in
  !! the AWI Finite-Volume model.
  !!
  SUBROUTINE calc_density_jm_eos(p_patch_3d, tracer, rho)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: p_patch_3d
    REAL(wp),    INTENT(in)                     :: tracer(:,:,:,:)
    REAL(wp), INTENT(inout)                       :: rho(:,:,:)
    
    ! local variables:
    REAL(wp) :: z_t
    REAL(wp) :: z_s
    REAL(wp) :: z_rhopot, z_bulk, pz !,z_in_situ_temp
    !INTEGER  :: slev, end_lev
    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3d%p_patch_2d(1)
    !---------------------------------------------------------------------------
    all_cells => p_patch%cells%ALL
    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
      
      DO jk=1, n_zlev
        DO jc = i_startidx, i_endidx
          IF ( p_patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            
            pz  = p_patch_3d%p_patch_1d(1)%zlev_m(jk)
            
            z_t = tracer(jc,jk,jb,1)
            z_s = tracer(jc,jk,jb,2)
            
            !compute secant bulk modulus
            z_bulk = 19092.56_wp + z_t*(209.8925_wp             &
              & - z_t*(3.041638_wp - z_t*(-1.852732e-3_wp          &
              & - z_t*(1.361629e-5_wp))))                          &
              & + z_s*(104.4077_wp - z_t*(6.500517_wp              &
              & - z_t*(0.1553190_wp - z_t*(-2.326469e-4_wp))))      &
              & + SQRT(z_s**3)*(-5.587545_wp                       &
              & + z_t*(0.7390729_wp - z_t*(1.909078e-2_wp)))       &
              & - pz *(4.721788e-1_wp + z_t*(1.028859e-2_wp        &
              & + z_t*(-2.512549e-4_wp - z_t*(5.939910e-7_wp))))   &
              & - pz*z_s*(-1.571896e-2_wp                          &
              & - z_t*(2.598241e-4_wp + z_t*(-7.267926e-6_wp)))    &
              & - pz*SQRT(z_s**3)                                  &
              & *2.042967e-3_wp + pz*pz*(1.045941e-5_wp            &
              & - z_t*(5.782165e-10_wp - z_t*(1.296821e-7_wp)))    &
              & + pz*pz*z_s                                        &
              & *(-2.595994e-7_wp                                  &
              & + z_t*(-1.248266e-9_wp + z_t*(-3.508914e-9_wp)))
            
            z_rhopot = ( 999.842594_wp                     &
              & + z_t*( 6.793952e-2_wp                        &
              & + z_t*(-9.095290e-3_wp                        &
              & + z_t*( 1.001685e-4_wp                        &
              & + z_t*(-1.120083e-6_wp                        &
              & + z_t*( 6.536332e-9_wp)))))                   &
              & + z_s*( 0.824493_wp                           &
              & + z_t *(-4.08990e-3_wp                        &
              & + z_t *( 7.64380e-5_wp                        &
              & + z_t *(-8.24670e-7_wp                        &
              & + z_t *( 5.38750e-9_wp)))))                   &
              & + SQRT(z_s**3)*(-5.72466e-3_wp                &
              & + z_t*( 1.02270e-4_wp                         &
              & + z_t*(-1.65460e-6_wp)))                      &
              & + 4.8314e-4_wp*z_s**2)
            
            rho(jc,jk,jb) = z_rhopot/(1.0_wp + 0.1_wp*pz/z_bulk)&
              & - rho_ref
            ! write(*,*)'density ',jc,jk,jb,rho(jc,jk,jb)
            
          ENDIF
        END DO
      END DO
    END DO
    STOP
  END SUBROUTINE calc_density_jm_eos
  !----------------------------------------------------------------
  
  !----------------------------------------------------------------
  !
  ! !DESCRIPTION:
  !
  !!  Calculates density as a function of potential temperature and salinity
  !! using the equation of state as described in Gill, Atmosphere-Ocean Dynamics, Appendix 3
  !! The code below is copied from MPIOM
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2011)
  !!
  SUBROUTINE calc_density_mpiom(p_patch_3d, tracer, rho)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: p_patch_3d
    REAL(wp), INTENT(in)                   :: tracer(:,:,:,:)
    REAL(wp), INTENT(inout)                :: rho(:,:,:)       !< density
    
    ! !LOCAL VARIABLES:
    ! loop indices
    REAL(wp):: z_p(n_zlev)
    INTEGER :: jc, jk, jb
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    all_cells => p_patch%cells%ALL
    !i_len      = SIZE(dz)
    
    ! compute pressure first
    DO jk=1, n_zlev
      z_p(jk) = p_patch_3d%p_patch_1d(1)%zlev_m(jk)*rho_ref*sitodbar
    END DO
    !  tracer 1: potential temperature
    !  tracer 2: salinity
    IF(no_tracer==2)THEN
      
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
        DO jc = i_startidx, i_endidx
          DO jk=1, p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb) ! operate on wet ocean points only
            rho(jc,jk,jb) = calc_density_mpiom_func( tracer(jc,jk,jb,1), tracer(jc,jk,jb,2), z_p(jk))
          END DO
        END DO
      END DO
      
    ELSEIF(no_tracer==1)THEN
      
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
        DO jc = i_startidx, i_endidx
          DO jk=1, p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb) ! operate on wet ocean points only
            rho(jc,jk,jb) = calc_density_mpiom_func( tracer(jc,jk,jb,1), sal_ref, z_p(jk))
          END DO
        END DO
      END DO
      
    ENDIF
  END SUBROUTINE calc_density_mpiom
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Calculates density as a function of potential temperature and salinity
  !! using the equation of state as described in Gill, Atmosphere-Ocean Dynamics, Appendix 3
  !! The code below is copied from MPIOM. Note that within the sbr the potential temperature
  !! is converted into in-situ temperature !
  !! The code was checked with testvalues from Gill's book.
  !! For testing insert values here:
  !!   s = 35.0_wp,  t = 25.0_wp, p = 1000.0_wp, s3h = SQRT(s**3)
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2011)
  !!
  FUNCTION calc_density_mpiom_func(tpot, sal, p) result(rho)
    REAL(wp), INTENT(in) :: tpot, sal, p
    REAL(wp)             :: rho
    
    REAL(wp) :: dvs, fne, fst, qn3, qnq, qvs, s, s3h, t, denom
    REAL(wp), PARAMETER :: z_sref = 35.0_wp
    
    !This is the adisit part, that transforms potential in in-situ temperature
    qnq = -p * (-a_a3 + p * a_c3)
    qn3 = -p * a_a4
    qvs = (p * (a_b1 - a_d * p))*(sal - z_sref) + p * (a_a1 + p * (a_c1 - a_e1 * p))
    dvs = (a_b2 * p)*(sal - z_sref) + 1.0_wp + p * (-a_a2 + p * (a_c2 - a_e2 * p))
    t   = (tpot + qvs)/dvs
    fne = - qvs + t*(dvs + t*(qnq + t*qn3)) - tpot
    fst = dvs + t*(2._wp*qnq + 3._wp*qn3*t)
    t   = t - fne/fst
    s   = MAX(sal, 0.0_wp)
    s3h = SQRT(s**3)
    
    rho = r_a0 + t * (r_a1 + t * (r_a2 + t * (r_a3 + t * (r_a4 + t * r_a5))))&
      & + s * (r_b0 + t * (r_b1 + t * (r_b2 + t * (r_b3 + t * r_b4))))    &
      & + r_d0 * s**2                                                     &
      & + s3h * (r_c0 + t * (r_c1 + r_c2 * t))
    denom = 1._wp                                                            &
      & - p / (p * (r_h0 + t * (r_h1 + t * (r_h2 + t * r_h3))            &
      & + s * (r_ai0 + t * (r_ai1 + r_ai2 * t))              &
      & + r_aj0 * s3h                                        &
      & + (r_ak0 + t * (r_ak1 + t * r_ak2)                   &
      & + s * (r_am0 + t * (r_am1 + t * r_am2))) * p)        &
      & + r_e0 + t * (r_e1 + t * (r_e2 + t * (r_e3 + t * r_e4)))  &
      & + s * (r_f0 + t * (r_f1 + t * (r_f2 + t * r_f3)))         &
      & + s3h * (r_g0 + t * (r_g1 + r_g2 * t)))
    rho = rho/denom
    
  END FUNCTION calc_density_mpiom_func
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  ! !DESCRIPTION:
  !
  !!  Calculates density as a function of potential temperature and salinity
  !! using the equation of state as described in Gill, Atmosphere-Ocean Dynamics, Appendix 3
  !! The code below is copied from MPIOM
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2011)
  !!
  SUBROUTINE calc_potential_density_mpiom(p_patch_3d, tracer, rhopot)
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: p_patch_3d
    REAL(wp), INTENT(in)                   :: tracer(:,:,:,:)
    REAL(wp), INTENT(inout)                :: rhopot(:,:,:)
    
    ! !LOCAL VARIABLES:
    ! loop indices
    REAL(wp):: z_p(n_zlev)
    INTEGER :: jc, jk, jb
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3d%p_patch_2d(1)
    all_cells => p_patch%cells%ALL
    !-------------------------------------------------------------------------
    
    !  tracer 1: potential temperature
    !  tracer 2: salinity
    IF(no_tracer==2)THEN
      rhopot(:,:,:) = calc_potential_density_mpiom_elemental( tracer(:,:,:,1), tracer(:,:,:,2))
    ELSEIF(no_tracer==1)THEN
      rhopot(:,:,:) = calc_potential_density_mpiom_elemental( tracer(:,:,:,1), 0.0_wp*tracer(:,:,:,1)+sal_ref)
    ENDIF
  END SUBROUTINE calc_potential_density_mpiom
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  ! elemental version of the above computation
  ELEMENTAL FUNCTION calc_density_mpiom_elemental(tpot, sal, p) result(rho)
    REAL(wp), INTENT(in) :: tpot, sal, p
    REAL(wp)             :: rho
    
    REAL(wp) :: dvs, fne, fst, qn3, qnq, qvs, s, s3h, t, denom
    REAL(wp), PARAMETER :: z_sref = 35.0_wp
    
    
    !This is the adisit part, that transforms potential in in-situ temperature
    qnq = -p * (-a_a3 + p * a_c3)
    qn3 = -p * a_a4
    qvs = (p * (a_b1 - a_d * p))*(sal - z_sref) + p * (a_a1 + p * (a_c1 - a_e1 * p))
    dvs = (a_b2 * p)*(sal - z_sref) + 1.0_wp + p * (-a_a2 + p * (a_c2 - a_e2 * p))
    t   = (tpot + qvs)/dvs
    fne = - qvs + t*(dvs + t*(qnq + t*qn3)) - tpot
    fst = dvs + t*(2._wp*qnq + 3._wp*qn3*t)
    t   = t - fne/fst
    s   = MAX(sal, 0.0_wp)
    s3h = SQRT(s**3)
    
    rho = r_a0 + t * (r_a1 + t * (r_a2 + t * (r_a3 + t * (r_a4 + t * r_a5))))&
      & + s * (r_b0 + t * (r_b1 + t * (r_b2 + t * (r_b3 + t * r_b4))))    &
      & + r_d0 * s**2                                                     &
      & + s3h * (r_c0 + t * (r_c1 + r_c2 * t))
    denom = 1._wp                                                            &
      & - p / (p * (r_h0 + t * (r_h1 + t * (r_h2 + t * r_h3))            &
      & + s * (r_ai0 + t * (r_ai1 + r_ai2 * t))              &
      & + r_aj0 * s3h                                        &
      & + (r_ak0 + t * (r_ak1 + t * r_ak2)                   &
      & + s * (r_am0 + t * (r_am1 + t * r_am2))) * p)        &
      & + r_e0 + t * (r_e1 + t * (r_e2 + t * (r_e3 + t * r_e4)))  &
      & + s * (r_f0 + t * (r_f1 + t * (r_f2 + t * r_f3)))         &
      & + s3h * (r_g0 + t * (r_g1 + r_g2 * t)))
    rho = rho/denom
    
  END FUNCTION calc_density_mpiom_elemental
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  ! potential density wrt to surface
  ELEMENTAL FUNCTION calc_potential_density_mpiom_elemental(tpot, sal) result(rho)
    REAL(wp), INTENT(in) :: tpot, sal
    REAL(wp)             :: rho
    
    REAL(wp), PARAMETER :: p = 0.0_wp
    
    REAL(wp)             :: dvs, fne, fst, qn3, qnq, qvs, s, s3h, t, denom
    REAL(wp), PARAMETER :: z_sref = 35.0_wp
    
    !This is the adisit part, that transforms potential in in-situ temperature
    qnq = -p * (-a_a3 + p * a_c3)
    qn3 = -p * a_a4
    qvs = (p * (a_b1 - a_d * p))*(sal - z_sref) + p * (a_a1 + p * (a_c1 - a_e1 * p))
    dvs = (a_b2 * p)*(sal - z_sref) + 1.0_wp + p * (-a_a2 + p * (a_c2 - a_e2 * p))
    t   = (tpot + qvs)/dvs
    fne = - qvs + t*(dvs + t*(qnq + t*qn3)) - tpot
    fst = dvs + t*(2._wp*qnq + 3._wp*qn3*t)
    t   = t - fne/fst
    s   = MAX(sal, 0.0_wp)
    s3h = SQRT(s**3)
    
    rho = r_a0 + t * (r_a1 + t * (r_a2 + t * (r_a3 + t * (r_a4 + t * r_a5))))&
      & + s * (r_b0 + t * (r_b1 + t * (r_b2 + t * (r_b3 + t * r_b4))))    &
      & + r_d0 * s**2                                                     &
      & + s3h * (r_c0 + t * (r_c1 + r_c2 * t))
    denom = 1._wp                                                            &
      & - p / (p * (r_h0 + t * (r_h1 + t * (r_h2 + t * r_h3))            &
      & + s * (r_ai0 + t * (r_ai1 + r_ai2 * t))              &
      & + r_aj0 * s3h                                        &
      & + (r_ak0 + t * (r_ak1 + t * r_ak2)                   &
      & + s * (r_am0 + t * (r_am1 + t * r_am2))) * p)        &
      & + r_e0 + t * (r_e1 + t * (r_e2 + t * (r_e3 + t * r_e4)))  &
      & + s * (r_f0 + t * (r_f1 + t * (r_f2 + t * r_f3)))         &
      & + s3h * (r_g0 + t * (r_g1 + r_g2 * t)))
    rho = rho/denom
    
  END FUNCTION calc_potential_density_mpiom_elemental
  !-------------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------------
  !>
  ! !DESCRIPTION:
  !
  !!  Calculates potential tempertaure from in-situ temperature.
  !!  Formula described in Gill, Atmosphere-Ocean Dynamics, Appendix 3
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2011)
  !!
  SUBROUTINE convert_insitu2pot_temp(p_patch_3d, rho_ref, temp_insitu, sal, temp_pot)
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: p_patch_3d
    REAL(wp)                      :: rho_ref
    REAL(wp)                      :: temp_insitu(:,:,:)
    REAL(wp)                      :: sal(:,:,:)
    REAL(wp)                      :: temp_pot(:,:,:)
    
    ! !LOCAL VARIABLES:
    ! loop indices
    REAL(wp):: z_press
    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------------------------------------
    all_cells => p_patch%cells%ALL
    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
      DO jk=1, n_zlev
        z_press=sfc_press_bar ! rho_ref*grav*v_base%zlev_m(jk)
        DO jc = i_startidx, i_endidx
          ! operate on wet ocean points only
          IF(p_patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            
            temp_pot(jc,jk,jb) = convert_insitu2pot_temp_func(temp_insitu(jc,jk,jb),&
              & sal(jc,jk,jb),&
              & z_press)
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE convert_insitu2pot_temp
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !!  Calculates density as a function of potential temperature and salinity
  !! using the equation of state as described in Gill, Atmosphere-Ocean Dynamics, Appendix 3
  !! The code below is copied from MPIOM
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2011)
  !!
  FUNCTION convert_insitu2pot_temp_func(t, s, p) result(temp_pot)
    REAL(wp), INTENT(in) :: t, s, p
    REAL(wp)             :: temp_pot
    
    REAL(wp) :: z_s_ref
    !!---------------------------------------------------------------------------
    !z_s_ref= s_ref    !  s_ref is initial salinity, not reference
    z_s_ref= sal_ref  !  sal_ref = 35.0 = constant salinity reference
    temp_pot=t-p*(a_a1+ a_a2*t-a_a3*t*t+a_a4*t*t*t) &
      & -p*(s-z_s_ref)*(a_b1 -a_b2*t)           &
      & -p*p*(a_c1 -a_c2*t + a_c3*t*t)          &
      & +a_d*(s-z_s_ref)*p*p                    &
      & -p*p*p*(-a_e1 + a_e2*t)
    
  END FUNCTION convert_insitu2pot_temp_func
  !-------------------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------------------
  SUBROUTINE convert_pot_temp2insitu(p_patch_3d,trac_t, trac_s, temp_insitu)
    !
    ! !DESCRIPTION:
    !
    !!  Calculates potential tempertaure from in-situ temperature.
    !!  Formula described in Gill, Atmosphere-Ocean Dynamics, Appendix 3
    !! @par Revision History
    !! Initial version by Peter Korn, MPI-M (2011)
    !!
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: p_patch_3d
    REAL(wp)                                    :: trac_t(:,:,:)
    REAL(wp)                                    :: trac_s(:,:,:)
    REAL(wp)                                    :: temp_insitu(:,:,:)
    
    ! !LOCAL VARIABLES:
    ! loop indices
    REAL(wp):: z_press
    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------------------------------------
    all_cells => p_patch%cells%ALL
    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
      DO jk=1, n_zlev
        z_press=rho_ref*p_patch_3d%p_patch_1d(1)%zlev_m(jk)*sitodbar ! grav
        DO jc = i_startidx, i_endidx
          ! operate on wet ocean points only
          IF(p_patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            
            temp_insitu(jc,jk,jb) = adisit(trac_t(jc,jk,jb), trac_s(jc,jk,jb), z_press)
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE convert_pot_temp2insitu
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  FUNCTION adisit(th, sh, pa) result(temp_insitu)
    !
    !**** *ADISIT*  - TRANSFORMS POTENTIAL TO IN-SITU TEMPERATURE.
    !
    !     MODIFIED
    !     --------
    !     O. BOEHRINGER     *DKRZ*                   95
    !        - THIS VERSION USES ONLY 44 % OF THE CPU OF THE ORIGINAL HOPC VERSION
    !     UWE MIKOLAJEWICZ 2/99
    !     ==>ONE-DIMENSIONAL ARRAY, MERGE LOOPS
    !
    !     METHOD.
    !     --------
    !     TRANSFORMATION FROM POTENTIAL TO IN SITU TEMPERATURE
    !
    !**   INTERFACE.
    !     ----------
    !     *CALL* *ADISIT(TH,SH,PA)*       CALLED FROM *OCTHER*.
    !
    !     *COMMON*    *"PARAM1*            - OCEAN GRID DIMENSIONS.
    !
    !
    !     INPUT:
    !     -----
    !     *TH*        POTENTIAL TEMPERATURE [DEG C]
    !     *SH*        SALINITY  [PSU.]
    !     *PA*        PRESSURE  [PA]
    !
    !     OUTPUT:
    !     ------
    !     *TH*        IN-SITU  TEMPERATURE [DEG C]
    !
    ! ------------------------------------------------------------------------------
    !
    !
    REAL(wp), INTENT(in) :: pa, sh
    REAL(wp), INTENT(inout) :: th
    REAL(wp) temp_insitu
    REAL(wp) :: pr, dc, dv, dvs, fne, fst, qc, qn3, qnq, qv, qvs, t, tpo
    REAL(wp), PARAMETER :: z_sref=35.0_wp
    
    pr=pa
    !
    !  CHECK VALUES
    !     TH(1)=8.4678516
    !     SH(1)= 25.
    !     PR=1000.
    !
    qc = pr * (a_a1 + pr * (a_c1 - a_e1 * pr))
    qv = pr * (a_b1 - a_d * pr)
    dc = 1._wp + pr * (-a_a2 + pr * (a_c2 - a_e2 * pr))
    dv = a_b2 * pr
    qnq  = -pr * (-a_a3 + pr * a_c3)
    qn3  = -pr * a_a4
    !
    !DO i = 1, len
    !
    tpo = th
    qvs = qv*(sh - z_sref) + qc
    dvs = dv*(sh - z_sref) + dc
    t   = (tpo + qvs)/dvs
    fne = - qvs + t*(dvs + t*(qnq + t*qn3)) - tpo
    fst = dvs + t*(2._wp*qnq + 3._wp*qn3*t)
    !th = t - fne/fst
    temp_insitu=t - fne/fst
    !ENDDO
  END FUNCTION adisit
  !------------------------------------------------------------------------------


END MODULE mo_oce_thermodyn

