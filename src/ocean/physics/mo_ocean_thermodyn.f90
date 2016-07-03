
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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ocean_thermodyn
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_ocean_nml,           ONLY: n_zlev, eos_type, no_tracer, fast_performance_level,l_partial_cells, &
    & LinearThermoExpansionCoefficient, LinearHalineContractionCoefficient,OceanReferenceDensity
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_impl_constants,      ONLY: sea_boundary, sea_boundary, min_dolic !, &
  USE mo_exception,           ONLY: finish
  USE mo_loopindices,         ONLY: get_indices_c!, get_indices_e, get_indices_v
  USE mo_physical_constants,  ONLY: grav, sal_ref, rho_inv,  &
    & sitodbar, sfc_press_bar
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_parallel_config,     ONLY: nproma
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  IMPLICIT NONE
  
  PRIVATE
  
 
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_thermodyn'
! CHARACTER(len=12)           :: str_module    = 'oce_thermody'  ! Output of module for 1 line debug
  
  ! PUBLIC :: ocean_correct_ThermoExpansion
  PUBLIC :: calc_internal_press
  PUBLIC :: calculate_density,calc_potential_density
  PUBLIC :: calculate_density_onColumn
  PUBLIC :: calc_internal_press_grad

  PUBLIC :: calculate_density_mpiom_onColumn, calculate_density_jmdwfg06_onColumn
  !each specific EOS comes as a sbr and as a function. The sbr version is private as it is
  !only used in "calc_internal_press", whilethe function version is used in mo_ocean_physics
  !(sbr "update_ho_params") to calculate the local Richardson number.
  !PUBLIC :: density_linear_function
  !PUBLIC :: density_jmdwfg06_function
  !PUBLIC :: density_mpiom_function
  PUBLIC :: convert_insitu2pot_temp_func

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

  ! for the TEOS10
  ! copied from MOM5.1
  REAL(wp), PARAMETER :: mbfj_rho   = 1.017775176234136d+3
  REAL(wp), PARAMETER :: mbfj_alpha = 2.435473441547041d-4
  REAL(wp), PARAMETER :: mbfj_beta  = 7.284367916939847d-4
  REAL(wp), PARAMETER :: mb_neutralrho=1033.093610463980

  REAL(wp), PARAMETER :: v01 =  9.998420897506056d+2
  REAL(wp), PARAMETER :: v02 =  2.839940833161907
  REAL(wp), PARAMETER :: v03 = -3.147759265588511d-2
  REAL(wp), PARAMETER :: v04 =  1.181805545074306d-3
  REAL(wp), PARAMETER :: v05 = -6.698001071123802
  REAL(wp), PARAMETER :: v06 = -2.986498947203215d-2
  REAL(wp), PARAMETER :: v07 =  2.327859407479162d-4
  REAL(wp), PARAMETER :: v08 = -3.988822378968490d-2
  REAL(wp), PARAMETER :: v09 =  5.095422573880500d-4
  REAL(wp), PARAMETER :: v10 = -1.426984671633621d-5
  REAL(wp), PARAMETER :: v11 =  1.645039373682922d-7
  REAL(wp), PARAMETER :: v12 = -2.233269627352527d-2
  REAL(wp), PARAMETER :: v13 = -3.436090079851880d-4
  REAL(wp), PARAMETER :: v14 =  3.726050720345733d-6
  REAL(wp), PARAMETER :: v15 = -1.806789763745328d-4
  REAL(wp), PARAMETER :: v16 =  6.876837219536232d-7
  REAL(wp), PARAMETER :: v17 = -3.087032500374211d-7
  REAL(wp), PARAMETER :: v18 = -1.988366587925593d-8
  REAL(wp), PARAMETER :: v19 = -1.061519070296458d-11
  REAL(wp), PARAMETER :: v20 =  1.550932729220080d-10
  REAL(wp), PARAMETER :: v21 =  1.0
  REAL(wp), PARAMETER :: v22 =  2.775927747785646d-3
  REAL(wp), PARAMETER :: v23 = -2.349607444135925d-5
  REAL(wp), PARAMETER :: v24 =  1.119513357486743d-6
  REAL(wp), PARAMETER :: v25 =  6.743689325042773d-10
  REAL(wp), PARAMETER :: v26 = -7.521448093615448d-3
  REAL(wp), PARAMETER :: v27 = -2.764306979894411d-5
  REAL(wp), PARAMETER :: v28 =  1.262937315098546d-7
  REAL(wp), PARAMETER :: v29 =  9.527875081696435d-10
  REAL(wp), PARAMETER :: v30 = -1.811147201949891d-11
  REAL(wp), PARAMETER :: v31 = -3.303308871386421d-5
  REAL(wp), PARAMETER :: v32 =  3.801564588876298d-7
  REAL(wp), PARAMETER :: v33 = -7.672876869259043d-9
  REAL(wp), PARAMETER :: v34 = -4.634182341116144d-11
  REAL(wp), PARAMETER :: v35 =  2.681097235569143d-12
  REAL(wp), PARAMETER :: v36 =  5.419326551148740d-6
  REAL(wp), PARAMETER :: v37 = -2.742185394906099d-5
  REAL(wp), PARAMETER :: v38 = -3.212746477974189d-7
  REAL(wp), PARAMETER :: v39 =  3.191413910561627d-9
  REAL(wp), PARAMETER :: v40 = -1.931012931541776d-12
  REAL(wp), PARAMETER :: v41 = -1.105097577149576d-7
  REAL(wp), PARAMETER :: v42 =  6.211426728363857d-10
  REAL(wp), PARAMETER :: v43 = -1.119011592875110d-10
  REAL(wp), PARAMETER :: v44 = -1.941660213148725d-11
  REAL(wp), PARAMETER :: v45 = -1.864826425365600d-14
  REAL(wp), PARAMETER :: v46 =  1.119522344879478d-14
  REAL(wp), PARAMETER :: v47 = -1.200507748551599d-15
  REAL(wp), PARAMETER :: v48 =  6.057902487546866d-17

  
CONTAINS
    !-------------------------------------------------------------------------
    !>
    !! Calculation the hydrostatic pressure gradient at edges by computing the pressure of the
    !! two adjacent cell fluid column as weight of the fluid column above a certain level.
    !! In this routine this level is given by the edge-level, this level goes down from the surface
    !! to the the deepest edge. At the deepest edge-level the horizontal gradient is taken by using
    !! the adjacent pressure values at this level. It might happen that one of the adjacent column
    !! reaches down deeper, but this does not influence the hydrostatic pressure. 
    !! This routine does calculathe the pressure only temporarily. 
    !! It overcomes a difficulty with the pressure gradient calculation for partial cells that arises
    !! in the subroutine "calc_internal_press" (see below). The calc_internal_press should not be used
    !! with partial cells.
    !!
    !! @par Revision History
    !! Initial version by Peter Korn, MPI-M (2014)
    !!
  !<Optimize:inUse>
  SUBROUTINE calc_internal_press_grad(patch_3d, rho, pressure_hyd, grad_coeff, press_grad)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(in)                 :: rho          (nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)  !< density
    REAL(wp), INTENT(inout)              :: pressure_hyd (nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp), INTENT(in), TARGET      :: prism_thick_e(1:nproma,1:n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)                 :: grad_coeff(:,:,:)
    REAL(wp), INTENT(inout)              :: press_grad    (nproma,n_zlev, patch_3d%p_patch_2d(1)%nblks_e)  !< hydrostatic pressure gradient

    ! local variables:
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !       & routine = (this_mod_name//':calc_internal_pressure')
    INTEGER :: je, jk, jb, jc, ic1,ic2,ib1,ib2
    INTEGER :: i_startblk, i_endblk, start_index, end_index
    REAL(wp) :: z_full_c1, z_box_c1, z_full_c2, z_box_c2
    REAL(wp) :: z_grav_rho_inv
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
    REAL(wp), POINTER :: prism_thick_e(:,:,:)
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: prism_center_dist !distance between prism centers without surface elevation
    !-----------------------------------------------------------------------
    z_grav_rho_inv = OceanReferenceDensity_inv * grav
    patch_2D        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    all_cells       => patch_2D%cells%ALL
    prism_thick_e   => patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_e

    iidx => patch_3D%p_patch_2D(1)%edges%cell_idx
    iblk => patch_3D%p_patch_2D(1)%edges%cell_blk

    pressure_hyd (1:nproma,1:n_zlev, 1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    !-------------------------------------------------------------------------

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_index, end_index)

      DO jc = start_index, end_index

        pressure_hyd(jc,1,jb) = rho(jc,1,jb)*z_grav_rho_inv*patch_3D%p_patch_1d(1)%constantPrismCenters_Zdistance(jc,1,jb)

        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

          pressure_hyd(jc,jk,jb) = pressure_hyd(jc,jk-1,jb) + 0.5_wp*(rho(jc,jk,jb)+rho(jc,jk-1,jb))&
            &*z_grav_rho_inv*patch_3D%p_patch_1d(1)%constantPrismCenters_Zdistance(jc,jk,jb)

        END DO
      END DO
    END DO
!ICON_OMP_END_DO

!ICON_OMP_DO PRIVATE(start_index, end_index, je, ic1, ib1, ic2, ib2, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)

      DO je = start_index, end_index

          ic1=patch_2D%edges%cell_idx(je,jb,1)
          ib1=patch_2D%edges%cell_blk(je,jb,1)
          ic2=patch_2D%edges%cell_idx(je,jb,2)
          ib2=patch_2D%edges%cell_blk(je,jb,2)

        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,jb)

          press_grad(je,jk,jb)=(pressure_hyd(ic2,jk,ib2)-pressure_hyd(ic1,jk,ib1))*grad_coeff(je,jk,jb)
          ! write(1234,*)'pressure', je,jk,jb,press_grad(je,jk,jb),pressure_hyd(ic2,jk,ib2),pressure_hyd(ic1,jk,ib1)
        END DO
      END DO
    END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    
    END SUBROUTINE calc_internal_press_grad
    !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Calculation the hydrostatic pressure gradient at edges by computing the pressure of the
  !! two adjacent cell fluid column as weight of the fluid column above a certain level.
  !! In this routine this level is given by the edge-level, this level goes down from the surface
  !! to the the deepest edge. At the deepest edge-level the horizontal gradient is taken by using
  !! the adjacent pressure values at this level. It might happen that one of the adjacent column
  !! reaches down deeper, but this does not influence the hydrostatic pressure. 
  !! This routine does calculathe the pressure only temporarily. 
  !! It overcomes a difficulty with the pressure gradient calculation for partial cells that arises
  !! in the subroutine "calc_internal_press" (see below). The calc_internal_press should not be used
  !! with partial cells.
  !!
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2014)
  !!
!<Optimize:inUse>
  SUBROUTINE calc_internal_press_grad0(patch_3d, rho, grad_coeff, press_grad)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)              :: rho          (nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)  !< density
    !REAL(wp), INTENT(in), TARGET      :: prism_thick_e(1:nproma,1:n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)              :: grad_coeff(:,:,:)   
    REAL(wp), INTENT(inout)           :: press_grad    (nproma,n_zlev, patch_3d%p_patch_2d(1)%nblks_e)  !< hydrostatic pressure gradient    
    
    ! local variables:
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !       & routine = (this_mod_name//':calc_internal_pressure')
    INTEGER :: je, jk, jb, jc, ic1,ic2,ib1,ib2
    INTEGER :: i_startblk, i_endblk, start_index, end_index
    REAL(wp) :: z_full_c1, z_box_c1, z_full_c2, z_box_c2
    REAL(wp) :: z_grav_rho_inv
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
    REAL(wp), POINTER :: prism_thick_e(:,:,:),prism_thick_c(:,:,:)
    REAL(wp) :: press_c1
    REAL(wp) :: press_c2
    !-----------------------------------------------------------------------
    z_grav_rho_inv = OceanReferenceDensity_inv * grav
    patch_2D        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain   
    prism_thick_e   => patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_e
    prism_thick_c   => patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c 
    
    iidx => patch_3D%p_patch_2D(1)%edges%cell_idx
    iblk => patch_3D%p_patch_2D(1)%edges%cell_blk
    !-------------------------------------------------------------------------

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, je, ic1, ib1, ic2, ib2, jk, z_full_c1,  &
!ICON_OMP z_full_c2, z_box_c1, z_box_c2, press_c1, press_c2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      
      DO je = start_index, end_index 
     
        ic1=patch_2D%edges%cell_idx(je,jb,1)
        ib1=patch_2D%edges%cell_blk(je,jb,1)
        ic2=patch_2D%edges%cell_idx(je,jb,2)
        ib2=patch_2D%edges%cell_blk(je,jb,2)
        
        z_full_c1 = 0.0_wp
        z_full_c2 = 0.0_wp
        
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,jb)

          z_box_c1 = prism_thick_c(ic1,jk,ib1) * rho(ic1,jk,ib1)
          z_box_c2 = prism_thick_c(ic2,jk,ib2) * rho(ic2,jk,ib2)
          
          press_c1 = ( z_full_c1 + 0.5_wp * z_box_c1 ) * z_grav_rho_inv
          press_c2 = ( z_full_c2 + 0.5_wp * z_box_c2 ) * z_grav_rho_inv
          
          press_grad(je,jk,jb)=(press_c2-press_c1)*grad_coeff(je,jk,jb)
          
          z_full_c1 = z_full_c1 + z_box_c1
          z_full_c2 = z_full_c2 + z_box_c2         

        END DO
!        ENDIF
      END DO
    END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    
  END SUBROUTINE calc_internal_press_grad0
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Calculation the hydrostatic pressure by computing the weight of the
  !! fluid column above a certain level.
  !! IMPORTANT: Do not use this with partial cells !!
  !!
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2009)
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !! Modified by Stephan Lorenz,        MPI-M (2010-10-22)
  !!  - division by OceanReferenceDensity included
  !!
!<Optimize:inUse>
  SUBROUTINE calc_internal_press(patch_3d, rho, prism_thick_c, h, press_hyd)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)              :: rho          (1:nproma,1:n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)  !< density
    REAL(wp), INTENT(in), TARGET :: prism_thick_c(1:nproma,1:n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)              :: h            (1:nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)     !< surface elevation at cells
    REAL(wp), INTENT(inout)           :: press_hyd    (1:nproma,1:n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)  !< hydrostatic pressure
    
    ! local variables:
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !       & routine = (this_mod_name//':calc_internal_pressure')
    INTEGER :: jc, jk, jb
    INTEGER :: i_startblk, i_endblk, start_index, end_index
    REAL(wp) :: z_full, z_box
    !   REAL(wp), POINTER :: del_zlev_m(:)
    REAL(wp) :: z_grav_rho_inv
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    z_grav_rho_inv = OceanReferenceDensity_inv * grav
    patch_2D   => patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    IF(l_partial_cells)THEN
      CALL finish('mo_ocean_thermodyn: This pressure calculation does NOT work with partial cells!','!!')
    ENDIF
    
    !CALL message (TRIM(routine), 'start')
    ! #slo# due to nag -nan compiler-option set intent(inout) variables to zero
    !press_hyd(:,:,:) = 0.0_wp
    all_cells => patch_2D%cells%ALL
    ! press_hyd(:,:,:) = 0.0_wp
    
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk, z_full, z_box) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_index, end_index)
      
      DO jc = start_index, end_index
        !
        !  #slo# calculation of pressure due to elevation is done here
        !  including actual density of surface water
        !
        !z_full = grav * rho(jc,toplev,jb) * h(jc,jb)
        !  #slo# 2011-01-19 - elevation not considered:
        !   - in SWM ok, since density is constant
        !   - check to include h if tracers (T, S) are active
        z_full  = 0.0_wp
        
!        IF(end_lev>=min_dolic)THEN
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

          z_box = prism_thick_c(jc,jk,jb) * rho(jc,jk,jb)      !-OceanReferenceDensity!&!     pressure in single box at layer jk

          press_hyd(jc,jk,jb) = ( z_full + 0.5_wp * z_box ) * z_grav_rho_inv
          ! OceanReferenceDensity_inv*grav  !hydrostatic press at level jk
          ! =half of pressure at actual box+ sum of all boxes above
          z_full              = z_full + z_box
          
        END DO
!        ENDIF
      END DO
    END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    
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
!<Optimize:inUse>
  SUBROUTINE calculate_density(patch_3d,tracer, rho)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp),    INTENT(in), TARGET :: tracer(:,:,:,:)     !< input of S and T
    REAL(wp), INTENT(inout), TARGET :: rho   (:,:,:)       !< density
    
    ! local variables:
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !      & routine = (this_mod_name//':calculate_density')
    ! TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    ! patch_2D   => patch_3d%p_patch_2d(1)
    !---------------------------------------------------------------------
    ! CALL message (TRIM(routine), 'start')
    
    !For calculate_density_lin_EOS and calculate_density_MPIOM the conversion to in-situ temperature is done
    !internally.
    SELECT CASE (eos_type)
    CASE(1)
      CALL calculate_density_linear(patch_3d, tracer, rho)
    CASE(2)
      CALL calculate_density_mpiom(patch_3d, tracer, rho)
    CASE(3)
      CALL calculate_density_jmdwfg06(patch_3d, tracer, rho)
      !CALL calculate_density_JM_EOS(patch_2D, tracer, rho)
    CASE(10)
      CALL calculate_density_EOS10(patch_3d, tracer, rho)
    CASE default
      
    END SELECT
    
  END SUBROUTINE calculate_density
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! pressure is in dbars !
  FUNCTION calculate_density_onColumn(temperatute, salinity, p, levels) result(rho)
    INTEGER,  INTENT(in)       :: levels
    REAL(wp), INTENT(in)       :: temperatute(1:levels)
    REAL(wp), INTENT(in)       :: salinity(1:levels)
    REAL(wp), INTENT(in)       :: p(1:levels)
    REAL(wp)                   :: rho (1:levels)      !< density

    SELECT CASE (eos_type)
    CASE(1)
      rho(1:levels) = calculate_density_linear_onColumn( &
        & temperatute(1:levels),  salinity(1:levels), p(1:levels), levels)
    CASE(2)
      ! pressure for density_mpiom should be in bars 
      rho(1:levels) = calculate_density_mpiom_onColumn( &
        & temperatute(1:levels),  salinity(1:levels), p(1:levels)*0.1_wp, levels)
    CASE(3)
      rho(1:levels) = calculate_density_jmdwfg06_onColumn( &
        & temperatute(1:levels),  salinity(1:levels), p(1:levels), levels)
    CASE(10)
      rho(1:levels) = calculate_density_EOS10_onColumn( &
        & temperatute(1:levels),  salinity(1:levels), p(1:levels), levels)
    CASE default

    END SELECT

  END FUNCTION calculate_density_onColumn
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
!<Optimize:inUse>
  SUBROUTINE calc_potential_density(patch_3d,tracer, rhopot)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp),    INTENT(in), TARGET :: tracer(:,:,:,:)     !< input of S and T
    REAL(wp), INTENT(inout), TARGET :: rhopot(:,:,:)       !< density
    
    !For calculate_density_lin_EOS and calculate_density_MPIOM the conversion to in-situ temperature is done
    !internally.
    !  SELECT CASE (EOS_TYPE)
    !    CASE(1)
    !      CALL calculate_density_lin_EOS(patch_3d, tracer, rhopot)
    !    CASE(2)
    CALL calc_potential_density_mpiom(patch_3d, tracer, rhopot)
    !    CASE(3)
    !      CALL calculate_density_JMDWFG06_EOS(patch_3d, tracer, rhopot)
    !      !CALL calculate_density_JM_EOS(patch_2D, tracer, rho)
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
  SUBROUTINE calculate_density_linear(patch_3d, tracer, rho)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp),    INTENT(in)       :: tracer(:,:,:,:)     !< input of S and T
    REAL(wp), INTENT(inout)       :: rho   (:,:,:)       !< density
    
    INTEGER :: jc, jk, jb
    INTEGER :: start_index, end_index
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    all_cells => patch_2D%cells%ALL
    
    IF(no_tracer==2)THEN
      
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        !  tracer 1: potential temperature
        !  tracer 2: salinity
        rho(:,:,jb) = OceanReferenceDensity   !  plotting purpose
        DO jc = start_index, end_index
          DO jk=1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
!            IF(patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
              rho(jc,jk,jb) = OceanReferenceDensity          &
                & - LinearThermoExpansionCoefficient * tracer(jc,jk,jb,1)   &
                & + LinearHalineContractionCoefficient * tracer(jc,jk,jb,2)
              !write(123,*)'density',jk,jc,jb,OceanReferenceDensity, tracer(jc,jk,jb,1),&
              ! &tracer(jc,jk,jb,2),rho(jc,jk,jb), a_T, b_S
!            ELSE
!              ! rho(jc,jk,jb) = 0.0_wp
!              rho(jc,jk,jb) = OceanReferenceDensity   !  plotting purpose
!            ENDIF
          END DO
        END DO
      END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
      
    ELSEIF(no_tracer==1)THEN
      
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        rho(:,:,jb) = OceanReferenceDensity   !  plotting purpose
        DO jc = start_index, end_index
          DO jk=1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
!            IF(patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
              rho(jc,jk,jb) = OceanReferenceDensity - LinearThermoExpansionCoefficient * tracer(jc,jk,jb,1) &
              &+ LinearHalineContractionCoefficient * sal_ref
              !write(123,*)'density',jk,jc,jb,rho(jc,jk,jb), tracer(jc,jk,jb,1),a_T
!            ELSE
!              rho(jc,jk,jb) = OceanReferenceDensity   !  plotting purpose
!            ENDIF
          END DO
        END DO
      END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
      
    ENDIF
    
  END SUBROUTINE calculate_density_linear
  !-------------------------------------------------------------------------
 
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
  SUBROUTINE calculate_density_jmdwfg06(patch_3d, tracer, rho)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)                   :: tracer(:,:,:,:)
    REAL(wp), INTENT(inout)                :: rho(:,:,:)       !< density

    ! !LOCAL VARIABLES:
    ! loop indices
    REAL(wp):: z_p(n_zlev), salinityReference_column(n_zlev)
    INTEGER :: jc, jk, jb
    INTEGER :: levels
    INTEGER :: i_startblk, i_endblk, start_index, end_index
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    all_cells => patch_2D%cells%ALL

    !NOTE: here we use the Boussinesq approximation
    IF(no_tracer==2)THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, jc, levels, z_p) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        DO jc = start_index, end_index
           levels = patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
           z_p(1:levels) = patch_3d%p_patch_1d(1)%depth_CellMiddle(jc,1:levels,jb) * ReferencePressureIndbars
           rho(jc,1:levels,jb) = calculate_density_jmdwfg06_onColumn( &
             & tracer(jc,1:levels,jb,1),  tracer(jc,1:levels,jb,2), z_p(1:levels), levels)
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    ELSEIF(no_tracer==1)THEN
      salinityReference_column(1:n_zlev) = sal_ref

!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, jc, levels, z_p) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        DO jc = start_index, end_index
          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          z_p(1:levels) = patch_3d%p_patch_1d(1)%depth_CellMiddle(jc,1:levels,jb) * ReferencePressureIndbars
          rho(jc,1:levels,jb) = calculate_density_jmdwfg06_onColumn( &
             & tracer(jc,1:levels,jb,1),  salinityReference_column(1:levels), z_p(1:levels), levels)
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    ENDIF

    CALL dbg_print('calculate_density_jmdwfg06: rho', rho , "" ,5, patch_2D%cells%in_domain)

  END SUBROUTINE calculate_density_jmdwfg06
  !-------------------------------------------------------------------------

  SUBROUTINE calculate_density_EOS10(patch_3d, tracer, rho)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)                   :: tracer(:,:,:,:)
    REAL(wp), INTENT(inout)                :: rho(:,:,:)       !< density

    ! !LOCAL VARIABLES:
    ! loop indices
    REAL(wp):: z_p(n_zlev), salinityReference_column(n_zlev)
    INTEGER :: jc, jk, jb
    INTEGER :: levels
    INTEGER :: i_startblk, i_endblk, start_index, end_index
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    all_cells => patch_2D%cells%ALL

    !NOTE: here we use the Boussinesq approximation
    IF(no_tracer>=2)THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, jc, levels, z_p) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        DO jc = start_index, end_index
           levels = patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
           z_p(1:levels) = patch_3d%p_patch_1d(1)%depth_CellMiddle(jc,1:levels,jb) * ReferencePressureIndbars
           rho(jc,1:levels,jb) = calculate_density_EOS10_onColumn( &
             & tracer(jc,1:levels,jb,1),  tracer(jc,1:levels,jb,2), z_p(1:levels), levels)
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    ELSE
      CALL finish("calculate_density_jmdwfg06", "no_tracer < 2")
    ENDIF

    CALL dbg_print('calculate_density_EOS10: rho', rho , "" ,5, patch_2D%cells%in_domain)

  END SUBROUTINE calculate_density_EOS10
  !-------------------------------------------------------------------------

  !----------------------------------------------------------------
  !>
  !!  Calculates density as a function of potential temperature and salinity
  !! using the equation of state as described in Gill, Atmosphere-Ocean Dynamics, Appendix 3
  !! The code below is copied from MPIOM
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2011)
  !!
!<Optimize:inUse>
  SUBROUTINE calculate_density_mpiom(patch_3d, tracer, rho)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)                   :: tracer(:,:,:,:)
    REAL(wp), INTENT(inout)                :: rho(:,:,:)       !< density

    ! !LOCAL VARIABLES:
    ! loop indices
    REAL(wp):: z_p(n_zlev), salinityReference_column(n_zlev)
    INTEGER :: jc, jk, jb
    INTEGER :: levels
    INTEGER :: i_startblk, i_endblk, start_index, end_index
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3d%p_patch_2d(1)
    all_cells => patch_2D%cells%ALL
    salinityReference_column(1:n_zlev) = sal_ref
    !-------------------------------------------------------------------------

    !  tracer 1: potential temperature
    !  tracer 2: salinity
    IF (no_tracer == 2) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, jc, levels, z_p) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index) 
        DO jc = start_index, end_index
          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          z_p(1:levels) = patch_3d%p_patch_1d(1)%depth_CellMiddle(jc,1:levels,jb) * OceanReferenceDensity * sitodbar
          rho(jc,1:levels,jb) = calculate_density_mpiom_onColumn( &
            & tracer(jc,1:levels,jb,1),  tracer(jc,1:levels,jb,2), z_p(1:levels), levels)
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    ELSE
      
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, jc, levels, z_p) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        DO jc = start_index, end_index
          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          z_p(1:levels) = patch_3d%p_patch_1d(1)%depth_CellMiddle(jc,1:levels,jb) * OceanReferenceDensity * sitodbar
          rho(jc,1:levels,jb) = calculate_density_mpiom_onColumn( &
             & tracer(jc,1:levels,jb,1),  salinityReference_column(1:levels), z_p(1:levels), levels)
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    ENDIF ! no_tracer==2
    

    CALL dbg_print('calculate_density_mpiom: rho', rho , "" ,5, patch_2D%cells%in_domain)

  END SUBROUTINE calculate_density_mpiom
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!Calculates the density via a linear equation-of-state.
  !!
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2009)
  !!
  FUNCTION density_linear_function(t,s,p) result(rho)
    !
    REAL(wp),INTENT(in) :: t
    REAL(wp),INTENT(in) :: s
    REAL(wp),INTENT(in) :: p     !  pressure is unused
    REAL(wp)            :: rho   !< density

    rho = OceanReferenceDensity - LinearThermoExpansionCoefficient * t  &
    &+ LinearHalineContractionCoefficient * s

  END FUNCTION density_linear_function
  !---------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!Calculates the density via a linear equation-of-state.
  !!
  !! @par Revision History
  !! Initial version by Peter Korn, MPI-M (2009)
  !!
  FUNCTION calculate_density_linear_onColumn(t,s,p,levels) result(rho)
    !
    INTEGER,INTENT(in)  :: levels
    REAL(wp),INTENT(in) :: t(1:levels)
    REAL(wp),INTENT(in) :: s(1:levels)
    REAL(wp),INTENT(in) :: p(1:levels)     !  pressure is unused
    REAL(wp)            :: rho(1:levels)   !< density

    rho(1:levels) = OceanReferenceDensity - LinearThermoExpansionCoefficient * t(1:levels) &
    & + LinearHalineContractionCoefficient * s(1:levels)

  END FUNCTION calculate_density_linear_onColumn
  !---------------------------------------------------------------------------

  !----------------------------------------------------------------
  !>
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
  FUNCTION calculate_density_jmdwfg06_onColumn(temperatute, salinity, p, levels) result(rho)
    INTEGER,  INTENT(in)       :: levels
    REAL(wp), INTENT(in)       :: temperatute(1:levels)
    REAL(wp), INTENT(in)       :: salinity(1:levels)
    REAL(wp), INTENT(in)       :: p(1:levels)
    REAL(wp)                   :: rho (1:levels)      !< density

    ! EOS variables, following the naming of the MITgcm implementation
    REAL (wp)  :: t1(1:levels), t2(1:levels), s1(1:levels), p1(1:levels), &
      & rhonum(1:levels), sp5(1:levels), den(1:levels)! , rhoden(1:levels)
    !-------------------------------------------------------------------------------------------------------
    !write(*,*)'inside EOS 06'

    ! temperatute=25.0_wp
    ! salinity=35.0_wp
    ! p=2000.0_wp

    ! abbreviations
    s1(1:levels) = MAX(salinity(1:levels),0.0_wp)
    t1(1:levels) = temperatute(1:levels)
    t2(1:levels) = t1(1:levels) * t1(1:levels)
    p1(1:levels) = p(1:levels)
    sp5(1:levels) = SQRT(s1(1:levels))

    rhonum(1:levels) = eosjmdwfgnum(0)                         &
      & + t1*(eosjmdwfgnum(1)                                  &
      & +     t1*(eosjmdwfgnum(2) + eosjmdwfgnum(3)*t1) )      &
      & + s1*(eosjmdwfgnum(4)                                  &
      & +     eosjmdwfgnum(5)*t1  + eosjmdwfgnum(6)*s1)        &
      & + p1*(eosjmdwfgnum(7) + eosjmdwfgnum(8)*t2             &
      & +     eosjmdwfgnum(9)*s1                               &
      & +     p1*(eosjmdwfgnum(10) + eosjmdwfgnum(11)*t2) )

    ! calculate the denominator of the Jackett et al.
    ! equation of state
    den(1:levels) = eosjmdwfgden(0)                                                       &
      & + t1(1:levels) * (eosjmdwfgden(1)                                                 &
      & +     t1(1:levels) * (eosjmdwfgden(2)                                             &
      & +         t1(1:levels) * (eosjmdwfgden(3) + t1(1:levels) * eosjmdwfgden(4) ) ) )  &
      & + s1(1:levels) * (eosjmdwfgden(5)                                                 &
      & +     t1(1:levels) * (eosjmdwfgden(6)                                             &
      & +         eosjmdwfgden(7) * t2(1:levels))                                         &
      & +     sp5(1:levels) * (eosjmdwfgden(8) + eosjmdwfgden(9) * t2(1:levels)) )        &
      & + p1(1:levels) * (eosjmdwfgden(10)                                                &
      & +     p1(1:levels) * t1(1:levels) * (eosjmdwfgden(11) * t2(1:levels) + eosjmdwfgden(12)*p1(1:levels)) )

    ! rhoden = 1.0_wp / (dbl_eps+den)

    !rhoLoc  = rhoNum*rhoDen - OceanReferenceDensity
    ! rho     = rhonum * rhoden
    rho(1:levels)     = rhonum(1:levels) / den(1:levels)

    ! &rhoConst*9.80665_wp*dz*SItodBar,rhoConst*9.80616_wp*dz*SItodBar,dz,&
    ! &locPres*SItodBar
  END FUNCTION calculate_density_jmdwfg06_onColumn
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! EOS10 implementation from MOM5.1
  FUNCTION calculate_density_EOS10_onColumn(temperatute, salinity, pressure, levels) result(rho)
    INTEGER,  INTENT(in)       :: levels
    REAL(wp), INTENT(in)       :: temperatute(1:levels)
    REAL(wp), INTENT(in)       :: salinity(1:levels)
    REAL(wp), INTENT(in)       :: pressure(1:levels)
    REAL(wp)                   :: rho (1:levels)      !< density

    REAL (wp)  :: t1(1:levels), t2(1:levels), s1(1:levels), p1(1:levels), &
      & num(1:levels), sp5(1:levels), den(1:levels), p1t1(1:levels)! , rhoden(1:levels)
    ! standard atmospheric pressure (dbar).
    ! Should be set to 0.0 if assume zero pressure for the
    ! overlying atmosphere.  But if have a realistic atmospheric
    ! pressure loading, then set press_standard=10.1325.
    real :: press_standard(1:levels)
    !-------------------------------------------------------------------------------------------------------
    press_standard(1:levels) = 0.0
    t1(1:levels)  = temperatute(1:levels)
    s1(1:levels)  = salinity(1:levels)
    p1(1:levels)  = pressure(1:levels) - press_standard(1:levels)

    sp5(1:levels) = sqrt(s1)
    t2(1:levels)  = t1*t1
    p1t1(1:levels) = p1*t1

    num(1:levels) = v01 + t1*(v02 + t1*(v03 + v04*t1))                &
          + s1*(v05 + t1*(v06 + v07*t1)                     &
          + sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1))))       &
          + p1*(v12 + t1*(v13 + v14*t1) + s1*(v15 + v16*t1) &
          + p1*(v17 + t1*(v18 + v19*t1) + v20*s1))

    den(1:levels) =  v21 + t1*(v22 + t1*(v23 + t1*(v24 + v25*t1)))                &
          + s1*(v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + v36*s1 &
          + sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))       &
          + p1*(v37 + t1*(v38 + t1*(v39 + v40*t1))                     &
          + s1*(v41 + v42*t1)                                          &
          + p1*(v43 + t1*(v44 + v45*t1 + v46*s1)                       &
          + p1*(v47 + v48*t1)))


!              rho(1:levels) = num(1:levels)/(epsln+den(1:levels))
    rho(1:levels) = num(1:levels)/den(1:levels)

  END FUNCTION calculate_density_EOS10_onColumn
  !-------------------------------------------------------------------------

  !----------------------------------------------------------------
  !>
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
  FUNCTION density_jmdwfg06_function(temperatute, salinity, p) result(rho)
    REAL(wp), INTENT(in)       :: temperatute
    REAL(wp), INTENT(in)       :: salinity
    REAL(wp), INTENT(in)       :: p
    REAL(wp)                   :: rho       !< density

    ! EOS variables, following the naming of the MITgcm implementation
    REAL (wp)  :: t1, t2, s1, p1, rhonum, sp5, p1t1, den
    !-------------------------------------------------------------------------------------------------------
    !write(*,*)'inside EOS 06'

    ! temperatute=25.0_wp
    ! salinity=35.0_wp
    ! p=2000.0_wp

    ! abbreviations
    t1 = temperatute
    t2 = t1*t1
    s1 = salinity
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

    ! rhoden = 1.0_wp / (dbl_eps+den)

    !rhoLoc  = rhoNum*rhoDen - OceanReferenceDensity
    ! rho     = rhonum * rhoden
    rho     = rhonum / den

    ! &rhoConst*9.80665_wp*dz*SItodBar,rhoConst*9.80616_wp*dz*SItodBar,dz,&
    ! &locPres*SItodBar
  END FUNCTION density_jmdwfg06_function
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
  FUNCTION density_mpiom_function(tpot, sal, p) result(rho)
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
    
  END FUNCTION density_mpiom_function
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
!<Optimize:inUse>
  FUNCTION calculate_density_mpiom_onColumn(tpot, sal, p, levels) result(rho)
    INTEGER, INTENT(in)  :: levels
    REAL(wp), INTENT(in) :: tpot(1:levels), sal(1:levels), p(1:levels)
    REAL(wp)             :: rho(1:levels)

    REAL(wp) :: dvs(1:levels), fne(1:levels), fst(1:levels), qn3(1:levels), &
      & qnq(1:levels), qvs(1:levels), s(1:levels), s3h(1:levels),           &
      & t(1:levels), denom(1:levels), s__2(1:levels)

    REAL(wp), PARAMETER :: z_sref = 35.0_wp

    !This is the adisit part, that transforms potential in in-situ temperature
    qnq(1:levels) = -p(1:levels) * (-a_a3 + p(1:levels) * a_c3)
    qn3(1:levels) = -p(1:levels) * a_a4
    qvs(1:levels) = (p(1:levels) * (a_b1 - a_d * p(1:levels))) * &
      & (sal(1:levels) - z_sref) + p(1:levels) * (a_a1 + p(1:levels) * (a_c1 - a_e1 * p(1:levels)))

    dvs(1:levels) = (a_b2 * p(1:levels)) * (sal(1:levels) - z_sref) + &
      & 1.0_wp + p(1:levels) * (-a_a2 + p(1:levels) * (a_c2 - a_e2 * p(1:levels)))

    t(1:levels)   = (tpot(1:levels) + qvs(1:levels)) / dvs(1:levels)
    fne(1:levels) = - qvs(1:levels) + t(1:levels) *    &
      & (dvs(1:levels) + t(1:levels) * (qnq(1:levels) + t(1:levels) * qn3(1:levels))) &
      & - tpot(1:levels)

    fst(1:levels) = dvs(1:levels) + t(1:levels) * (2._wp * qnq(1:levels) + &
      & 3._wp * qn3(1:levels) * t(1:levels))

    t(1:levels)   = t(1:levels) - fne(1:levels) / fst(1:levels)
    s(1:levels)    = MAX(sal(1:levels), 0.0_wp)
    s__2(1:levels) = s(1:levels)**2
!     s3h(1:levels)  = SQRT(s__2(1:levels) * s(1:levels))
    s3h(1:levels)  = s(1:levels) * SQRT(s(1:levels))

    rho(1:levels) = r_a0 + t(1:levels) * &
      & (r_a1 + t(1:levels) * (r_a2 + t(1:levels) *                   &
      &   (r_a3 + t(1:levels) * (r_a4 + t(1:levels) * r_a5))))        &
      & + s(1:levels) * (r_b0 + t(1:levels) * (r_b1 + t(1:levels) &
      &     * (r_b2 + t(1:levels) * (r_b3 + t(1:levels) * r_b4))))    &
      & + r_d0 * s__2(1:levels)                                           &
      & + s3h(1:levels) * (r_c0 + t(1:levels) * (r_c1 + r_c2 * t(1:levels)))

    denom(1:levels) = 1._wp                                                           &
      & - p(1:levels) / (p(1:levels) * (r_h0 + t(1:levels) *                  &
      &     (r_h1 + t(1:levels) * (r_h2 + t(1:levels) * r_h3))                    &
      & + s(1:levels) * (r_ai0 + t(1:levels) * (r_ai1 + r_ai2 * t(1:levels))) &
      & + r_aj0 * s3h(1:levels)                                                       &
      & + (r_ak0 + t(1:levels) * (r_ak1 + t(1:levels) * r_ak2)                    &
      & + s(1:levels) * (r_am0 + t(1:levels) *                                    &
      &     (r_am1 + t(1:levels) * r_am2))) * p(1:levels))                        &
      & + r_e0 + t(1:levels) * (r_e1 + t(1:levels) *                              &
      &      (r_e2 + t(1:levels) * (r_e3 + t(1:levels) * r_e4)))                  &
      & + s(1:levels) * (r_f0 + t(1:levels) * (r_f1 + t(1:levels) *           &
      &       (r_f2 + t(1:levels) * r_f3)))                                           &
      & + s3h(1:levels) * (r_g0 + t(1:levels) * (r_g1 + r_g2 * t(1:levels))))

    rho(1:levels) = rho(1:levels)/denom(1:levels)

  END FUNCTION calculate_density_mpiom_onColumn
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
!<Optimize:inUse>
  SUBROUTINE calc_potential_density_mpiom(patch_3d, tracer, rhopot)
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)                   :: tracer(:,:,:,:)
    REAL(wp), INTENT(inout)                :: rhopot(:,:,:)
    
    ! !LOCAL VARIABLES:
    ! loop indices
    REAL(wp):: z_p(n_zlev)
    INTEGER :: jc, jk, jb

    INTEGER :: start_index, end_index
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    !-------------------------------------------------------------------------
    
    !  tracer 1: potential temperature
    !  tracer 2: salinity
    IF(no_tracer==2)THEN
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        DO jc = start_index, end_index
          DO jk=1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb) ! operate on wet ocean points only
            rhopot(jc,jk,jb) = calc_potential_density_mpiom_elemental( &
              & tracer(jc,jk,jb,1), tracer(jc,jk,jb,2))
          END DO
        END DO
      END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    ELSEIF(no_tracer==1)THEN
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        DO jc = start_index, end_index
          DO jk=1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb) ! operate on wet ocean points only
            rhopot(jc,jk,jb) = calc_potential_density_mpiom_elemental( &
              & tracer(jc,jk,jb,1), sal_ref)
          END DO
        END DO
      END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    ENDIF
  END SUBROUTINE calc_potential_density_mpiom
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! potential density wrt to surface
!<Optimize:inUse>
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
  SUBROUTINE convert_insitu2pot_temp(patch_3d, OceanReferenceDensity, temp_insitu, sal, temp_pot)
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp)                      :: OceanReferenceDensity
    REAL(wp)                      :: temp_insitu(:,:,:)
    REAL(wp)                      :: sal(:,:,:)
    REAL(wp)                      :: temp_pot(:,:,:)
    
    ! !LOCAL VARIABLES:
    ! loop indices
    REAL(wp):: z_press
    INTEGER :: jc, jk, jb
    INTEGER :: start_index, end_index
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL

    ! z_press is used uninitialized
    ! z_press=OceanReferenceDensity*patch_3d%p_patch_1d(1)%zlev_m(jk)*sitodbar ! grav
    
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_index, end_index)
      DO jc = start_index, end_index
        DO jk=1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb) ! operate on wet ocean points only
          temp_pot(jc,jk,jb) = convert_insitu2pot_temp_func(&
            & temp_insitu(jc,jk,jb),&
            & sal(jc,jk,jb),   &
            & z_press)
        END DO
      END DO
    END DO
 !ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
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
  
 
!  !-------------------------------------------------------------------------
!  !>
!  !!  Calculates density as a function of potential temperature and salinity
!  !! using the Jackett and McDougall equation of state
!  !! @par Revision History
!  !! Initial version by Peter Korn, MPI-M (2011)
!  !! Code below is an adaption of Sergey Danilov's implementation in
!  !! the AWI Finite-Volume model.
!  !!
!  SUBROUTINE calculate_density_jm_eos(patch_3d, tracer, rho)
!    !
!    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
!    REAL(wp),    INTENT(in)                     :: tracer(:,:,:,:)
!    REAL(wp), INTENT(inout)                       :: rho(:,:,:)
!
!    ! local variables:
!    REAL(wp) :: z_t
!    REAL(wp) :: z_s
!    REAL(wp) :: z_rhopot, z_bulk, pz !,z_in_situ_temp
!    !INTEGER  :: slev, end_lev
!    INTEGER :: jc, jk, jb
!    INTEGER :: start_index, end_index
!    TYPE(t_subset_range), POINTER :: all_cells
!    TYPE(t_patch), POINTER :: patch_2D
!    !-----------------------------------------------------------------------
!    patch_2D   => patch_3d%p_patch_2d(1)
!    !---------------------------------------------------------------------------
!    all_cells => patch_2D%cells%ALL
!
!    DO jb = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, jb, start_index, end_index)
!      DO jc = start_index, end_index
!        DO jk=1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
!          IF ( patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
!
!          pz  = patch_3d%p_patch_1d(1)%zlev_m(jk)
!
!          z_t = tracer(jc,jk,jb,1)
!          z_s = tracer(jc,jk,jb,2)
!
!          !compute secant bulk modulus
!          z_bulk = 19092.56_wp + z_t*(209.8925_wp             &
!            & - z_t*(3.041638_wp - z_t*(-1.852732e-3_wp          &
!            & - z_t*(1.361629e-5_wp))))                          &
!            & + z_s*(104.4077_wp - z_t*(6.500517_wp              &
!            & - z_t*(0.1553190_wp - z_t*(-2.326469e-4_wp))))      &
!            & + SQRT(z_s**3)*(-5.587545_wp                       &
!            & + z_t*(0.7390729_wp - z_t*(1.909078e-2_wp)))       &
!            & - pz *(4.721788e-1_wp + z_t*(1.028859e-2_wp        &
!            & + z_t*(-2.512549e-4_wp - z_t*(5.939910e-7_wp))))   &
!            & - pz*z_s*(-1.571896e-2_wp                          &
!            & - z_t*(2.598241e-4_wp + z_t*(-7.267926e-6_wp)))    &
!            & - pz*SQRT(z_s**3)                                  &
!            & *2.042967e-3_wp + pz*pz*(1.045941e-5_wp            &
!            & - z_t*(5.782165e-10_wp - z_t*(1.296821e-7_wp)))    &
!            & + pz*pz*z_s                                        &
!            & *(-2.595994e-7_wp                                  &
!            & + z_t*(-1.248266e-9_wp + z_t*(-3.508914e-9_wp)))
!
!          z_rhopot = ( 999.842594_wp                     &
!            & + z_t*( 6.793952e-2_wp                        &
!            & + z_t*(-9.095290e-3_wp                        &
!            & + z_t*( 1.001685e-4_wp                        &
!            & + z_t*(-1.120083e-6_wp                        &
!            & + z_t*( 6.536332e-9_wp)))))                   &
!            & + z_s*( 0.824493_wp                           &
!            & + z_t *(-4.08990e-3_wp                        &
!            & + z_t *( 7.64380e-5_wp                        &
!            & + z_t *(-8.24670e-7_wp                        &
!            & + z_t *( 5.38750e-9_wp)))))                   &
!            & + SQRT(z_s**3)*(-5.72466e-3_wp                &
!            & + z_t*( 1.02270e-4_wp                         &
!            & + z_t*(-1.65460e-6_wp)))                      &
!            & + 4.8314e-4_wp*z_s**2)
!
!          rho(jc,jk,jb) = z_rhopot/(1.0_wp + 0.1_wp*pz/z_bulk)&
!            & - OceanReferenceDensity
!          ! write(*,*)'density ',jc,jk,jb,rho(jc,jk,jb)
!
!          ! ENDIF
!        END DO
!      END DO
!    END DO
!
!    STOP
!  END SUBROUTINE calculate_density_jm_eos
!  !----------------------------------------------------------------
  
  !-------------------------------------------------------------------------------------
!  SUBROUTINE convert_pot_temp2insitu(patch_3d,trac_t, trac_s, temp_insitu)
!    !
!    ! !DESCRIPTION:
!    !
!    !!  Calculates potential tempertaure from in-situ temperature.
!    !!  Formula described in Gill, Atmosphere-Ocean Dynamics, Appendix 3
!    !! @par Revision History
!    !! Initial version by Peter Korn, MPI-M (2011)
!    !!
!
!    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
!    REAL(wp)                                    :: trac_t(:,:,:)
!    REAL(wp)                                    :: trac_s(:,:,:)
!    REAL(wp)                                    :: temp_insitu(:,:,:)
!
!    ! !LOCAL VARIABLES:
!    ! loop indices
!    REAL(wp):: z_press
!    INTEGER :: jc, jk, jb
!    INTEGER :: start_index, end_index
!    TYPE(t_subset_range), POINTER :: all_cells
!    TYPE(t_patch), POINTER :: patch_2D
!    !-----------------------------------------------------------------------
!    patch_2D   => patch_3d%p_patch_2d(1)
!    !-------------------------------------------------------------------------------------------------------
!    all_cells => patch_2D%cells%ALL
!
!    DO jb = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, jb, start_index, end_index)
!      DO jk=1, n_zlev
!        z_press=OceanReferenceDensity*patch_3d%p_patch_1d(1)%zlev_m(jk)*sitodbar ! grav
!        DO jc = start_index, end_index
!          ! operate on wet ocean points only
!          IF(patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
!
!            temp_insitu(jc,jk,jb) = adisit(trac_t(jc,jk,jb), trac_s(jc,jk,jb), z_press)
!          END IF
!        END DO
!      END DO
!    END DO
!  END SUBROUTINE convert_pot_temp2insitu
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
!  FUNCTION adisit(th, sh, pa) result(temp_insitu)
!    !
!    !**** *ADISIT*  - TRANSFORMS POTENTIAL TO IN-SITU TEMPERATURE.
!    !
!    !     MODIFIED
!    !     --------
!    !     O. BOEHRINGER     *DKRZ*                   95
!    !        - THIS VERSION USES ONLY 44 % OF THE CPU OF THE ORIGINAL HOPC VERSION
!    !     UWE MIKOLAJEWICZ 2/99
!    !     ==>ONE-DIMENSIONAL ARRAY, MERGE LOOPS
!    !
!    !     METHOD.
!    !     --------
!    !     TRANSFORMATION FROM POTENTIAL TO IN SITU TEMPERATURE
!    !
!    !**   INTERFACE.
!    !     ----------
!    !     *CALL* *ADISIT(TH,SH,PA)*       CALLED FROM *OCTHER*.
!    !
!    !     *COMMON*    *"PARAM1*            - OCEAN GRID DIMENSIONS.
!    !
!    !
!    !     INPUT:
!    !     -----
!    !     *TH*        POTENTIAL TEMPERATURE [DEG C]
!    !     *SH*        SALINITY  [PSU.]
!    !     *PA*        PRESSURE  [PA]
!    !
!    !     OUTPUT:
!    !     ------
!    !     *TH*        IN-SITU  TEMPERATURE [DEG C]
!    !
!    ! ------------------------------------------------------------------------------
!    !
!    !
!    REAL(wp), INTENT(in) :: pa, sh
!    REAL(wp), INTENT(inout) :: th
!    REAL(wp) temp_insitu
!    REAL(wp) :: pr, dc, dv, dvs, fne, fst, qc, qn3, qnq, qv, qvs, t, tpo
!    REAL(wp), PARAMETER :: z_sref=35.0_wp
!
!    pr=pa
!    !
!    !  CHECK VALUES
!    !     TH(1)=8.4678516
!    !     SH(1)= 25.
!    !     PR=1000.
!    !
!    qc = pr * (a_a1 + pr * (a_c1 - a_e1 * pr))
!    qv = pr * (a_b1 - a_d * pr)
!    dc = 1._wp + pr * (-a_a2 + pr * (a_c2 - a_e2 * pr))
!    dv = a_b2 * pr
!    qnq  = -pr * (-a_a3 + pr * a_c3)
!    qn3  = -pr * a_a4
!    !
!    !DO i = 1, len
!    !
!    tpo = th
!    qvs = qv*(sh - z_sref) + qc
!    dvs = dv*(sh - z_sref) + dc
!    t   = (tpo + qvs)/dvs
!    fne = - qvs + t*(dvs + t*(qnq + t*qn3)) - tpo
!    fst = dvs + t*(2._wp*qnq + 3._wp*qn3*t)
!    !th = t - fne/fst
!    temp_insitu=t - fne/fst
!    !ENDDO
!  END FUNCTION adisit
  !------------------------------------------------------------------------------
  
    
END MODULE mo_ocean_thermodyn

