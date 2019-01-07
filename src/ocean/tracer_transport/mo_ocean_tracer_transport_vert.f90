!>
!! Contains the implementation of the vertical tracer transport routines for the ICON ocean model.
!! This comprises vertical advection.
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/01)
!!
!! @par Copyright and License
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!----------------------------
#include "omp_definitions.inc"
#include "icon_definitions.inc"
!----------------------------
MODULE mo_ocean_tracer_transport_vert
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  !USE mo_math_utilities,            ONLY: t_cartesian_coordinates
  USE mo_impl_constants,            ONLY: sea_boundary, max_char_length
  USE mo_math_constants,            ONLY: dbl_eps
  USE mo_ocean_nml,                 ONLY: n_zlev,                    &
    & upwind, central,fct_vert_adpo,fct_vert_ppm,fct_vert_zalesak,  &
    & fct_vert_minmod, l_adpo_flowstrength,                         &
    & flux_calculation_vert, fast_performance_level
  USE mo_parallel_config,           ONLY: nproma
  USE mo_run_config,                ONLY: dtime, ltimer
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_adv_vert, timer_ppm_slim, &
    & timer_adpo_vert, timers_level  !, timer_dif_vert,
  USE mo_ocean_types,                 ONLY: t_verticalAdvection_ppm_coefficients, &
    & t_operator_coeff
  USE mo_model_domain,              ONLY: t_patch,t_patch_3d, t_patch_vert
  USE mo_exception,                 ONLY: finish !, message_text, message
  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_ocean_physics
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_c, sync_patch_array
  USE mo_ocean_limiter,             ONLY: v_ppm_slimiter_mo_onblock 
  USE mo_ocean_tracer_transport_types,  ONLY: t_ocean_transport_state

IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=12)           :: str_module    = 'oceTracVert '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  PUBLIC :: advect_flux_vertical
 
CONTAINS
  
  !-------------------------------------------------------------------------
  !! SUBROUTINE advects vertically the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
  !! mpi parallelized, sync required: trac_out
!<Optimize:inUse>
  SUBROUTINE advect_flux_vertical( patch_3d,   &
    & trac_old,             &
    & transport_state,          &
    & operators_coeff,      &
    & flux_div_vert)
    
    TYPE(t_patch_3d ),TARGET :: patch_3d
    REAL(wp), INTENT(inout)           :: trac_old(:,:,:) ! (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_ocean_transport_state), TARGET :: transport_state
    TYPE(t_operator_coeff), TARGET    :: operators_coeff
    REAL(wp), INTENT(inout)           :: flux_div_vert(:,:,:) ! (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !new tracer
    
    !Local variables
    REAL(wp) :: deriv_fst, deriv_sec, adpo_weight_upw_cntr, adpo_weight_cntr_upw
    REAL(wp) :: prism_volume, transport_in, transport_out, adpo_r1, adpo_r2
    INTEGER :: startIndex, endIndex
    INTEGER :: jc, jk, jb
    ! vertical advective tracer fluxes:
    REAL(wp), ALLOCATABLE :: z_adv_flux_v2(:,:,:) ! resulting flux
    REAL(wp), ALLOCATABLE :: z_adv_flux_v (:,:,:)  ! resulting flux
    REAL(wp), ALLOCATABLE :: z_adv_flux_vu(:,:,:)  ! upwind flux
    
    REAL(wp), ALLOCATABLE :: adpo_weight(:,:,:)
    REAL(wp) :: z_flux_div_upw, z_flux_div_cnt
    REAL(wp), ALLOCATABLE :: a_v(:,:,:)
    
    TYPE(t_patch), POINTER :: patch_2D
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_tracer_advection:advect_flux_vertical')
    !-------------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-------------------------------------------------------------------------------
    start_timer(timer_adv_vert,2)

    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    
    ! CALL sync_patch_array(sync_c, patch_2D, trac_old)
    
    ! This is already synced in  edges_in_domain !
    ! CALL sync_patch_array(SYNC_C, patch_2D, transport_state%w)
    

    IF (flux_calculation_vert == fct_vert_ppm) THEN

      ! Vertical advection scheme: piecewise parabolic method (ppm) inUse
      CALL upwind_vflux_ppm( patch_3d,              &
        & trac_old,                                 &
        & transport_state%w,              &
        & dtime, 1,                                &
        & patch_3d%p_patch_1d(1)%prism_thick_c,     &
        & patch_3d%p_patch_1d(1)%inv_prism_thick_c, &
        & operators_coeff%verticalAdvectionPPMcoeffs, &
        & flux_div_vert)

        stop_timer(timer_adv_vert,2)
        RETURN

    ENDIF

    
  END SUBROUTINE advect_flux_vertical
  !-------------------------------------------------------------------------
  
    
  !------------------------------------------------------------------------
  !! Otpimized version of the third order PPM scheme
  !!
!<Optimize:inUse>
  SUBROUTINE upwind_vflux_ppm(       &
    & patch_3d, tracer,                   &
    & w, dtime, vertical_limiter_type,    &
    & cell_thickeness,  cell_invheight,   &
    & verticalAdvection_ppm_coefficients, &
    & flux_div_vert)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d   
    REAL(wp), INTENT(inout)           :: tracer(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)  !< in: advected cell centered variable
    REAL(wp), INTENT(inout)           :: w(nproma,n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !<  in : vertical velocity
    REAL(wp), INTENT(in)              :: dtime  !< time step
    REAL(wp), INTENT(inout)           :: cell_thickeness(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< in: layer thickness at cell center at time n
    REAL(wp), INTENT(inout)           :: cell_invheight(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< layer thickness at cell center at time n
    TYPE(t_verticalAdvection_ppm_coefficients) :: verticalAdvection_ppm_coefficients(patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)           :: flux_div_vert(:,:,:) ! (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !new tracer
    INTEGER, INTENT(in)               :: vertical_limiter_type                                  !< parameter to select limiter
    !
    !-----------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    INTEGER                       :: startIndex, endIndex, jb
    !-----------------------------------------------------------------------
    cells_in_domain => patch_3d%p_patch_2d(1)%cells%in_domain
#ifdef NAGFOR
    flux_div_vert(:,:,:) = 0.0_wp
#endif
   
!ICON_OMP_PARALLEL_DO PRIVATE(startIndex, endIndex) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
      CALL  upwind_vflux_ppm_onBlock(         &
            & tracer(:,:,jb),                 &
            & w(:,:,jb),                      &
            & dtime, vertical_limiter_type,   &
            & cell_thickeness(:,:,jb),        &
            & cell_invheight(:,:,jb),         &
            & verticalAdvection_ppm_coefficients(jb), &
            & flux_div_vert(:,:,jb),          & ! out
            & startIndex, endIndex,           &
            & patch_3d%p_patch_1d(1)%dolic_c(:,jb))
    ENDDO ! end loop over blocks
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE upwind_vflux_ppm
  !-------------------------------------------------------------------------

  !------------------------------------------------------------------------
  !! Otpimized version of the third order PPM scheme
  !!
  !! Calculation of time averaged vertical tracer fluxes using the third
  !! order PPM scheme.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2009-08-12)
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized to height based vertical coordinate systems. Included
  !!   parameter coeff_grid, in order to apply the same code to either a
  !!   pressure based or height based vertical coordinate system.
  !! Modification by Daniel Reinert, DWD (2011-01-17)
  !! - added optional parameter opt_lout_edge which will provide the
  !!   reconstructed 'edge' value.
  !!
  !! mpi parallelized, sync required: p_upflux
  !
  ! !LITERATURE
  ! - Colella and Woodward (1984), JCP, 54, 174-201
  ! - Carpenter et al. (1989), MWR, 118, 586-612
  ! - Lin and Rood (1996), MWR, 124, 2046-2070
  !
  ! Optimized version of upwind_vflux_ppm
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE upwind_vflux_ppm_onBlock(  &
    & tracer,                            &
    & w, dtime, vertical_limiter_type,   &
    & cell_thickeness, cell_invheight,   &
    & ppmCoeffs,                         &
    & flux_div_vert,                     &
    & startIndex, endIndex, cells_noOfLevels)

    REAL(wp), INTENT(in)           :: tracer(nproma,n_zlev)      !< advected cell centered variable
    REAL(wp), INTENT(in)           :: w(nproma,n_zlev+1)         !<  in : vertical velocity
    REAL(wp), INTENT(in)           :: dtime                      !< time step
    REAL(wp), INTENT(in)           :: cell_thickeness(nproma,n_zlev) !< layer thickness at cell center at time n
    REAL(wp), INTENT(in)           :: cell_invheight(nproma,n_zlev)!< layer thickness at cell center at time n
    TYPE(t_verticalAdvection_ppm_coefficients) :: ppmCoeffs
    INTEGER, INTENT(in)            :: vertical_limiter_type                    !< parameter to select limiter
    REAL(wp), INTENT(inout)        :: flux_div_vert(nproma, n_zlev) !new tracer
    INTEGER, INTENT(in)            :: startIndex, endIndex
    INTEGER, INTENT(in)            :: cells_noOfLevels(nproma)
    !
    REAL(wp) :: upward_tracer_flux(nproma,n_zlev+1)      !< tracer flux
    REAL(wp) :: z_face(nproma,n_zlev+1)   !< face values of transported field
    REAL(wp) :: z_face_up(nproma,n_zlev)  !< face value (upper face)
    REAL(wp) :: z_face_low(nproma,n_zlev) !< face value (lower face)
    REAL(wp) :: z_lext_1(n_zlev+1)                 !< linear extrapolation value 1
    REAL(wp) :: z_lext_2(n_zlev+1)                 !< linear extrapolation value 2
    REAL(wp) :: z_cfl_m, z_cfl_p !< CFL number (weta>0, w<0), CFL number (weta<0, w>0)
    REAL(wp) :: z_slope(nproma,n_zlev+1)  !< monotonized slope
    REAL(wp) :: z_slope_u, z_slope_l                            !< one-sided slopes
    !< for weta >0 and weta <0
    REAL(wp) :: z_a11, z_a12                                    !< 1/6 * a6,i (see Colella and Woodward (1984))
    REAL(wp) :: z_weta_dt                                       !< weta times dtime
    INTEGER :: firstLevel, secondLevel                                    !< vertical start thisLevel and start thisLevel +1
    INTEGER :: levelAbove, levelBelow,level2Below                        !< vertical thisLevel minus and plus one, plus two
    INTEGER :: jc, thisLevel, cell_levels
    !LOGICAL  :: opt_lout_edge !< optional: output edge value (.TRUE.),
    !                          !< or the flux across the edge   !< (.FALSE./not specified)
    !REAL(wp) :: opt_topflx_tra(nproma,patch_3D%p_patch_2D(1)%alloc_cell_blocks)  !< vertical tracer flux at upper boundary
    INTEGER, PARAMETER :: islopel_vsm = 1
    
    REAL(wp), POINTER ::  cellHeightRatio_This_toBelow(:,:)
    REAL(wp), POINTER ::  cellHeightRatio_This_toThisBelow(:,:)
    REAL(wp), POINTER ::  cellHeight_2xBelow_x_RatioThis_toThisBelow(:,:)
    REAL(wp), POINTER ::  cellHeightRatio_This_toThisAboveBelow(:,:)
    REAL(wp), POINTER ::  cellHeightRatio_2xAboveplusThis_toThisBelow(:,:)
    REAL(wp), POINTER ::  cellHeightRatio_2xBelowplusThis_toThisAbove(:,:)
    REAL(wp), POINTER ::  cellHeightRatio_ThisAbove_to2xThisplusBelow(:,:)
    REAL(wp), POINTER ::  cellHeightRatio_ThisBelow_to2xThisplusAbove(:,:)
    REAL(wp), POINTER ::  cellHeight_inv_ThisAboveBelow2Below(:,:)
    !-----------------------------------------------------------------------
    cellHeightRatio_This_toBelow                   =>  ppmCoeffs%cellHeightRatio_This_toBelow
    cellHeightRatio_This_toThisBelow               =>  ppmCoeffs%cellHeightRatio_This_toThisBelow
    cellHeight_2xBelow_x_RatioThis_toThisBelow     =>  ppmCoeffs%cellHeight_2xBelow_x_RatioThis_toThisBelow
    cellHeightRatio_This_toThisAboveBelow          =>  ppmCoeffs%cellHeightRatio_This_toThisAboveBelow
    cellHeightRatio_2xAboveplusThis_toThisBelow    =>  ppmCoeffs%cellHeightRatio_2xAboveplusThis_toThisBelow
    cellHeightRatio_2xBelowplusThis_toThisAbove    =>  ppmCoeffs%cellHeightRatio_2xBelowplusThis_toThisAbove
    cellHeightRatio_ThisAbove_to2xThisplusBelow    =>  ppmCoeffs%cellHeightRatio_ThisAbove_to2xThisplusBelow
    cellHeightRatio_ThisBelow_to2xThisplusAbove    =>  ppmCoeffs%cellHeightRatio_ThisBelow_to2xThisplusAbove
    cellHeight_inv_ThisAboveBelow2Below            =>  ppmCoeffs%cellHeight_inv_ThisAboveBelow2Below

    firstLevel  = 1
    secondLevel = 2

    ! advection is done with an upwind scheme and a piecwise parabolic
    ! approx. of the subgrid distribution is used.
    ! 3 options:  standard without limiter
    !             standard with semi-monotone or monotone limiter
    !             special version with limiter which handles CFL >1
    !
    !------------------------------------------------
    ! 1. Calculate Courant number for weta>0 (w<0) and weta<0 (w>0)
    ! ..............

    ! 2. Calculate monotonized slope
    !
    z_slope(:, :) = 0._wp

! !CDIR NODEP
    DO jc = startIndex, endIndex

! !CDIR NODEP
      DO thisLevel = secondLevel, cells_noOfLevels(jc) - 1

        levelAbove    = thisLevel - 1                 ! index of top half thisLevel
        levelBelow    = thisLevel + 1  ! index of bottom half thisLevel

        z_slope_u = 2._wp * (tracer(jc,thisLevel)  - tracer(jc,levelAbove))
        z_slope_l = 2._wp * (tracer(jc,levelBelow) - tracer(jc,thisLevel))

        IF ((z_slope_u * z_slope_l) > 0._wp) THEN

          z_slope(jc,thisLevel) = &
            & ( cellHeightRatio_This_toThisAboveBelow(jc,thisLevel)  )           &
            & * ( &
            & (cellHeightRatio_2xAboveplusThis_toThisBelow(jc,thisLevel))        &
            & &
            & * (tracer(jc,levelBelow) - tracer(jc,thisLevel))                   &
            & &
            & + (cellHeightRatio_2xBelowplusThis_toThisAbove(jc,thisLevel))      &
            & &
            & * (tracer(jc,thisLevel) - tracer(jc,levelAbove)) )

          z_slope(jc,thisLevel) = SIGN(                                           &
            & MIN( ABS(z_slope(jc,thisLevel)), ABS(z_slope_u), ABS(z_slope_l) ),  &
            & z_slope(jc,thisLevel))

        ENDIF

      END DO ! jc = startIndex, endIndex
    END DO ! end loop over vertical levels

    !
    ! 3. reconstruct face values at vertical half-levels
    !
    ! Boundary values for two highest and lowest half-levels
    !
    ! for faces k=secondLevel and k=nlevp1-1 reconstructed face values are calculated by
    ! interpolating a quadratic (instead of quartic) polynomial through 3
    ! values of the indefinite integral A=\int_{\eta_{0}}^{\eta}q\,\mathrm{d}\eta
    !
    ! for faces k=firstLevel and k=nlevp1 a zero gradient condition is assumed and the
    ! face values are set to the tracer values of the corresponding cell centers
    !
    z_face(:, :)  = 0.0_wp

! !CDIR NODEP
    DO jc = startIndex, endIndex
      cell_levels = cells_noOfLevels(jc)

! !CDIR NODEP
      DO thisLevel = secondLevel, cell_levels - 2
        levelAbove  = thisLevel - 1
        levelBelow  = thisLevel + 1
        level2Below = thisLevel + 2

        z_face(jc,levelBelow) = tracer(jc,thisLevel) &
          & + (cellHeightRatio_This_toThisBelow(jc,thisLevel))                   &
          & &
          & * (tracer(jc,levelBelow) - tracer(jc,thisLevel))                     &
          & &
          & + cellHeight_inv_ThisAboveBelow2Below(jc,thisLevel)                  &
          & * &
          & ( &
          &   (cellHeight_2xBelow_x_RatioThis_toThisBelow(jc,thisLevel))         &
          & &
          &  * (cellHeightRatio_ThisAbove_to2xThisplusBelow(jc,thisLevel)        &
          &     -  cellHeightRatio_ThisBelow_to2xThisplusAbove(jc,thisLevel))    &
          & &
          &  * (tracer(jc,levelBelow) - tracer(jc,thisLevel))                    &
          &  - z_slope(jc,levelBelow)                                            &
          &  * cell_thickeness(jc,thisLevel)                                     &
          &  * cellHeightRatio_ThisAbove_to2xThisplusBelow(jc,thisLevel)         &
          & &
          &  +  z_slope(jc,thisLevel) *                                          &
          & &
          &  cell_thickeness(jc,levelBelow) * &
          &  cellHeightRatio_ThisBelow_to2xThisplusAbove(jc,levelBelow)          &
          & )

      END DO ! end loop over vertical levels

      ! compute top 2 levels
      IF ( cells_noOfLevels(jc) >= firstLevel  ) THEN
      
        z_face(jc,firstLevel) = tracer(jc,firstLevel)
        
        IF ( cell_levels >= secondLevel ) THEN

          z_face(jc,secondLevel) = &
            & tracer(jc,firstLevel) *                                        &
            & (1._wp - cellHeightRatio_This_toBelow(jc,firstLevel))          &
            & &
            & + (cellHeightRatio_This_toThisBelow(jc,firstLevel))            &
            & * &
            & ( cellHeightRatio_This_toBelow(jc,firstLevel)                  &
            & * tracer(jc,firstLevel)                                        & 
            & + tracer(jc,secondLevel))
            
        ENDIF
      ENDIF

      ! compute bottom thisLevel
      IF ( cells_noOfLevels(jc) > secondLevel ) THEN

        z_face(jc, cell_levels) =                                                       &
          & tracer(jc,cell_levels-1) *                                                  &
          & ( 1._wp - cellHeightRatio_This_toBelow(jc,cell_levels-1))                   &
          & +  &
          & (cell_thickeness(jc,cell_levels-1) / (cell_thickeness(jc,cell_levels-1)     &
          &   + cell_thickeness(jc,cell_levels)))                                       &
          &  * &
          &  (cellHeightRatio_This_toBelow(jc,cell_levels-1)                            &
          &  * tracer(jc,cell_levels-1)                                                 &
          &  + tracer(jc,cell_levels))

      ENDIF

    ENDDO

    ! 4. Limitation of first guess parabola (which is based on z_face)
    ! Note that z_face_up(k) does not need to equal z_face_low(k-1) after
    ! the limitation procedure.
    ! Therefore 2 additional fields z_face_up and z_face_low are
    ! introduced.
    z_face_low(1:nproma,1:n_zlev) = 0.0_wp
    z_face_up (1:nproma,1:n_zlev) = 0.0_wp

    IF (vertical_limiter_type == islopel_vsm) THEN
      !     ! monotonic (mo) limiter
      CALL v_ppm_slimiter_mo_onBlock( &
        & tracer(:,:),       &
        & z_face(:,:),       &
        & z_slope(:,:),      &
        & z_face_up,         &
        & z_face_low,        &
        & startIndex,        &
        & endIndex,          &
        & cells_noOfLevels)

    ELSE
        ! simply copy face values to 'face_up' and 'face_low' arrays
        DO jc = startIndex, endIndex
          DO thisLevel = secondLevel, cells_noOfLevels(jc)-1
            z_face_up(jc, thisLevel)  = z_face(jc, thisLevel    )
            z_face_low(jc,thisLevel)  = z_face(jc, thisLevel + 1)
          ENDDO
        END DO

    ENDIF  !  p_ityp_vlimit

    upward_tracer_flux(:,:) = 0.0_wp

! !CDIR NODEP
    DO jc = startIndex, endIndex
! !CDIR NODEP
      DO thisLevel = secondLevel, cells_noOfLevels(jc)
        ! index of top half thisLevel
        levelAbove = thisLevel - 1
        ! linear extrapolated values
        ! for the height based coordinate system multiplication by coeff_grid
        ! is not necessary due to compensating (-) signs.
        ! first (of cell above) (case of w < 0; weta > 0)
        ! z_delta_m = (z_face_low(jc,levelAbove) - z_face_up(jc,levelAbove))
        z_a11     = tracer(jc,levelAbove)                                  &
          & - 0.5_wp * (z_face_low(jc,levelAbove) + z_face_up(jc,levelAbove))

        ! Calculate local Courant number at half levels
        ! z_cfl_m for weta >0 (w <0)
        ! z_cfl_p for weta <0 (w >0)
        ! z_weta_dt = 0.0_wp
        z_weta_dt = ABS(w(jc,thisLevel)) * dtime
        z_cfl_p = z_weta_dt * cell_invheight(jc, thisLevel)
        z_cfl_m = z_weta_dt * cell_invheight(jc, levelAbove)

        z_lext_1(thisLevel) = tracer(jc,levelAbove)                             &
          & + 0.5_wp * (z_face_low(jc,levelAbove) - z_face_up(jc,levelAbove))   &
          &          * ( 1.0_wp - z_cfl_m )                                     &
          & - z_a11 - z_a11 * z_cfl_m * ( -3._wp  + 2._wp * z_cfl_m)

        ! second (of cell below) (case of w > 0; weta < 0)
        ! z_delta_p = (z_face_low(jc,thisLevel) - z_face_up(jc,thisLevel))
        z_a12     = tracer(jc,thisLevel)                                      &
          & - 0.5_wp * (z_face_low(jc,thisLevel) + z_face_up(jc,thisLevel))

        z_lext_2(thisLevel) = tracer(jc,thisLevel)                            &
          & - 0.5_wp * (z_face_low(jc,thisLevel) - z_face_up(jc,thisLevel))   &
          &          * ( 1.0_wp - z_cfl_p )                                   &
          & - z_a12 + z_a12 * z_cfl_p * (- 3._wp + 2._wp * z_cfl_p)
        !
        ! calculate vertical tracer flux
        !
        !upward_tracer_flux(jc,thisLevel) =                                    &
        !  & laxfr_upflux_v( w(jc,thisLevel),                                  &
        !  & z_lext_1(jc,thisLevel), z_lext_2(jc,thisLevel))
        ! copy of   FUNCTION laxfr_upflux_v( p_vn, p_psi1, p_psi2 )  RESULT(upward_tracer_flux)
        upward_tracer_flux(jc,thisLevel) = 0.5_wp *  &
          & ( w(jc,thisLevel)  * ( z_lext_1(thisLevel) + z_lext_2(thisLevel) )           &
          & +  ABS( w(jc,thisLevel) ) * ( z_lext_2(thisLevel) - z_lext_1(thisLevel) )    &
          & )

      END DO ! end loop over cells
    ENDDO ! end loop over vertical levels
    !
    ! set lower boundary condition
    !
    ! upward_tracer_flux(startIndex:endIndex,nlevp1) = 0.0_wp

    DO jc = startIndex, endIndex
! !CDIR NODEP
      DO thisLevel = firstLevel, cells_noOfLevels(jc)
        ! positive vertical divergence in direction of w (upward positive)
        flux_div_vert(jc,thisLevel) = upward_tracer_flux(jc, thisLevel) &
          & - upward_tracer_flux(jc, thisLevel+1)
      ENDDO
      DO thisLevel = cells_noOfLevels(jc)+1, n_zlev
        ! positive vertical divergence in direction of w (upward positive)
        flux_div_vert(jc,thisLevel) = 0.0_wp
      ENDDO
    END DO

  END SUBROUTINE upwind_vflux_ppm_onBlock
  !-------------------------------------------------------------------------
 
   
END MODULE mo_ocean_tracer_transport_vert

