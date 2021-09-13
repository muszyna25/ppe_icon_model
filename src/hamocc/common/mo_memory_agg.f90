!>
!! @brief variables for marine aggregate settling
!!
!! Definition of variables and allocation of memory
!!
MODULE mo_memory_agg

  USE mo_kind, ONLY : wp
  USE mo_control_bgc, ONLY: bgc_nproma, bgc_zlevs

  IMPLICIT NONE

  PUBLIC

  REAL(wp) :: AJ1, AJ2, AJ3, BJ1, BJ2, BJ3 ! constants for CD 

! for fix b (numbers distribution slope of aggregates)
  ! primary particle diameter for POM & PIM species involved in parametrized aggregation (m) 
  REAL(wp) :: dp_dust ! primary particle diameter dust
  REAL(wp) :: dp_det  ! primary particle diameter detritus
  REAL(wp) :: dp_calc ! primary particle diameter calc
  REAL(wp) :: dp_opal ! primary particle diameter opal
  REAL(wp) :: stickiness_tep  ! stickiness of TEP (related to opal frustules)
  REAL(wp) :: stickiness_det  ! normal detritus stickiness
  REAL(wp) :: stickiness_opal ! stickiness of opal (without TEP - just normal coating)
  REAL(wp) :: stickiness_calc ! stickiness of calc particles (coated with organics)
  REAL(wp) :: stickiness_dust ! stickiness of dust particles (coated with organics)
  REAL(wp) :: agg_df_max      ! maximum fractal dimension of aggregates (~2.5)
  REAL(wp) :: agg_df_min      ! minimum fractal dimension of aggregates (~1.2 - 1.6)
  REAL(wp) :: rho_tep         ! density of TEP particles
  REAL(wp) :: agg_Re_crit ! critical particle Reynolds number for nr-distribution limiting
  ! organic detritus density (alternative to orgdens to avoid negative ws)
  REAL(wp) :: agg_org_dens
  REAL(wp) :: det_mol2mass ! mol detritus P/m^3 to kg POM /m^3 (according to stoichiometry)

  REAL(wp),ALLOCATABLE :: av_dp(:,:),              &  ! mean primary particle diameter
                       &  av_rho_p(:,:),           &  ! mean primary particle density
                       &  df_agg(:,:),             &  ! fractal dimension of aggregates
                       &  b_agg(:,:),              &  ! aggregate number distribution slope
                       &  Lmax_agg(:,:),           &  ! maximum diamater of aggregates
                       &  ws_agg(:,:),             &  ! aggregate mean sinking velocity
                       &  stickiness_agg(:,:),     &  ! mean aggregate stickiness
                       &  stickiness_frustule(:,:),&  ! frustule stickiness
                       &  dynvis(:,:),             &  ! molecular dynamic viscosity
                       &  av_rhof_V(:,:) !,          &  ! volume-weighted aggregate density
!                       &  av_d_c(:,:),             &  ! concentration-weighted mean diameter of aggs
!                       &  av_por_V(:,:)               ! volume-weighted aggregate porosity

  INTEGER, PARAMETER :: &
       kavdp               =  1, &
       kavrhop             =  2, &
       kdfagg              =  3, &
       ksticka             =  4, &
       kLmaxagg            =  5, &
       kavrhof             =  6, &
!       kwsagg              =  7, &
!       kbagg               =  8, &
!       kstickf             =  9, &
!       kdynvis             = 10, &
!       kavdc               = 11, &
!       kavpor              = 12, &
       naggdiag            = 6

  REAL(wp), DIMENSION (:,:,:), ALLOCATABLE, TARGET :: aggdiag    ! 3d concentration EU

CONTAINS

  SUBROUTINE ALLOC_MEM_AGGREGATES
     !-----------------------------------------------------------------------
     !>
     !! Initialization/allocation fields
     !! Called in ini_bgc after read_namelist
     !!

     ! allocate memory space for aggregate properties
     ALLOCATE(av_dp(bgc_nproma,bgc_zlevs))
     ALLOCATE(av_rho_p(bgc_nproma,bgc_zlevs))
     ALLOCATE(df_agg(bgc_nproma,bgc_zlevs))
     ALLOCATE(b_agg(bgc_nproma,bgc_zlevs))
     ALLOCATE(Lmax_agg(bgc_nproma,bgc_zlevs))
     ALLOCATE(stickiness_agg(bgc_nproma,bgc_zlevs))
     ALLOCATE(stickiness_frustule(bgc_nproma,bgc_zlevs))
     ALLOCATE(av_rhof_V(bgc_nproma,bgc_zlevs))
!     ALLOCATE(av_d_C(bgc_nproma,bgc_zlevs))
!     ALLOCATE(av_por_V(bgc_nproma,bgc_zlevs))

     ALLOCATE(aggdiag(bgc_nproma,bgc_zlevs,naggdiag))

     ! mean sinking velocity
     ALLOCATE(ws_agg(bgc_nproma,bgc_zlevs))

     ! molecular dynamic viscosity
     ALLOCATE(dynvis(bgc_nproma,bgc_zlevs))

     av_dp = 0._wp
     av_rho_p = 0._wp
     df_agg = 0._wp
     b_agg = 0._wp
     Lmax_agg = 0._wp
     stickiness_agg = 0._wp
     stickiness_frustule = 0._wp
     av_rhof_V = 0._wp
!     av_d_C = 0._wp
!     av_por_V = 0._wp
     aggdiag = 0._wp

  END SUBROUTINE ALLOC_MEM_AGGREGATES

  SUBROUTINE CLEANUP_MEM_AGGREGATES

     DEALLOCATE(av_dp)
     DEALLOCATE(av_rho_p)
     DEALLOCATE(df_agg)
     DEALLOCATE(b_agg)
     DEALLOCATE(Lmax_agg)
     DEALLOCATE(stickiness_agg)
     DEALLOCATE(stickiness_frustule)
     DEALLOCATE(aggdiag)
     DEALLOCATE(ws_agg)
     DEALLOCATE(dynvis)
     DEALLOCATE(av_rhof_V)
!     DEALLOCATE(av_d_C)
!     DEALLOCATE(av_por_V)

  END SUBROUTINE CLEANUP_MEM_AGGREGATES

END MODULE mo_memory_agg
