!>
!! @brief Subroutine echam_phy_main calls all the parameterization schemes
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
!!  Modified for ICONAM by Marco Giorgetta (2014)
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

MODULE mo_interface_echam_radiation

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: pi
  USE mo_impl_constants      ,ONLY: min_rlcell_int, grf_bdywidth_c
  USE mo_run_config,          ONLY: ntracer, nlev, nlevp1
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, timer_radiation
  USE mtime,                  ONLY: datetime
  USE mo_psrad_radiation,     ONLY: psrad_radiation

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_model_domain        ,ONLY: t_patch
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_radiation

CONTAINS

  !----------------------------------------------------------------
  SUBROUTINE interface_echam_radiation( patch, field, this_datetime)
    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(datetime),  POINTER     :: this_datetime  !< time step

    INTEGER  :: itype(nproma)              !< type of convection
    LOGICAL  :: lglac(nproma)

    INTEGER  :: jg             
    INTEGER  :: i_nchdom, rl_start, rl_end
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block

    jg         = patch%id
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int
 
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    IF (ltimer) CALL timer_start(timer_radiation)

!$OMP PARALLEL DO PRIVATE(jcs,jce, lglac, itype)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      lglac(jcs:jce) = field%lfland(jcs:jce,jb) .AND.field%glac(jcs:jce,jb).GT.0._wp
      itype(jcs:jce) = NINT(field%rtype(jcs:jce,jb))

        CALL psrad_radiation(                       &
        & jg                                       ,&!< in  domain index
        & jb                                       ,&!< in  block index
        & kproma         = jce                     ,&!< in  end index for loop over block
        & kbdim          = nproma                   ,&!< in  dimension of block over cells
        & klev           = nlev                    ,&!< in  number of full levels = number of layers
        & klevp1         = nlevp1                  ,&!< in  number of half levels = number of layer interfaces
        & ktype          = itype(:)                ,&!< in  type of convection
        & loland         = field%lfland(:,jb)      ,&!< in  land-sea mask. (logical)
        & loglac         = lglac                   ,&!< in  glacier mask (logical)
        & this_datetime  = this_datetime           ,&!< in  actual time step
        & pcos_mu0       = field%cosmu0_rt(:,jb)   ,&!< in  solar zenith angle
        & alb_vis_dir    = field%albvisdir(:,jb)   ,&!< in  surface albedo for visible range, direct
        & alb_nir_dir    = field%albnirdir(:,jb)   ,&!< in  surface albedo for near IR range, direct
        & alb_vis_dif    = field%albvisdif(:,jb)   ,&!< in  surface albedo for visible range, diffuse
        & alb_nir_dif    = field%albnirdif(:,jb)   ,&!< in  surface albedo for near IR range, diffuse
        & tk_sfc         = field%ts_rad_rt(:,jb)   ,&!< in  grid box mean surface temperature
        & zf             = field%zf(:,:,jb)        ,&!< in  geometric height at full level      [m]
        & zh             = field%zh(:,:,jb)        ,&!< in  geometric height at half level      [m]
        & dz             = field%dz(:,:,jb)        ,&!< in  geometric height thickness of layer [m]
        & pp_hl          = field%presi_old(:,:,jb) ,&!< in  pressure at half levels at t-dt [Pa]
        & pp_fl          = field%presm_old(:,:,jb) ,&!< in  pressure at full levels at t-dt [Pa]
        & tk_fl          = field%ta(:,:,jb)        ,&!< in  tk_fl  = temperature at full level at t-dt
        & xm_dry         = field%mdry(:,:,jb)      ,&!< in  dry air mass in layer [kg/m2]
        & xm_trc         = field%mtrc(:,:,jb,:)    ,&!< in  tracer  mass in layer [kg/m2]
        & xm_ozn         = field%o3(:,:,jb)        ,&!< inout  ozone  mass mixing ratio [kg/kg]
        !
        & cdnc           = field% acdnc(:,:,jb)    ,&!< in   cloud droplet number conc
        & cld_frc        = field% aclc(:,:,jb)     ,&!< in   cloud fraction [m2/m2]
        & cld_cvr        = field%aclcov(:,jb)      ,&!< out  total cloud cover
        !
        & lw_dnw_clr     = field%rldcs_rt(:,:,jb)  ,&!< out  Clear-sky net longwave  at all levels
        & lw_upw_clr     = field%rlucs_rt(:,:,jb)  ,&!< out  Clear-sky net longwave  at all levels
        & sw_dnw_clr     = field%rsdcs_rt(:,:,jb)  ,&!< out  Clear-sky net shortwave at all levels
        & sw_upw_clr     = field%rsucs_rt(:,:,jb)  ,&!< out  Clear-sky net shortwave at all levels
        & lw_dnw         = field%rld_rt  (:,:,jb)  ,&!< out  All-sky net longwave  at all levels
        & lw_upw         = field%rlu_rt  (:,:,jb)  ,&!< out  All-sky net longwave  at all levels
        & sw_dnw         = field%rsd_rt  (:,:,jb)  ,&!< out  All-sky net longwave  at all levels
        & sw_upw         = field%rsu_rt  (:,:,jb)  ,&!< out  All-sky net longwave  at all levels
        !
        & vis_dn_dir_sfc = field%rvds_dir_rt(:,jb) ,&!< out  all-sky downward direct visible radiation at surface
        & par_dn_dir_sfc = field%rpds_dir_rt(:,jb) ,&!< all-sky downward direct PAR     radiation at surface
        & nir_dn_dir_sfc = field%rnds_dir_rt(:,jb) ,&!< all-sky downward direct near-IR radiation at surface
        & vis_dn_dff_sfc = field%rvds_dif_rt(:,jb) ,&!< all-sky downward diffuse visible radiation at surface
        & par_dn_dff_sfc = field%rpds_dif_rt(:,jb) ,&!< all-sky downward diffuse PAR     radiation at surface
        & nir_dn_dff_sfc = field%rnds_dif_rt(:,jb) ,&!< all-sky downward diffuse near-IR radiation at surface
        & vis_up_sfc     = field%rvus_rt    (:,jb) ,&!< all-sky upward visible radiation at surface
        & par_up_sfc     = field%rpus_rt    (:,jb) ,&!< all-sky upward PAR     radiation at surfac
        & nir_up_sfc     = field%rnus_rt    (:,jb) ) !< all-sky upward near-IR radiation at surface

    ENDDO
!$OMP END PARALLEL DO 
        
    IF (ltimer) CALL timer_stop(timer_radiation)

  END SUBROUTINE interface_echam_radiation
  !----------------------------------------------------------------

END MODULE mo_interface_echam_radiation
