!>
!! Contains the implementation of the vertical tracer transport routines for the ICON ocean model.
!! This comprises vertical advection.
!! 
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/01) 
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
MODULE mo_oce_tracer_transport_vert
!-------------------------------------------------------------------------  
!
!    ProTeX FORTRAN source: Style 2  
!    modified for ICON project, DWD/MPI-M 2006
!  
!-------------------------------------------------------------------------  
!
!
!
USE mo_kind,                      ONLY: wp
!USE mo_math_utilities,            ONLY: t_cartesian_coordinates
USE mo_impl_constants,            ONLY: sea_boundary
USE mo_math_constants,            ONLY: dbl_eps
USE mo_ocean_nml,                 ONLY: n_zlev,  &
  &                                     UPWIND, CENTRAL, MIMETIC_MIURA, ADPO, &
  &                                     FLUX_CALCULATION_VERT, l_vert_limiter_advection
USE mo_parallel_config,           ONLY: nproma
!USE mo_dynamics_config,           ONLY: nold, nnew
USE mo_run_config,                ONLY: dtime, ltimer
USE mo_timer,                     ONLY: timer_start, timer_stop, timer_adv_vert, timer_ppm_slim!, &
  !&                                     timer_dif_vert
USE mo_oce_state,                 ONLY: t_hydro_ocean_state!,v_base
USE mo_model_domain,              ONLY: t_patch,t_patch_3D
USE mo_exception,                 ONLY: finish !, message_text, message
USE mo_util_dbg_prnt,             ONLY: dbg_print
USE mo_oce_physics
USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
USE mo_sync,                      ONLY: SYNC_C, sync_patch_array
IMPLICIT NONE

PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'
CHARACTER(len=12)           :: str_module    = 'oceTracVert '  ! Output of module for 1 line debug
INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug

!
! PUBLIC INTERFACE
!
PUBLIC :: advect_flux_vertical
! Private implemenation
!
PRIVATE :: upwind_vflux_oce
PRIVATE :: central_vflux_oce
!PRIVATE :: mimetic_vflux_oce
PRIVATE :: upwind_vflux_ppm
PRIVATE :: v_ppm_slimiter_mo
PRIVATE :: laxfr_upflux_v

CONTAINS

  !-------------------------------------------------------------------------
  !! SUBROUTINE advects vertically the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
  !! mpi parallelized, sync required: trac_out
  SUBROUTINE advect_flux_vertical( p_patch_3D,   &
                                 & trac_old,             &
                                 & p_os,                 &
                                 & bc_top_tracer,        &
                                 & bc_bot_tracer,        &
                                 & flux_div_vert,        &
                                 & tracer_id)

    !TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
    TYPE(t_patch_3D ),TARGET          :: p_patch_3D
    REAL(wp), INTENT(INOUT)           :: trac_old(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET :: p_os
    REAL(wp)                          :: bc_top_tracer(nproma, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)                          :: bc_bot_tracer(nproma, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), INTENT(INOUT)           :: flux_div_vert(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) !new tracer
    INTEGER, INTENT(IN)               :: tracer_id

    !Local variables
    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: jc, jk, jb
    INTEGER  :: z_dolic
    REAL(wp) :: z_adv_flux_v (nproma, n_zlev+1, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) ! vertical advective tracer flux
    REAL(wp) :: z_adv_flux_vu(nproma, n_zlev+1, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) ! vertical upwind tracer flux
    REAL(wp) :: z_adv_flux_vc(nproma, n_zlev+1, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) ! vertical central tracer flux
    REAL(wp) :: trac_tst     (nproma, n_zlev  , p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) ! vertical advective tracer flux
    TYPE(t_patch), POINTER :: p_patch
    REAL(wp),      POINTER :: cell_area(:,:)

    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_tracer_advection:advect_individual_tracer')
    !-------------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-------------------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    cells_in_domain => p_patch%cells%in_domain
    cell_area       => p_patch%cells%area

    z_adv_flux_v  (1:nproma, 1:n_zlev+1, 1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp
    z_adv_flux_vu (1:nproma, 1:n_zlev+1, 1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp
    z_adv_flux_vc (1:nproma, 1:n_zlev+1, 1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp
    trac_tst      (1:nproma, 1:n_zlev  , 1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp

    CALL sync_patch_array(SYNC_C, p_patch, trac_old)

    ! This is already synced in  edges_in_domain !
    ! CALL sync_patch_array(SYNC_C, p_patch, p_os%p_diag%w_time_weighted)

    ! Initialize timer for vertical advection
    IF (ltimer) CALL timer_start(timer_adv_vert)

    !SELECT CASE(FLUX_CALCULATION_VERT)
    !CASE(UPWIND)

    IF (FLUX_CALCULATION_VERT == UPWIND .OR. FLUX_CALCULATION_VERT == ADPO) THEN

      CALL upwind_vflux_oce( p_patch_3D,                 &
                            & trac_old,                   &
                            & p_os%p_diag%w_time_weighted,&  ! in
                            & bc_top_tracer,              &
                            & z_adv_flux_vu,tracer_id )

      z_adv_flux_v (1:nproma, 1:n_zlev+1, 1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = &
      &  z_adv_flux_vu(1:nproma, 1:n_zlev+1, 1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    ENDIF

    !CASE(CENTRAL)
    IF (FLUX_CALCULATION_VERT == CENTRAL .OR. FLUX_CALCULATION_VERT == ADPO) THEN

      CALL central_vflux_oce(p_patch_3D,                 &
                           & trac_old,                   &
                           & p_os%p_diag%w_time_weighted,& ! in
                           & z_adv_flux_vc, tracer_id)

      !IF (l_vert_limiter_advection) &
      !& CALL vflx_limiter_pd_oce( p_patch,               &
      !                        & dtime,                   &
      !                        & trac_old,                &
      !                        & p_patch_3D%p_patch_1D(1)%prism_thick_c, &
      !                        & z_adv_flux_vc)

      z_adv_flux_v (1:nproma, 1:n_zlev+1, 1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = &
      &  z_adv_flux_vc(1:nproma, 1:n_zlev+1, 1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    ENDIF

    ! Vertical advection scheme: piecewise parabolic method (ppm)
    !CASE(MIMETIC_MIURA)
    IF (FLUX_CALCULATION_VERT == MIMETIC_MIURA) THEN

#ifndef __SX__
      !TODO review: here p_patch_3D%p_patch_1D(1)%prism_thick_c is used as a
      !data variable, whereas it should only contain geometrics constant in
      !time
      CALL upwind_vflux_ppm( p_patch_3D,                 &
                           & trac_old,                   &
                           & p_os%p_diag%w_time_weighted,& ! in
                           & dtime, 1 ,                  & 
                           & p_patch_3D%p_patch_1D(1)%prism_thick_c,  &
                           & z_adv_flux_v, tracer_id)
#endif
    ENDIF

    ! Vertical advection scheme: ADPO, adapted from MPIOM (Ernst Meier-Reimer)
    !   The implementation follows the MPIOM Draft documentation, section 5.2.13 describing ocadpo.f90
    IF (FLUX_CALCULATION_VERT == ADPO) THEN

      CALL adpo_vflux_oce( p_patch_3D,                                &
        &                    trac_old,                                &
        &                    p_os%p_diag%w_time_weighted,             &
        &                    dtime,                                   & 
        &                    p_patch_3D%p_patch_1D(1)%prism_thick_c,  &
        &                    trac_tst)

      z_adv_flux_v(1:nproma, 1:n_zlev+1, 1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) =                             &
        & z_adv_flux_vu(1:nproma, 1:n_zlev+1, 1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
  !     & z_adv_flux_vu(1:nproma, 1:n_zlev+1, 1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) * adpo_weight_upw_cntr + &
  !     & z_adv_flux_vc(1:nproma, 1:n_zlev+1, 1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) * adpo_weight_cntr_upw

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('AdvVertAdpo: adv_flux_v' ,z_adv_flux_v  ,str_module,idt_src)
    !---------------------------------------------------------------------

    ENDIF

    !CASE DEFAULT
    !  CALL finish('TRIM(advect_diffuse_flux_vert)',"This flux option is not supported")
    !END SELECT

    CALL sync_patch_array(SYNC_C, p_patch, z_adv_flux_v)

    !divergence is calculated for advective fluxes
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        DO jk = 1,p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
            ! positive vertical divergence in direction of w (upward positive)
            flux_div_vert(jc,jk,jb) = z_adv_flux_v(jc,jk,  jb) &
                                    &-z_adv_flux_v(jc,jk+1,jb)
        ENDDO
      END DO
    END DO

    IF (ltimer) CALL timer_stop(timer_adv_vert)

    CALL sync_patch_array(SYNC_C, p_patch, flux_div_vert)


    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('AdvVert: w_time_weighted',p_os%p_diag%w_time_weighted ,str_module,idt_src)
    !CALL dbg_print('AdvifVert: div_mass_flx_c' ,p_os%p_diag%div_mass_flx_c  ,str_module,idt_src)
    CALL dbg_print('AdvVert: adv_flux_v'     ,z_adv_flux_v                ,str_module,idt_src)
    CALL dbg_print('AdvVert: flux_div_vert'  ,flux_div_vert                ,str_module,idt_src)
    !---------------------------------------------------------------------
  END SUBROUTINE advect_flux_vertical
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !! First order upwind scheme for vertical tracer advection
  !!
  !! Calculation of vertical tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !! Initial revision by Jochen Foerstner, DWD (2008-05-15)
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized to height based vertical coordinate systems
  !! Modification by Stephan Lorenz, MPI (2010-09-07)
  !! - adapted to hydrostatic ocean core
  !!
  !! mpi parallelized, no sync
  SUBROUTINE upwind_vflux_oce( p_patch_3D, pvar_c, pw_c,top_bc_t, pupflux_i, tracer_id )

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    REAL(wp), INTENT(INOUT)           :: pvar_c(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)     !< advected cell centered variable
    REAL(wp), INTENT(INOUT)           :: pw_c(nproma,n_zlev+1, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)     !< in: vertical velocity on cells
    REAL(wp), INTENT(INOUT)           :: top_bc_t(nproma,      p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)     !< top boundary condition traver
    REAL(wp), INTENT(INOUT)           :: pupflux_i(nproma,n_zlev+1, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) !< variable in which the upwind flux is stored
    INTEGER, INTENT(IN)               :: tracer_id
    ! local variables
    ! height based but reversed (downward increasing depth) coordinate system,
    ! grid coefficient is negative (same as pressure based atmospheric coordinate system
    !REAL(wp), PARAMETER :: zcoeff_grid = -1.0_wp
    INTEGER             :: z_dolic
    INTEGER             :: i_startidx_c, i_endidx_c
    INTEGER             :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER             :: jkm1                     !< jk - 1
    INTEGER             :: k

    TYPE(t_patch), POINTER :: p_patch
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-------------------------------------------------------------------------
#ifndef __SX__
    p_patch         => p_patch_3D%p_patch_2D(1)
    cells_in_domain => p_patch%cells%in_domain

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        DO jk = 2, p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
            jkm1 = jk - 1
            ! calculate vertical tracer flux using upwind method
             pupflux_i(jc,jk,jb) =                 &
               &  laxfr_upflux_v( pw_c(jc,jk,jb),  &
               &                  pvar_c(jc,jkm1,jb), pvar_c(jc,jk,jb) )
        ENDDO
      END DO
    END DO 

#endif
  END SUBROUTINE upwind_vflux_oce
  !-------------------------------------------------------------------------
!   !>
!   !! First order upwind scheme for vertical tracer advection
!   !!
!   !! Calculation of vertical tracer fluxes 
!   !!
!   !! @par Revision History
!   !!!! mpi parallelized, no sync
!   SUBROUTINE mimetic_vflux_oce( p_patch, pvar_c, pw_c, pupflux_i, tracer_id )
! 
!     TYPE(t_patch), TARGET, INTENT(IN) :: p_patch           !< patch on which computation is performed
!     REAL(wp), INTENT(INOUT)           :: pvar_c(:,:,:)    !< advected cell centered variable
!     REAL(wp), INTENT(INOUT)           :: pw_c(:,:,:)      !< vertical velocity on cells
!     REAL(wp), INTENT(INOUT)           :: pupflux_i(:,:,:) !< variable in which the upwind flux is stored
!                                                           !< dim: (nproma,n_zlev+1,alloc_cell_blocks)
!     INTEGER, INTENT(IN)               :: tracer_id
!     ! local variables
!     INTEGER  :: i_startidx_c, i_endidx_c
!     INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
!     INTEGER  :: jkm1, jkp1, z_dolic                    !< jk - 1
!     REAL(wp) :: w_ave(n_zlev)
!     REAL(wp) :: w_avep1(n_zlev)
!     !-------------------------------------------------------------------------
!     TYPE(t_subset_range), POINTER :: cells_in_domain
!     !-------------------------------------------------------------------------
!     cells_in_domain => p_patch%cells%in_domain
! 
!     w_ave(:)  = 0.0_wp
!     w_avep1(:)= 0.0_wp
! 
!     DO jb = cells_in_domain%start_block, cells_in_domain%end_block
!       CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
!         DO jc = i_startidx_c, i_endidx_c
!           z_dolic = v_base%dolic_c(jc,jb)
!           IF(z_dolic>=MIN_DOLIC)THEN
!             DO jk = 2, z_dolic
!               jkm1 = jk - 1
!               jkp1 = jk + 1
!               IF(pw_c(jc,jk,jb)>pw_c(jc,jkm1,jb))THEN
! 
!                 w_ave(jk)=pw_c(jc,jk,jb)*pvar_c(jc,jk,jb)
!               ELSE
!                 w_ave(jk)=pw_c(jc,jk,jb)*pvar_c(jc,jkm1,jb)
!               ENDIF
!              !w_avep1(jk)=0.5_wp*(pw_c(jc,jk,jb)+pw_c(jc,jkp1,jb))&
!              !       &*pvar_c(jc,jk,jb)
! 
!             !pupflux_i(jc,jk,jb) = 0.5_wp*(w_ave(jk) +w_avep1(jk))
!             !pupflux_i(jc,jk,jb) = (w_ave(jk)*v_base%del_zlev_m(jk) &
!             !                    & +w_avep1(jk)*v_base%del_zlev_m(jkp1))&
!             !                   & /(v_base%del_zlev_m(jk)+v_base%del_zlev_m(jkp1))
!             END DO
!             DO jk=1,z_dolic-1
!               jkm1 = jk - 1
!               jkp1 = jk + 1
!               pupflux_i(jc,jk,jb) = (w_ave(jk)*v_base%del_zlev_m(jk) &
!                                  & +w_avep1(jk)*v_base%del_zlev_m(jkp1))&
!                                  & /(v_base%del_zlev_m(jk)+v_base%del_zlev_m(jkp1))
!             END DO
!             !w_ave(z_dolic)=0.5_wp*(pw_c(jc,z_dolic,jb)+pw_c(jc,z_dolic-1,jb))&
!             !        &*pvar_c(jc,z_dolic,jb)
!             !pupflux_i(jc,z_dolic,jb) = w_ave(z_dolic)*v_base%del_zlev_m(z_dolic)&
!             !                          &/v_base%del_zlev_i(z_dolic)
!             ! no fluxes at bottom boundary
!             pupflux_i(jc,z_dolic+1,jb) = 0.0_wp
!         ENDIF
!       END DO
!     END DO
! !     DO jb = all_cells%start_block, all_cells%end_block
! !       CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
! !         DO jc = i_startidx, i_endidx
! !           DO jk = 2, n_zlev
! !            ! index of top & bottom half level
! !            jkm1 = jk - 1
! !            jkp1 = jk + 1
! !              w_ave(jk)=0.5_wp*(pw_c(jc,jk,jb)+pw_c(jc,jkm1,jb))&
! !                    &*pvar_c(jc,jk,jb)
! !            END DO ! end cell loop
! !        END DO ! end level loop
! !      END DO ! end block loop
! !     DO jb = all_cells%start_block, all_cells%end_block
! !       CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
! !         DO jc = i_startidx, i_endidx
! !          DO jk = 1, n_zlev-1
! !            ! index of top & bottom half level
! !            jkm1 = jk - 1
! !            jkp1 = jk + 1
! !             w_avep1(jk)=0.5_wp*(pw_c(jc,jkp1,jb)+pw_c(jc,jk,jb))&
! !                    &*pvar_c(jc,jkp1,jb)
! !            END DO ! end cell loop
! !        END DO ! end level loop
! !      END DO ! end block loop
! !     DO jb = all_cells%start_block, all_cells%end_block
! !       CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
! !         DO jc = i_startidx, i_endidx 
! !          DO jk = 1, n_zlev-1
! !            ! index of top & bottom half level
! !            jkm1 = jk - 1
! !            jkp1 = jk + 1
! !            pupflux_i(jc,jk,jb) = 0.5_wp*(w_ave(jk)  *v_base%del_zlev_m(jk) &
! !                                &        +w_avep1(jk)*v_base%del_zlev_m(jkp1))/v_base%del_zlev_i(jk)
! !            END DO ! end cell loop
! !        END DO ! end level loop
! !      END DO ! end block loop
!   END SUBROUTINE mimetic_vflux_oce
  !-------------------------------------------------------------------------
  !>
  !!
  !! Calculation of central vertical tracer fluxes
  !!
  !! @par Revision History
  !! Petter Korn, MPI-M
  !!
  !! mpi parallelized, no sync
  SUBROUTINE central_vflux_oce( p_patch_3D, pvar_c, pw_c, c_flux_i, tracer_id )

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    REAL(wp), INTENT(INOUT)  :: pvar_c(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)     !< advected cell centered variable
    REAL(wp), INTENT(INOUT)  :: pw_c(nproma,n_zlev+1, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)     !< in: vertical velocity on cells
    REAL(wp), INTENT(INOUT)  :: c_flux_i(nproma,n_zlev+1, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) !< variable in which the central flux is stored
    INTEGER, INTENT(IN)      :: tracer_id
    ! local variables
    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: jkm1                     !< jk - 1
    INTEGER  :: z_dolic
    TYPE(t_patch), POINTER :: p_patch
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-------------------------------------------------------------------------
#ifndef __SX__
    p_patch         => p_patch_3D%p_patch_2D(1)
    cells_in_domain => p_patch%cells%in_domain

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)!v_base%dolic_c(jc,jb)
          !IF(z_dolic>=MIN_DOLIC)THEN
          DO jk = 2, z_dolic
            IF ( p_patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
              ! index of top half level
              jkm1 = jk - 1

              ! calculate vertical tracer flux using central method
              c_flux_i(jc,jk,jb) = 0.5_wp*pw_c(jc,jk,jb)              &
                &  * (pvar_c(jc,jkm1,jb)+ pvar_c(jc,jk,jb) )          !&
                !&          +0.5_wp*pw_c(jc,jk,jb)*pw_c(jc,jk,jb)*dtime&
                !&          / v_base%del_zlev_i(jk)                    &
                !&          *( pvar_c(jc,jkm1,jb) -pvar_c(jc,jk,jb))
            ENDIF 
           !zero advective flux at bottom boundary
           !c_flux_i(jc,z_dolic+1,jb)=0.0_wp

        ENDDO
      END DO
    END DO

#endif
  END SUBROUTINE central_vflux_oce

  !-------------------------------------------------------------------------
  !>
  !!
  !! Calculation of vertical tracer advection using ADPO-scheme of MPIOM
  !!
  !! This scheme uses a weighted average of upwind and central tracer advection
  !! schemes. The weights 'r' and 'R' are dependent on the first and second
  !! derivative of the tracer and the strength of the flow.
  !!
  !! The code is adopted from the MPIOM-model, written by Ernst Meier-Reimer
  !!
  !! LITERATURE
  !!  - Sweby (1984), SIAM-JNA, 21, 995-1011
  !!  - MPIOM Draft Documentation (Wetzel, Haak, Jungclaus, Meier-Reimer)
  !!
  !! @par Revision History
  !! Initial revision by Stephan Lorenz, MPI-M (2013-10)
  !! Stephan Lorenz, MPI-M
  !!
  !! not yet mpi parallelized, no sync
  SUBROUTINE adpo_vflux_oce ( p_patch_3D,  &
    &                         p_cc,        & ! advected cell centered variable
    &                         p_w,         & ! vertical velocity
    &                         p_dtime,     & ! time step
    &                         p_thick_c,   & ! layer thickness
    &                         p_trac_out)    ! output tracer field, tracer flux

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    REAL(wp), INTENT(INOUT)  :: p_cc        (nproma,n_zlev  , p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), INTENT(INOUT)  :: p_w         (nproma,n_zlev+1, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), INTENT(IN)     :: p_dtime                                                                  
    REAL(wp), INTENT(INOUT)  :: p_thick_c   (nproma,n_zlev, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)  
    REAL(wp), INTENT(INOUT)  :: p_trac_out  (nproma,n_zlev, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)  

    ! local variables
    INTEGER  :: jc, jk, jb, i_startidx_c, i_endidx_c
    REAL(wp) :: deriv_fst, deriv_sec!, adpo_weight_upw_cntr, adpo_weight_cntr_upw
    REAL(wp) :: prism_volume, new_volume, adpo_r1, adpo_r2, adpo_r3
    REAL(wp) :: advvol_in, advvol_out, advtrc_in, advtrc_out, wupw_in, wupw_out
    REAL(wp) :: z_adpo_vol_in (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: z_adpo_vol_out(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: z_adpo_trc_in (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: z_adpo_trc_out(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_patch),        POINTER :: p_patch
    REAL(wp),             POINTER :: cell_area(:,:)
    TYPE(t_subset_range), POINTER :: cells_in_domain

    p_patch         => p_patch_3D%p_patch_2D(1)
    cells_in_domain => p_patch%cells%in_domain
    cell_area       => p_patch%cells%area

    z_adpo_vol_in (1:nproma, 1:n_zlev,   1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp
    z_adpo_vol_out(1:nproma, 1:n_zlev,   1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp
    z_adpo_trc_in (1:nproma, 1:n_zlev,   1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp
    z_adpo_trc_out(1:nproma, 1:n_zlev,   1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp

    ! special treatment of top and bottom
    !  - no calculation of r and R
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c

        jk = 1
        IF ( p_patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN

          ! calculation of volume transport, same as below with top boundary condition
          wupw_in    =     p_w(jc,1,jb) ! zero in jk=1: + ABS(p_w(jc,jk-1,jb))
          wupw_out   = ABS(p_w(jc,2,jb)) - p_w(jc,2,jb)
          advvol_in  = wupw_in  * 0.5_wp * p_dtime * cell_area(jc,jb)
          advvol_out = wupw_out * 0.5_wp * p_dtime * cell_area(jc,jb)
          z_adpo_vol_in (jc,1,jb) = advvol_in
          z_adpo_vol_out(jc,1,jb) = advvol_out
    
          ! calculate tracer transport up/downward (TRP/TRM in MPIOM) at top level
       !  z_adpo_trc_in (jc, 1,jb) = 0.0_wp  !  not used, since no trc_in(jk-1) is needed below
          z_adpo_trc_out(jc, 1,jb) = advvol_out * p_cc(jc,1,jb) ! central part ignored at top

        ENDIF

        ! calculation of volume transport, same as below with bottom boundary condition
        jk = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
        IF ( jk > 1 ) THEN  !  sea-point

          wupw_in    =     p_w(jc,jk  ,jb)  + ABS(p_w(jc,jk-1,jb))
        ! wupw_in    =     p_w(jc,jk  ,jb)  + ABS(p_w(jc,jk  ,jb)) ! MPIOM
          wupw_out   = ABS(p_w(jc,jk+1,jb)) -     p_w(jc,jk+1,jb)
          advvol_in  = wupw_in  * 0.5_wp * p_dtime * cell_area(jc,jb)
          advvol_out = wupw_out * 0.5_wp * p_dtime * cell_area(jc,jb)
          z_adpo_vol_in (jc,jk,jb) = advvol_in
          z_adpo_vol_out(jc,jk,jb) = advvol_out
          z_adpo_trc_in (jc,jk,jb) = advvol_in  * p_cc(jc,jk,jb)  ! central part ignored at bottom
          z_adpo_trc_out(jc,jk,jb) = advvol_out * p_cc(jc,jk,jb)  ! central part ignored at bottom
          ! Remark: in MPIOM the central part is only ignored if ocean bottom is at k=kend

        ENDIF

      ENDDO
    ENDDO

    ! calculate weighting factors for part of upwind and central schemes: main vertical loop
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        DO jk = 2,p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)-1

          ! calculation of first and second spatial derivative
          prism_volume  = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,jk,jb) * cell_area(jc,jb)
          deriv_fst = p_cc(jc,jk-1,jb)-p_cc(jc,jk+1,jb)
          deriv_sec = p_cc(jc,jk-1,jb)-p_cc(jc,jk+1,jb) - 2.0_wp*p_cc(jc,jk,jb)

          ! weighting factor corresponding to "r" in MPIOM documentation
          adpo_r1 = MAX(0.0_wp, (deriv_fst-deriv_sec) / (deriv_fst + 1.E-20_wp))

          ! weighting factor corresponding to "R" in MPIOM documentation
          !  - this calculation is independent of tracers, see ocadpo_base in MPIOM
          !  - z_adpo_vol_in  is water transport wtp in mpiom, is U_in  in MPIOM documentation
          !  - z_adpo_vol_out is water transport wtm in mpiom, is U_out in MPIOM documentation
          !  - these transports can be replaced by upwind flux z_adv_flux_vu
          wupw_in    =     p_w(jc,jk  ,jb)  + ABS(p_w(jc,jk-1,jb))
          wupw_out   = ABS(p_w(jc,jk+1,jb)) -     p_w(jc,jk+1,jb)
          ! ATTENTION - implementation in MPIOM-model is different from MPIOM documentation:
          !  - wupw_in (wtk in MPIOM) is taken as simple upwind part using level k only
        ! wupw_in    =     p_w(jc,jk  ,jb)  + ABS(p_w(jc,jk  ,jb))

          ! multiplying with factors
          advvol_in  = wupw_in  * 0.5_wp * p_dtime * cell_area(jc,jb)
          advvol_out = wupw_out * 0.5_wp * p_dtime * cell_area(jc,jb)
          z_adpo_vol_in (jc,jk,jb) = advvol_in
          z_adpo_vol_out(jc,jk,jb) = advvol_out
          adpo_r2 = MIN(1.0_wp, &
            &       adpo_r1*prism_volume / (advvol_in + advvol_out + 1.E-20_wp)) * &
            &       p_patch_3D%wet_c(jc,jk+1,jb)  ! set zero on land
          adpo_r3 = 1.0_wp - adpo_r2

          ! calculate resulting vertical advection fluxes: main vertical loop
          !  - these transports can be replaced by weighted upwind and central flux z_adv_flux_vu/vc
          advtrc_in  = adpo_r3*advvol_in *p_cc(jc,jk,jb) + adpo_r2*0.5_wp*(p_cc(jc,jk,jb)+p_cc(jc,jk-1,jb))
          advtrc_out = adpo_r3*advvol_out*p_cc(jc,jk,jb) + adpo_r2*0.5_wp*(p_cc(jc,jk,jb)+p_cc(jc,jk+1,jb))
          z_adpo_trc_in (jc,jk,jb) = advtrc_in
          z_adpo_trc_out(jc,jk,jb) = advtrc_out

        ENDDO
      ENDDO
    ENDDO

    ! calculate resulting new tracer concentration
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c

        !  - main vertical loop
        DO jk = 2,p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)-1

       !  write(*,*) 'XXXXX indices: jb,jc,jk,dolic = ',jb,jc,jk,p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)-1
          prism_volume  = p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,jk,jb) * cell_area(jc,jb)
          new_volume    = prism_volume + z_adpo_vol_in  (jc,jk+1,jb) - z_adpo_vol_in  (jc,jk,jb) &
            &                          + z_adpo_vol_out (jc,jk-1,jb) - z_adpo_vol_out (jc,jk,jb)
       !  write(*,*) 'XXXXX volumes                 = ', prism_volume, new_volume
       !  write(*,*) 'XXXXX adpo_trcin (jk+1) = ', z_adpo_trc_in(jc,jk+1,jb)
       !  write(*,*) 'XXXXX adpo_trcin (jk)   = ', z_adpo_trc_in(jc,jk  ,jb)
       !  write(*,*) 'XXXXX adpo_trcout(jk)   = ', z_adpo_trc_out(jc,jk  ,jb)
       !  write(*,*) 'XXXXX adpo_trcout(jk-1) = ', z_adpo_trc_out(jc,jk-1,jb)

          ! - update new tracer values
          p_trac_out(jc,jk,jb) = (p_cc(jc,jk,jb)*prism_volume &
                                       + z_adpo_trc_in (jc,jk+1,jb) - z_adpo_trc_in (jc,jk,jb)  &
            &                          + z_adpo_trc_out(jc,jk-1,jb) - z_adpo_trc_out(jc,jk,jb)) &
            &                   /new_volume

        ENDDO
      ENDDO
    ENDDO

    ! special treatment of top and bottom levels
    !  - no calculation of r and R
    !  - this loop can be integrated above
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c

        jk = 1
        IF ( p_patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN

          prism_volume  = p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,jk,jb) * cell_area(jc,jb)
          ! new_volume consists of 3 summands: wtm(i,j,k-1)=z_adpo_vol_out(jc,jk-1,jb)=0
          new_volume    = prism_volume + z_adpo_vol_in  (jc,jk+1,jb) - z_adpo_vol_in  (jc,jk,jb) &
            &                                                        - z_adpo_vol_out (jc,jk,jb)

          ! - update new tracer values: T_in(k) and T_out(k-1) ignored - no remark in MPIOM documentation
          p_trac_out(jc,jk,jb) = (p_cc(jc,jk,jb)*prism_volume &
                                       + z_adpo_trc_in (jc,jk+1,jb) - z_adpo_trc_out(jc,jk,jb)) &
            &                   /new_volume

        ENDIF

        ! - special treatment of bottom
        jk = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
        IF ( jk > 1 ) THEN  !  sea-point

          prism_volume  = p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,jk,jb) * cell_area(jc,jb)
          ! new_volume consists of 3 summands: wtp(i,j,k+1)=z_adpo_vol_in(jc,jk+1,jb)=0
          new_volume    = prism_volume                               - z_adpo_vol_in  (jc,jk,jb) &
            &                          + z_adpo_vol_out (jc,jk-1,jb) - z_adpo_vol_out (jc,jk,jb)

          ! - update new tracer values: T_in(k+1) ignored
          p_trac_out(jc,jk,jb) = (p_cc(jc,jk,jb)*prism_volume &
                                                                    - z_adpo_trc_in (jc,jk,jb)  &
            &                          + z_adpo_trc_out(jc,jk-1,jb) - z_adpo_trc_out(jc,jk,jb)) &
            &                   /new_volume

        ENDIF

      ENDDO
    ENDDO

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('AdvVertAdpo: adpo_vol_in'   ,z_adpo_vol_in ,str_module,idt_src)
    CALL dbg_print('AdvVertAdpo: adpo_vol_out'  ,z_adpo_vol_out,str_module,idt_src)
    CALL dbg_print('AdvVertAdpo: adpo_trc_in '  ,z_adpo_trc_in ,str_module,idt_src)
    CALL dbg_print('AdvVertAdpo: adpo_trc_out'  ,z_adpo_trc_out,str_module,idt_src)
    CALL dbg_print('AdvVertAdpo: trac_out'      ,p_trac_out    ,str_module,idt_src)
    !---------------------------------------------------------------------

  END SUBROUTINE adpo_vflux_oce

  !------------------------------------------------------------------------
  !! The third order PPM scheme
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
  SUBROUTINE upwind_vflux_ppm( p_patch_3D, p_cc,        &
    &                      p_w, p_dtime, p_itype_vlimit,&
    &                      p_cellhgt_mc_now, &
    &                      p_upflux, tracer_id)

!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
!!$      &  routine = 'mo_advection_vflux: upwind_vflux_ppm'

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    REAL(wp), INTENT(INOUT)           :: p_cc(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)  !< advected cell centered variable
    REAL(wp), INTENT(INOUT)           :: p_w(nproma,n_zlev+1, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) !<  in, but synced: vertical velocity
    REAL(wp), INTENT(IN)              :: p_dtime  !< time step
    REAL(wp), INTENT(INOUT)           :: p_cellhgt_mc_now(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)!< layer thickness at cell center at time n
    REAL(wp), INTENT(INOUT)           :: p_upflux(nproma,n_zlev+1, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)      !< output field, tracer flux
    INTEGER, INTENT(IN)               :: p_itype_vlimit                                  !< parameter to select limiter
    INTEGER, INTENT(IN)               :: tracer_id
!
!local variables
    REAL(wp) :: z_face(nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)   !< face values of transported field
    REAL(wp) :: z_face_up(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)  !< face value (upper face)
    REAL(wp) :: z_face_low(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) !< face value (lower face)
    REAL(wp) :: z_lext_1(nproma,n_zlev+1)                 !< linear extrapolation value 1 
    REAL(wp) :: z_lext_2(nproma,n_zlev+1)                 !< linear extrapolation value 2
    REAL(wp) :: z_cfl_m(nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)  !< CFL number (weta>0, w<0)
    REAL(wp) :: z_cfl_p(nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)  !< CFL number (weta<0, w>0)
    REAL(wp) :: z_slope(nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)  !< monotonized slope
    REAL(wp) :: z_slope_u, z_slope_l                            !< one-sided slopes
    REAL(wp) :: z_delta_m, z_delta_p                            !< difference between lower and upper face value
                                                                !< for weta >0 and weta <0
    REAL(wp) :: z_a11, z_a12                                    !< 1/6 * a6,i (see Colella and Woodward (1984))
    REAL(wp) :: z_weta_dt                                       !< weta times p_dtime
    INTEGER  :: slev, slevp1                                    !< vertical start level and start level +1
    INTEGER  :: nlevp1                                          !< number of full and half levels
    INTEGER  :: ikm1, ikp1, ikp1_ic,ikp2                        !< vertical level minus and plus one, plus two
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: jc, jk, jb                                      !< index of cell, vertical level and block
    !LOGICAL  :: opt_lout_edge !< optional: output edge value (.TRUE.),
    !                          !< or the flux across the edge   !< (.FALSE./not specified)
    !REAL(wp) :: opt_topflx_tra(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)  !< vertical tracer flux at upper boundary
    INTEGER, PARAMETER :: islopel_vsm = 1   
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-----------------------------------------------------------------------
#ifndef __SX__
    p_patch         => p_patch_3D%p_patch_2D(1)
    cells_in_domain => p_patch%cells%in_domain

    slev   = 1
    slevp1 = 2
    nlevp1 = n_zlev+1

    z_cfl_m   (1:nproma,1:n_zlev+1,1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp
    z_cfl_p   (1:nproma,1:n_zlev+1,1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp
    z_face    (1:nproma,1:n_zlev+1,1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp
    z_face_low(1:nproma,1:n_zlev,  1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp
    z_face_up (1:nproma,1:n_zlev,  1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp
    z_slope   (1:nproma,1:n_zlev,  1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp
    z_lext_1  (1:nproma,1:n_zlev+1)=0.0_wp 
    z_lext_2  (1:nproma,1:n_zlev+1)=0.0_wp


    CALL sync_patch_array(SYNC_C, p_patch, p_cc)
    CALL sync_patch_array(SYNC_C, p_patch, p_cellhgt_mc_now)

    ! this is noit needed
    ! CALL sync_patch_array(SYNC_C, p_patch, p_w)


    ! advection is done with an upwind scheme and a piecwise parabolic
    ! approx. of the subgrid distribution is used.
    ! 3 options:  standard without limiter
    !             standard with semi-monotone or monotone limiter
    !             special version with limiter which handles CFL >1
    !
    ! 1. Calculate Courant number for weta>0 (w<0) and weta<0 (w>0)
    !
! !$OMP PARALLEL
! !$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,ikm1,z_weta_dt,ikp1_ic,ikp1, &
! !$OMP            z_slope_u,z_slope_l,ikp2)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)

      ! Courant number at top and bottom
      z_cfl_p(i_startidx:i_endidx,slev,jb)   = 0._wp
      z_cfl_m(i_startidx:i_endidx,slev,jb)   = 0._wp
      z_cfl_p(i_startidx:i_endidx,nlevp1,jb) = 0._wp
      z_cfl_m(i_startidx:i_endidx,nlevp1,jb) = 0._wp

      DO jk = slevp1, n_zlev
        ! index of top half level
        ikm1 = jk - 1
        DO jc = i_startidx, i_endidx
          ! Calculate local Courant number at half levels
          ! z_cfl_m for weta >0 (w <0)
          ! z_cfl_p for weta <0 (w >0)
          !z_weta_dt = 0.0_wp

            IF ( p_patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
               z_weta_dt = ABS(p_w(jc,jk,jb)) * p_dtime
 
              z_cfl_m(jc,jk,jb) = z_weta_dt &
              &/ p_cellhgt_mc_now(jc,ikm1,jb) !p_patch_3D%p_patch_1D(1)%del_zlev_m(ikm1)!
              z_cfl_p(jc,jk,jb) = z_weta_dt &
              &/ p_cellhgt_mc_now(jc,jk,jb)!p_patch_3D%p_patch_1D(1)%del_zlev_m(jk)!
            ENDIF
        END DO ! end loop over cells
      ENDDO ! end loop over vertical levels
      !
      ! 2. Calculate monotonized slope
      !
      z_slope(i_startidx:i_endidx,slev,jb) = 0._wp

      DO jk = slevp1, n_zlev
        ! index of top half level
        ikm1    = jk - 1
        ! index of bottom half level
        ikp1_ic = jk + 1
        ikp1    = MIN( ikp1_ic, n_zlev )

        DO jc = i_startidx, i_endidx
          !IF ( p_patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            z_slope_u = 2._wp * (p_cc(jc,jk,jb) - p_cc(jc,ikm1,jb))
          !ELSE
          !  z_slope_u = 0.0_wp
          !ENDIF
          !IF ( p_patch_3D%lsm_c(jc,ikp1,jb) <= sea_boundary ) THEN
            z_slope_l = 2._wp * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))
          !ELSE
          !  z_slope_l = 0.0_wp
          !ENDIF

          IF ((z_slope_u * z_slope_l) > 0._wp) THEN

            z_slope(jc,jk,jb) = ( p_cellhgt_mc_now(jc,jk,jb)                             &
              &  / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)            &
              &  + p_cellhgt_mc_now(jc,ikp1,jb)) )                                       &
              &  * ( (2._wp * p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)) &
              &  / (p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,jk,jb))           &
              &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))                                   &
              &  + (p_cellhgt_mc_now(jc,jk,jb) + 2._wp * p_cellhgt_mc_now(jc,ikp1,jb))   &
              &  / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb))           &
              &  * (p_cc(jc,jk,jb) - p_cc(jc,ikm1,jb)) )

            z_slope(jc,jk,jb) = SIGN(                                            &
              &  MIN( ABS(z_slope(jc,jk,jb)), ABS(z_slope_u), ABS(z_slope_l) ),  &
              &    z_slope(jc,jk,jb))
          ELSE
            z_slope(jc,jk,jb) = 0._wp
          ENDIF
        END DO ! end loop over cells
      END DO ! end loop over vertical levels

      !
      ! 3. reconstruct face values at vertical half-levels
      !
      ! Boundary values for two highest and lowest half-levels
      !
      ! for faces k=slevp1 and k=nlevp1-1 reconstructed face values are calculated by
      ! interpolating a quadratic (instead of quartic) polynomial through 3
      ! values of the indefinite integral A=\int_{\eta_{0}}^{\eta}q\,\mathrm{d}\eta
      !
      ! for faces k=slev and k=nlevp1 a zero gradient condition is assumed and the
      ! face values are set to the tracer values of the corresponding cell centers
      !
      DO jc = i_startidx, i_endidx
        IF ( p_patch_3D%lsm_c(jc,slevp1,jb) <= sea_boundary ) THEN
        z_face(jc,slevp1,jb) = p_cc(jc,slev,jb)*(1._wp - (p_cellhgt_mc_now(jc,slev,jb)&
          &       / p_cellhgt_mc_now(jc,slevp1,jb))) + (p_cellhgt_mc_now(jc,slev,jb)  &
          &       /(p_cellhgt_mc_now(jc,slev,jb) + p_cellhgt_mc_now(jc,slevp1,jb)))   &
          &       * ((p_cellhgt_mc_now(jc,slev,jb) / p_cellhgt_mc_now(jc,slevp1,jb))  &
          &       * p_cc(jc,slev,jb) + p_cc(jc,slevp1,jb))
        ELSE
        z_face(jc,slevp1,jb) = 0.0_wp
        ENDIF 
        IF ( p_patch_3D%lsm_c(jc,n_zlev,jb) <= sea_boundary ) THEN
        z_face(jc,n_zlev,jb) = p_cc(jc,n_zlev-1,jb)*( 1._wp                               &
          &       - (p_cellhgt_mc_now(jc,n_zlev-1,jb) / p_cellhgt_mc_now(jc,n_zlev,jb)))  &
          &       + (p_cellhgt_mc_now(jc,n_zlev-1,jb)/(p_cellhgt_mc_now(jc,n_zlev-1,jb)   &
          &       + p_cellhgt_mc_now(jc,n_zlev,jb))) * ((p_cellhgt_mc_now(jc,n_zlev-1,jb) &
          &       / p_cellhgt_mc_now(jc,n_zlev,jb)) * p_cc(jc,n_zlev-1,jb)                &
          &       + p_cc(jc,n_zlev,jb))
        ELSE
          z_face(jc,n_zlev,jb) = 0.0_wp
        ENDIF
        IF ( p_patch_3D%lsm_c(jc,slev,jb) <= sea_boundary ) THEN
          z_face(jc,slev,jb) = p_cc(jc,slev,jb)
        ELSE
          z_face(jc,slev,jb) = 0.0_wp
        ENDIF
        IF ( p_patch_3D%lsm_c(jc,n_zlev,jb) <= sea_boundary ) THEN
          z_face(jc,nlevp1,jb) = p_cc(jc,n_zlev,jb)
        ELSE
          z_face(jc,nlevp1,jb) = 0.0_wp
        ENDIF
      ENDDO


      DO jk = slevp1, n_zlev-2
        ! index of top half level
        ikm1 = jk - 1
        ! index of bottom half level
        ikp1 = jk + 1
        ikp2 = jk + 2
        DO jc = i_startidx, i_endidx
          IF ( p_patch_3D%lsm_c(jc,ikp2,jb) <= sea_boundary ) THEN
          z_face(jc,ikp1,jb) = p_cc(jc,jk,jb) &
            &  + (p_cellhgt_mc_now(jc,jk,jb)                                          &
            &  / (p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb)))         &
            &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))                                  &
            &  + (1._wp/(p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)    &
            &  + p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,ikp2,jb)))        &
            &  * ( (2._wp * p_cellhgt_mc_now(jc,ikp1,jb) * p_cellhgt_mc_now(jc,jk,jb) &
            &  / (p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb)))         &
            &  * ( (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb))        &
            &  / (2._wp*p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb))    &
            &  - (p_cellhgt_mc_now(jc,ikp2,jb) + p_cellhgt_mc_now(jc,ikp1,jb))        &
            &  / (2._wp*p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,jk,jb)) )  &
            &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb)) - p_cellhgt_mc_now(jc,jk,jb)     &
            &  * z_slope(jc,ikp1,jb) * (p_cellhgt_mc_now(jc,ikm1,jb)                  &
            &  + p_cellhgt_mc_now(jc,jk,jb)) / (2._wp*p_cellhgt_mc_now(jc,jk,jb)      &
            &  + p_cellhgt_mc_now(jc,ikp1,jb)) + p_cellhgt_mc_now(jc,ikp1,jb)         &
            &  * z_slope(jc,jk,jb) * (p_cellhgt_mc_now(jc,ikp1,jb)                    &
            &  + p_cellhgt_mc_now(jc,ikp2,jb)) / (p_cellhgt_mc_now(jc,jk,jb)          &
            &  + 2._wp*p_cellhgt_mc_now(jc,ikp1,jb)) )
           ELSE
           z_face(jc,ikp1,jb) = 0.0_wp
           ENDIF
        END DO ! end loop over cells
      END DO ! end loop over vertical levels

    END DO
! !$OMP END DO
! !$OMP END PARALLEL

    CALL sync_patch_array(SYNC_C, p_patch, z_face)
    CALL sync_patch_array(SYNC_C, p_patch, z_face_up)
    CALL sync_patch_array(SYNC_C, p_patch, z_face_low)
    !CALL sync_patch_array(SYNC_C, p_patch, z_slope)

    !
    ! 4. Limitation of first guess parabola (which is based on z_face)
    ! Note that z_face_up(k) does not need to equal z_face_low(k-1) after
    ! the limitation procedure.
    ! Therefore 2 additional fields z_face_up and z_face_low are
    ! introduced.
    !
    IF (l_vert_limiter_advection) THEN
     IF (p_itype_vlimit == islopel_vsm) THEN
!      ! monotonic (mo) limiter
      IF (ltimer) CALL timer_start(timer_ppm_slim)

      CALL v_ppm_slimiter_mo( p_patch_3D, &
                            & p_cc,       &
                            & z_face,     &
                            & z_slope,    &
                            & z_face_up,  &
                            & z_face_low)

      IF (ltimer) CALL timer_stop(timer_ppm_slim)
    ELSE
!      ! simply copy face values to 'face_up' and 'face_low' arrays
! !$OMP PARALLEL
! !$OMP DO PRIVATE(jk,ikp1,jb,i_startidx,i_endidx)
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
        DO jk = slev, n_zlev
          ! index of bottom half level
          ikp1 = jk + 1
          !IF ( v_base%lsm_c(jc,ikp1,jb) <= sea_boundary ) THEN
          IF ( p_patch_3D%lsm_c(jc,ikp1,jb) <= sea_boundary ) THEN
          z_face_up(i_startidx:i_endidx,jk,jb)  = z_face(i_startidx:i_endidx,jk,jb)
          z_face_low(i_startidx:i_endidx,jk,jb) = z_face(i_startidx:i_endidx,ikp1,jb)
          ELSE
          z_face_up(i_startidx:i_endidx,jk,jb)  = 0.0_wp
          z_face_low(i_startidx:i_endidx,jk,jb) = 0.0_wp
          ENDIF
        ENDDO
      ENDDO
! !$OMP ENDDO
! !$OMP END PARALLEL
     ENDIF  !  p_ityp_vlimit
    ENDIF  ! l_vert_limiter_advection


! !$OMP PARALLEL
    ! 5b. extrapolation using piecewise parabolic approx. of the transported
    ! quantity to the edge and finally, calculation of the upwind fluxes
    !

! !$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_lext_1,z_lext_2,ikm1,z_delta_m, &
! !$OMP            z_delta_p,z_a11,z_a12)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)

      z_lext_1(i_startidx:i_endidx,slev)   = p_cc(i_startidx:i_endidx,slev,jb)
      z_lext_2(i_startidx:i_endidx,slev)   = p_cc(i_startidx:i_endidx,slev,jb)
      z_lext_1(i_startidx:i_endidx,nlevp1) = p_cc(i_startidx:i_endidx,n_zlev,jb)
      z_lext_2(i_startidx:i_endidx,nlevp1) = p_cc(i_startidx:i_endidx,n_zlev,jb)

      DO jk = slevp1, n_zlev
        ! index of top half level
        ikm1 = jk -1
        DO jc = i_startidx, i_endidx
        IF ( p_patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN          
          ! linear extrapolated values
          ! for the height based coordinate system multiplication by coeff_grid
          ! is not necessary due to compensating (-) signs.
          ! first (of cell above) (case of w < 0; weta > 0)
          z_delta_m = z_face_low(jc,ikm1,jb) - z_face_up(jc,ikm1,jb)
          z_a11     = p_cc(jc,ikm1,jb)                                  &
            &       - 0.5_wp * (z_face_low(jc,ikm1,jb) + z_face_up(jc,ikm1,jb))

          z_lext_1(jc,jk) = p_cc(jc,ikm1,jb)                            &
            &  + (0.5_wp * z_delta_m * (1._wp - z_cfl_m(jc,jk,jb)))     &
            &  - z_a11*(1._wp - 3._wp*z_cfl_m(jc,jk,jb)                 &
            &  + 2._wp*z_cfl_m(jc,jk,jb)*z_cfl_m(jc,jk,jb))

          ! second (of cell below) (case of w > 0; weta < 0)
          z_delta_p = z_face_low(jc,jk,jb) - z_face_up(jc,jk,jb)
          z_a12     = p_cc(jc,jk,jb)                                    &
            &       - 0.5_wp * (z_face_low(jc,jk,jb) + z_face_up(jc,jk,jb))

          z_lext_2(jc,jk) = p_cc(jc,jk,jb)                              &
            &  - (0.5_wp * z_delta_p * (1._wp - z_cfl_p(jc,jk,jb)))     &
            &  - z_a12*(1._wp - 3._wp*z_cfl_p(jc,jk,jb)                 &
            &  + 2._wp*z_cfl_p(jc,jk,jb)*z_cfl_p(jc,jk,jb))
          !
          ! calculate vertical tracer flux
          !
          p_upflux(jc,jk,jb) =                                  &
            &  laxfr_upflux_v( p_w(jc,jk,jb),       &
            &                z_lext_1(jc,jk), z_lext_2(jc,jk))
        ELSE
          p_upflux(jc,jk,jb) = 0.0_wp
        ENDIF
        END DO ! end loop over cells
      ENDDO ! end loop over vertical levels
      !
      ! set lower boundary condition
      !
      p_upflux(i_startidx:i_endidx,nlevp1,jb) = 0.0_wp
      !
    ENDDO ! end loop over blocks
! !$OMP END DO
! !$OMP END PARALLEL


      CALL sync_patch_array(SYNC_C, p_patch, p_upflux)


      ! 6. If desired, apply a flux limiter to limit computed fluxes.
      !    These flux limiters are based on work by Zalesak (1979)
      !    positive-definite (pd) limiter      !
      !    IF (p_itype_vlimit == ifluxl_vpd) THEN
      !     IF (l_vert_limiter_advection) &
      !       & CALL vflx_limiter_pd_oce( p_patch, p_dtime, p_cc, p_cellhgt_mc_now, p_upflux)
      !   ENDIF

    CALL sync_patch_array(SYNC_C, p_patch, p_upflux)
#endif
  END SUBROUTINE upwind_vflux_ppm
 !-------------------------------------------------------------------------
  !>
  !! Limiter for PPM (3rd order) vertical advection (monotone version)
  !!
  !! Removes over- and undershoots in first guess parabola by resetting the
  !! upper or lower interface values.
  !! Avoids non-physical over/undershoots in advected fields.
  !!
  !! Note that this limiter was coded assuming a pressure based vertical
  !! coordinate system. Nevertheless this limiter works for a height based
  !! vertical system, too. This is due to a 'wrong' computation of z_delta
  !! in the case of a height based coordinate system (i.e. z_delta is
  !! implicity multiplied by -1)
  !!
  !! Literature
  !! Lin and Rood (1996), MWR, 124, 2046-2070
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-02-04)
  !!
  !! mpi parallelized, only cells_in_domain are computed, no sync
  SUBROUTINE v_ppm_slimiter_mo( p_patch_3D, p_cc, p_face, p_slope, p_face_up, p_face_low )

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    REAL(wp), INTENT(INOUT)           :: p_cc(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)      !< advected cell centered variable
    REAL(wp), INTENT(INOUT)           :: p_face(nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)  !< reconstructed face values of the advected field
    REAL(wp), INTENT(INOUT)           :: p_slope(nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) !< monotonized slope
    REAL(wp), INTENT(INOUT)           :: p_face_up(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) !< final face value (upper face, height based)
    REAL(wp), INTENT(INOUT)           :: p_face_low(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)!< final face value (lower face, height based)

    ! locals
    INTEGER  :: nlev                      !< number of full levels
    INTEGER  :: slev                      !< vertical start level
    INTEGER  :: jc, jk, jb                !< index of cell, vertical level and block
    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: ikp1                      !< vertical level plus one
    INTEGER  :: i_dolic
    REAL(wp) :: z_delta                   !< lower minus upper face value
    REAL(wp) :: z_a6i                     !< curvature of parabola
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    cells_in_domain => p_patch%cells%in_domain

    slev = 1
    nlev = n_zlev
! !$OMP PARALLEL
! !$OMP DO PRIVATE(jb,jk,jc,i_startidx_c,i_endidx_c,ikp1,z_delta,z_a6i)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c

        i_dolic=p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) !v_base%dolic_c(jc,jb)

        DO jk = slev, i_dolic
          ! index of bottom half level
          ikp1 = jk + 1

          IF ( p_patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            z_delta   = p_face(jc,ikp1,jb) - p_face(jc,jk,jb)
            z_a6i     = 6._wp * (p_cc(jc,jk,jb)                           &
              &       - 0.5_wp * (p_face(jc,jk,jb) + p_face(jc,ikp1,jb)))

            IF ( p_slope(jc,jk,jb) == 0._wp) THEN
              p_face_up(jc,jk,jb)  = p_cc(jc,jk,jb)
              p_face_low(jc,jk,jb) = p_cc(jc,jk,jb)

            ELSE IF (z_delta * z_a6i > z_delta * z_delta) THEN
              p_face_up(jc,jk,jb)  = 3._wp*p_cc(jc,jk,jb) - 2._wp*p_face(jc,ikp1,jb)
              p_face_low(jc,jk,jb) = p_face(jc,ikp1,jb)

            ELSE IF (z_delta * z_a6i < -1._wp * (z_delta * z_delta)) THEN
              p_face_up(jc,jk,jb)  = p_face(jc,jk,jb)
              p_face_low(jc,jk,jb) = 3._wp*p_cc(jc,jk,jb) - 2._wp*p_face(jc,jk,jb)

            ELSE
              p_face_up(jc,jk,jb)  = p_face(jc,jk,jb)
              p_face_low(jc,jk,jb) = p_face(jc,ikp1,jb)
            ENDIF
          ENDIF
        END DO
      END DO
    END DO
! !$OMP END DO
! !$OMP END PARALLEL

  END SUBROUTINE v_ppm_slimiter_mo
 !-------------------------------------------------------------------------
  !>
  !! Positive definite flux limiter for vertical advection
  !!
  !! Positive definite Zalesak Flux-Limiter (Flux corrected transport).
  !! for the hydrostatic core. Only outward fluxes are re-scaled, in 
  !! order to maintain positive definiteness.
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !! - Harris, L. M. and P. H. Lauritzen (2010): A flux-form version of the
  !!   Conservative Semi-Lagrangian Multi-tracer transport scheme (CSLAM) on
  !!   the cubed sphere grid.  J. Comput. Phys., 230, 1215-1237
  !! - Smolarkiewicz, P. K., 1989: Comment on "A positive definite advection 
  !!   scheme obtained by nonlinear renormalization of the advective fluxes.", 
  !!   Mon. Wea. Rev., 117, 2626-2632
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2011-01-07)
  !!
  SUBROUTINE vflx_limiter_pd_oce( p_patch_3D, p_dtime, p_cc, p_cellhgt_mc_now, p_flx_tracer_v, &
    &                            opt_slev, opt_elev )

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    REAL(wp), INTENT(INOUT)          :: p_cc(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)           !< advected cell centered variable at time (n)
    REAL(wp), INTENT(INOUT)          :: p_cellhgt_mc_now(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), INTENT(IN)             :: p_dtime
    REAL(wp), INTENT(INOUT)          :: p_flx_tracer_v(nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) !< calculated vertical tracer mass flux
    INTEGER,  INTENT(IN), OPTIONAL   :: opt_slev
    INTEGER,  INTENT(IN), OPTIONAL   :: opt_elev
!
!Local variables
    REAL(wp) :: r_m(nproma,n_zlev)      !< fraction which must multiply all
                                         !< outgoing fluxes of cell jc
                                         !< to guarantee positive definiteness
    REAL(wp) :: p_m(nproma)            !< sum of fluxes out of cell
    REAL(wp) :: z_signum(nproma)       !< sign of mass flux       !< >0: upward; <0: downward

    INTEGER  :: slev, elev             !< vertical start and end level
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: jk, jb, jc            !< index of edge, vert level, block, cell
    INTEGER  :: jkp1, jkm1
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    cells_in_domain => p_patch%cells%in_domain

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = n_zlev
    END IF

   r_m(1:nproma,1:n_zlev)= 0.0_wp
   p_m(1:nproma)         = 0.0_wp
   z_signum(1:nproma)    = 0.0_wp


! !$OMP PARALLEL
! !$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,jkp1,p_m,r_m,jkm1,z_signum)
   DO jb = cells_in_domain%start_block, cells_in_domain%end_block
     CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
     !
     ! 1. Compute total outward mass
     !
     DO jk = slev, elev-1
       jkp1 = jk+1

       DO jc = i_startidx, i_endidx

         IF ( p_patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
           ! Sum of all outgoing fluxes out of cell jk
           p_m(jc) = p_dtime                                  &
           &     * (MAX(0._wp,p_flx_tracer_v(jc,jkp1,jb))  &  ! upper half level
           &      - MIN(0._wp,p_flx_tracer_v(jc,jk,jb)) )     ! lower half level

           ! fraction which must multiply the fluxes out of cell jk to guarantee no
           ! undershoot
           ! Nominator: maximum allowable decrease \rho^n q^n
           r_m(jc,jk) = MIN(1._wp, (p_cc(jc,jk,jb)*p_cellhgt_mc_now(jc,jk,jb)) &
           &         /(p_m(jc) + dbl_eps) )
         ENDIF
       ENDDO

     ENDDO

 
     DO jc = i_startidx, i_endidx
 
       IF ( p_patch_3D%lsm_c(jc,elev,jb) <= sea_boundary ) THEN
         ! Sum of all outgoing fluxes out of cell jk
         p_m(jc) = p_dtime                                  &
         &     * (MAX(0._wp,p_flx_tracer_v(jc,elev+1,jb))   &  ! upper half level
         &      - MIN(0._wp,p_flx_tracer_v(jc,elev,jb)) )     ! lower half level

         ! fraction which must multiply the fluxes out of cell jk to guarantee no
         ! undershoot
         ! Nominator: maximum allowable decrease \rho^n q^n
         r_m(jc,elev) = MIN(1._wp, (p_cc(jc,elev,jb)*p_cellhgt_mc_now(jc,elev,jb)) &
         &         /(p_m(jc) + dbl_eps) )
       ENDIF
     ENDDO

     !
     ! 2. Limit outward fluxes (loop over half levels)
     !    Choose r_m depending on the sign of p_mflx_tracer_v
     !
     DO jk = slev+1, elev
     jkm1 = jk-1
       DO jc = i_startidx, i_endidx

         IF ( p_patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
           ! p_mflx_tracer_v(k-1/2) > 0: flux directed from cell k-1 -> k
           ! p_mflx_tracer_v(k-1/2) < 0: flux directed from cell k   -> k-1
           z_signum(jc) = SIGN(1._wp,p_flx_tracer_v(jc,jk,jb))

           p_flx_tracer_v(jc,jk,jb) =  p_flx_tracer_v(jc,jk,jb)  * 0.5_wp    &
             &                       * ( (1._wp + z_signum(jc)) * r_m(jc,jkm1) &
             &                       +   (1._wp - z_signum(jc)) * r_m(jc,jk) )
         ENDIF
       ENDDO
     ENDDO
   ENDDO
! !$OMP END DO NOWAIT
! !$OMP END PARALLEL


  END SUBROUTINE vflx_limiter_pd_oce

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Lax Friedrichs first order upwind flux for vertical advection,.
  !!
  !! Generalized Lax Friedrichs first order upwind flux,
  !! used in conservative vertical advection routines.
  !! For passive advection, equivalent to any other first
  !! order upwind flux.
  !! Applicable to both pressure based and height based vertical
  !! coordinate systems. Depending on the coordinate system chosen,
  !! the sign of the second term in the flux equation changes.
  !! - (-) for pressure based vertical coordinate systems
  !! - (+) for height based coordinate systems
  !! In order to get the correct sign, the variable p_coeff_grid
  !! has been introduced which is =1 for pressure based and =-1
  !! for height based coordinate systems.
  !!
  !! @par Revision History
  !! Developed  by L.Bonaventura  (2004).
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized for p- and z-based vertical coordinate systems
  !!
  FUNCTION laxfr_upflux_v( p_vn, p_psi1, p_psi2 )  RESULT(p_upflux)
    !

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: p_vn
    REAL(wp), INTENT(in) :: p_psi1, p_psi2

    REAL(wp) :: p_upflux

    !-----------------------------------------------------------------------
    !p_upflux = 0.5_wp * (                       p_vn  *( p_psi1 + p_psi2 )    &
    !  &                   - p_coeff_grid * ABS( p_vn )*( p_psi2 - p_psi1 ) )
    p_upflux = 0.5_wp * (                       p_vn  *( p_psi1 + p_psi2 )    &
      &                   +  ABS( p_vn )*( p_psi2 - p_psi1 ) )
  END FUNCTION laxfr_upflux_v


  !-------------------------------------------------------------------------
! !   !! SUBROUTINE advects vertically the tracers present in the ocean model.
! !   !!
! !   !! @par Revision History
! !   !! Developed  by  Peter Korn, MPI-M (2010).
! !   !!
! !   !! mpi parallelized, sync required: trac_out
! !   SUBROUTINE advect_diffuse_vertical(p_patch, trac_in, trac_old, &
! !                            & p_os,                           &
! !                            & bc_top_tracer, bc_bot_tracer,   &
! !                            & A_v,                            &
! !                            & trac_out, timestep, delta_t,    &
! !                            & cell_thick_intermed_c,          &
! !                            & FLUX_CALCULATION_VERT, tracer_id)
! ! 
! !     TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
! !     REAL(wp), INTENT(IN)              :: trac_in(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
! !     REAL(wp), INTENT(IN)              :: trac_old(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
! !     TYPE(t_hydro_ocean_state), TARGET :: p_os
! !     REAL(wp)                          :: bc_top_tracer(nproma, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
! !     REAL(wp)                          :: bc_bot_tracer(nproma, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
! !     REAL(wp)                          :: A_v(:,:,:)                                   !vertical mixing coeff
! !     REAL(wp), INTENT(inout)             :: trac_out(:,:,:)                              !new tracer
! !     INTEGER                           :: timestep
! !     REAL(wp)                          :: delta_t
! !     REAL(wp), INTENT(INOUT)           :: cell_thick_intermed_c(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
! !     INTEGER                           :: FLUX_CALCULATION_VERT
! !     INTEGER, INTENT(IN)               :: tracer_id
! ! 
! !     !Local variables
! !     REAL(wp) :: delta_z, delta_z2
! !     INTEGER  :: i_startidx_c, i_endidx_c
! !     INTEGER  :: jc, jk, jb!, je!jkp1        !< index of edge, vert level, block
! !     INTEGER  :: z_dolic
! !     REAL(wp) :: z_adv_flux_v (nproma, n_zlev+1, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)  ! vertical advective tracer flux
! !     REAL(wp) :: z_div_adv_v (nproma, n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)        ! vertical tracer divergence
! !     REAL(wp) :: z_div_diff_v (nproma, n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)        ! vertical tracer divergence
! !     REAL(wp) :: z_temp(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
! !     REAL(wp) :: z_diff_flux_v(nproma, n_zlev+1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)   ! vertical diffusive tracer flux
! ! 
! !     REAL(wp) :: z_h(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
! !     ! CHARACTER(len=max_char_length), PARAMETER :: &
! !     !        & routine = ('mo_tracer_advection:advect_individual_tracer')
! !     !-------------------------------------------------------------------------------
! !     TYPE(t_subset_range), POINTER :: cells_in_domain
! !     !-------------------------------------------------------------------------------
! !     cells_in_domain => p_patch%cells%in_domain
! ! 
! !     z_adv_flux_v   = 0.0_wp
! !     z_div_adv_v    = 0.0_wp
! !     z_div_diff_v   = 0.0_wp
! !     z_diff_flux_v  = 0.0_wp
! !     z_temp         = 0.0_wp
! !     z_h            = 1.0_wp  !  division by z_h in tracer_diffusion_vert_expl
! ! 
! !     jk = 1
! !     DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !       CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
! !       DO jc = i_startidx_c, i_endidx_c
! !         IF ( v_base%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
! !           delta_z = v_base%del_zlev_m(jk)+p_os%p_prog(nold(1))%h(jc,jb)&
! !                   &*v_base%wet_c(jc,jk,jb)
! !    
! !           p_os%p_diag%depth_c(jc,jk,jb) = delta_z
! !    
! !           cell_thick_intermed_c(jc,jk,jb)&
! !           & = delta_z-delta_t*p_os%p_diag%div_mass_flx_c(jc,jk,jb)&
! !           & -delta_t*(p_os%p_diag%w_time_weighted(jc,jk,jb)        &
! !           & - p_os%p_diag%w_time_weighted(jc,jk+1,jb))
! !         ENDIF
! !       END DO
! !     END DO
! !    
! !     DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !       CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
! !       DO jk = 2, n_zlev
! !         delta_z = v_base%del_zlev_m(jk)
! !         DO jc = i_startidx_c, i_endidx_c
! !           IF ( v_base%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
! !             p_os%p_diag%depth_c(jc,jk,jb) = delta_z
! !    
! !             cell_thick_intermed_c(jc,jk,jb)= &
! !             & delta_z-delta_t*p_os%p_diag%div_mass_flx_c(jc,jk,jb)&
! !             & -delta_t*(p_os%p_diag%w_time_weighted(jc,jk,jb)     &
! !             & - p_os%p_diag%w_time_weighted(jc,jk+1,jb))
! !           ENDIF
! !         END DO
! !       END DO
! !     END DO
! ! 
! !     !---------DEBUG DIAGNOSTICS-------------------------------------------
! !     idt_src=5  ! output print level (1-5, fix)
! !     CALL dbg_print('AdvVert: div_adv_v'        ,z_div_adv_v              ,str_module,idt_src)
! !     !---------------------------------------------------------------------
! ! 
! ! 
! !     ! Initialize timer for horizontal advection
! !     IF (ltimer) CALL timer_start(timer_adv_vert)
! ! 
! !     SELECT CASE(FLUX_CALCULATION_VERT)
! ! 
! !     CASE(UPWIND)
! ! 
! !       CALL upwind_vflux_oce( p_patch,                    &
! !                            & trac_old,                   &
! !                            & p_os%p_diag%w_time_weighted,& 
! !                            & bc_top_tracer,              &
! !                            & z_adv_flux_v,tracer_id )
! !     CASE(CENTRAL)
! !       CALL central_vflux_oce( p_patch,                   &
! !                            & trac_old,                   &
! !                            & p_os%p_diag%w_time_weighted,&
! !                            & z_adv_flux_v, tracer_id)
! !     CASE(MIMETIC,MIMETIC_MIURA)
! !       CALL upwind_vflux_ppm( p_patch, trac_old,          &
! !                            & p_os%p_diag%w_time_weighted,&
! !                            & dtime, 1 ,                  & !p_itype_vlimit,  &
! !                            & cell_thick_intermed_c,      &!p_cellhgt_mc_now, &
! !                            & z_adv_flux_v, tracer_id)
! !     END SELECT
! ! 
! !     IF (ltimer) CALL timer_stop(timer_adv_vert)
! ! 
! !     !divergence is calculated for advective fluxes
! !     DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !       CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
! !       DO jc = i_startidx_c, i_endidx_c
! !         z_dolic = v_base%dolic_c(jc,jb)
! !          !IF ( v_base%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
! !           DO jk = 1, z_dolic
! !             ! positive vertical divergence in direction of w (upward positive)
! !              z_div_adv_v(jc,jk,jb) = z_adv_flux_v(jc,jk,jb) &
! !                                    &- z_adv_flux_v(jc,jk+1,jb)
! !           END DO
! !         !ENDIF
! !       END DO
! !     END DO
! ! 
! !     !---------DEBUG DIAGNOSTICS-------------------------------------------
! !     idt_src=4  ! output print level (1-5, fix)
! !     CALL dbg_print('AdvVert: adv_flux_v'       ,z_adv_flux_v             ,str_module,idt_src)
! !     CALL dbg_print('AdvVert: div_adv_v'        ,z_div_adv_v              ,str_module,idt_src)
! !     !---------------------------------------------------------------------
! ! 
! !     IF (ltimer) CALL timer_start(timer_dif_vert)
! ! 
! !     !Case: Implicit Vertical diffusion
! !     IF(expl_vertical_tracer_diff==1)THEN
! ! 
! !       !Add advective part to old tracer
! !       !surface forcing applied as volume forcing at rhs, i.e.part of explicit term in momentum and tracer eqs.
! !       !in this case, top boundary ondition of vertical Laplacians are homogeneous
! !       jk=1
! !       DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !         CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
! !           DO jc = i_startidx_c, i_endidx_c
! !             z_dolic = v_base%dolic_c(jc,jb)
! !             !IF(z_dolic>=MIN_DOLIC)THEN
! !             IF ( v_base%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
! !                  delta_z = v_base%del_zlev_m(jk)+p_os%p_prog(nnew(1))%h(jc,jb)
! ! 
! !                z_temp(jc,jk,jb)= trac_in(jc,jk,jb)&
! !                & -(delta_t/delta_z)*(z_div_adv_v(jc,jk,jb)-bc_top_tracer(jc,jb))
! ! 
! !             ENDIF
! !         END DO
! !       END DO
! ! 
! !       DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !         CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
! !           DO jc = i_startidx_c, i_endidx_c
! !             z_dolic = v_base%dolic_c(jc,jb)
! !             !IF(z_dolic>=MIN_DOLIC)THEN
! !             DO jk = 2, z_dolic
! !               IF ( v_base%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
! !                 delta_z = v_base%del_zlev_m(jk)
! ! 
! !                 z_temp(jc,jk,jb)= trac_in(jc,jk,jb)&
! !                 & -(delta_t/delta_z)*z_div_adv_v(jc,jk,jb)
! ! 
! ! ! IF(v_base%lsm_c(jc,jk,jb)==-1)THEN
! ! !   write(204060,*)'data',jc,jk,v_base%dolic_c(jc,jb),jb,z_temp(jc,jk,jb),&
! ! ! !  &trac_in(jc,jk,jb),&
! ! !   &trac_in(jc,jk,jb)*cell_thick_intermed_c(jc,jk,jb)/delta_z,&
! ! !   &-delta_t*z_div_adv_v(jc,jk,jb)/delta_z,&
! ! ! !   &cell_thick_intermed_c(jc,jk,jb)/delta_z,&
! ! ! ! &v_base%dolic_c(jc,jb),&
! ! ! &z_adv_flux_v(jc,jk,jb), z_adv_flux_v(jc,jk+1,jb),&
! ! ! &p_os%p_diag%w_time_weighted(jc,jk,jb),&
! ! ! &p_os%p_diag%w_time_weighted(jc,jk+1,jb)
! ! ! Endif
! !               ENDIF
! !             ENDDO
! !         END DO
! !       END DO
! ! 
! !       !---------DEBUG DIAGNOSTICS-------------------------------------------
! !       idt_src=5  ! output print level (1-5, fix)
! !       CALL dbg_print('AdvVert: bef.impl.diff:temp',z_temp                   ,str_module,idt_src)
! !       CALL dbg_print('AdvVert: bef.impl.diff:flux',z_adv_flux_v             ,str_module,idt_src)
! !       !---------------------------------------------------------------------
! ! 
! !       !calculate vert diffusion impicit: result is stored in trac_out
! !       CALL tracer_diffusion_vert_impl_hom( p_patch,           &
! !                                      & z_temp,                &
! !                                      & p_os%p_prog(nnew(1))%h,&
! !                                      & A_v,                   &
! !                                      & trac_out(:,:,:))
! ! 
! !       !---------DEBUG DIAGNOSTICS-------------------------------------------
! !       idt_src=5  ! output print level (1-5, fix)
! !       CALL dbg_print('AdvVert: aft.impl.diff:trac',trac_out                 ,str_module,idt_src)
! !       !---------------------------------------------------------------------
! ! 
! !     !vertival diffusion is calculated explicitely
! !     ELSEIF(expl_vertical_tracer_diff==0)THEN
! ! 
! !       CALL finish("advect_diffuse_vertical", &
! !       &"xplicit vertical tracer mixing currently not supported -  change namlist option")
! ! 
! ! !      ! #slo# 2012-06-11: ATTENTION: z_h is used for division in explicit tracer diffusion
! ! !      !                              but z_h is not calculated!
! ! !      CALL tracer_diffusion_vert_expl( p_patch,       &
! ! !                                    & trac_in,        &! z_trac_c,      &
! ! !                                    &  z_h,           &
! ! !                                    &  bc_top_tracer, &
! ! !                                    &  bc_bot_tracer, & 
! ! !                                    &  A_v,           &
! ! !                                   &  z_div_diff_v)
! ! 
! !       ! #slo# 2012-06-11: ATTENTION: G_n/G_nimd/G_nm1_c_v not used any more
! ! 
! !   !   DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !   !     CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
! !   !     DO jc = i_startidx_c, i_endidx_c
! !   !       !interior: from one below surface to the ground
! !   !       z_dolic = v_base%dolic_c(jc,jb)
! !   !       IF(z_dolic>=MIN_DOLIC)THEN
! !   !         DO jk = 1, z_dolic
! !   !            delta_z  = v_base%del_zlev_m(jk)
! !   !            ! positive vertical divergence in direction of w (upward positive)
! !   !            G_n_c_v(jc,jk,jb) = z_div_adv_v(jc,jk,jb)/z_h(jc,jk,jb) - z_div_diff_v(jc,jk,jb)
! ! 
! !   !         END DO
! !   !       ENDIF
! !   !     END DO
! !   !   END DO
! ! 
! !   !   IF( is_initial_timestep(timestep))THEN
! !   !     G_nimd_c_v(:,:,:) = G_n_c_v(:,:,:)
! !   !   ELSE
! !   !     G_nimd_c_v(:,:,:) = (1.5_wp+AB_const)* G_n_c_v(:,:,:)   &
! !   !       &               - (0.5_wp+AB_const)*G_nm1_c_v(:,:,:)
! !   !   ENDIF
! ! 
! !       !Add advective and diffusive part to old tracer
! !   !   DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !   !     CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
! !   !     DO jc = i_startidx_c, i_endidx_c
! !   !       z_dolic = v_base%dolic_c(jc,jb)
! !   !       IF(z_dolic>=MIN_DOLIC)THEN
! !   !         !top  level
! !   !         jk=1
! !   !         delta_z = (p_os%p_prog(nold(1))%h(jc,jb)+v_base%del_zlev_m(jk))&
! !   !                 &/z_h(jc,jk,jb) 
! !   ! !         z_trac_c(jc,jk,jb) = (trac_in(jc,jk,jb)*delta_z      &
! !   ! !                            & -delta_t*z_div_adv_v(jc,jk,jb)) &
! !   ! !                            &/z_h(jc,jk,jb)
! !   ! !          trac_out(jc,jk,jb) = z_trac_c(jc,jk,jb)               &
! !   ! !                             & +delta_t*z_div_diff_v(jc,jk,jb)  
! !   ! !         trac_out(jc,jk,jb) = trac_in(jc,jk,jb)*delta_z      &
! !   ! !                            & -delta_t*(z_div_adv_v(jc,jk,jb) &
! !   ! !                            &/z_h(jc,jk,jb)                  &
! !   ! !                            &-z_div_diff_v(jc,jk,jb))
! !   !         trac_out(jc,jk,jb) = trac_in(jc,jk,jb)*delta_z      &
! !   !                            & -delta_t*G_nimd_c_v(jc,jk,jb) 
! ! 
! !   !         !interior down to bottom
! !   !         DO jk = 2, z_dolic
! !   !           delta_z = v_base%del_zlev_m(jk)/z_h(jc,jk,jb)
! !   ! !           z_trac_c(jc,jk,jb) = (trac_in(jc,jk,jb)*delta_z     &
! !   ! !                              & -delta_t*z_div_adv_v(jc,jk,jb))  &
! !   ! !                              &/z_h(jc,jk,jb) 
! !   ! !          trac_out(jc,jk,jb) = z_trac_c(jc,jk,jb)                &
! !   ! !                               & +delta_t*z_diff_flux_v(jc,jk,jb)
! !   ! !           trac_out(jc,jk,jb) = trac_in(jc,jk,jb)*delta_z      &
! !   ! !                              & -delta_t*(z_div_adv_v(jc,jk,jb) &
! !   ! !                              &/z_h(jc,jk,jb)                  &
! !   ! !                              &-z_div_diff_v(jc,jk,jb))
! !   !           trac_out(jc,jk,jb) = trac_in(jc,jk,jb)*delta_z      &
! !   !                              & -delta_t*G_nimd_c_v(jc,jk,jb) 
! !   !         END DO
! !   !       ELSE
! !   !         trac_out(jc,:,jb) = 0.0_wp
! !   !       ENDIF
! !   !     END DO
! !   !   END DO
! ! 
! !       !---------DEBUG DIAGNOSTICS-------------------------------------------
! !       idt_src=4  ! output print level (1-5, fix)
! !       CALL dbg_print('AdvVert:ExplDiff: trac_out',trac_out                 ,str_module,idt_src)
! !       idt_src=5  ! output print level (1-5, fix)
! !       CALL dbg_print('AdvVert:ExplDiff: diff_flx',z_diff_flux_v            ,str_module,idt_src)
! !       CALL dbg_print('AdvVert:ExplDiff: div_diff',z_div_diff_v             ,str_module,idt_src)
! !       !---------------------------------------------------------------------
! ! 
! !     ENDIF ! lvertical_diff_implicit
! ! 
! !     IF (ltimer) CALL timer_stop(timer_dif_vert)
! ! 
! !     CALL sync_patch_array(SYNC_C, p_patch, trac_out)
! ! 
! !   END SUBROUTINE advect_diffuse_vertical
END MODULE mo_oce_tracer_transport_vert
