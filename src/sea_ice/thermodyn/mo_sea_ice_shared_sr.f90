!>
!! Provide an implementation of the sea-ice model.
!!
!! Provide an implementation of the parameters of the surface module (sea ice)
!! used between the atmopshere and the hydrostatic ocean model.
!!
!! @author 
!! 
!! @par Revision History
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
MODULE mo_sea_ice_shared_sr
  
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_dynamics_config,     ONLY: nold
  USE mo_model_domain,        ONLY: t_patch
  USE mo_physical_constants,  ONLY: rho_ref, clw, Cd_io, Ch_io
  USE mo_ocean_types,           ONLY: t_hydro_ocean_state
  USE mo_sea_ice_types,       ONLY: t_sea_ice
  USE mo_io_units,            ONLY: nerr
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range 
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_sea_ice_nml,         ONLY: i_Qio_type
  USE mo_util_dbg_prnt,       ONLY: dbg_print


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: oce_ice_heatflx
  PUBLIC :: print_maxmin_si
  PUBLIC :: print_cells

  CHARACTER(len=12)           :: str_module    = 'SeaIceShared'  ! Output of module for 1 line debug

CONTAINS

  !-------------------------------------
  !
  ! oce_ice_heatflx
  !
  ! Calculates the heat flux from the uppermost water layer into the ice.
  !
  ! Currently (as in growth.f90): all energy available in upper ocean grid cell 
  ! is supplied to the ice and the upper ocean temperature is held at the 
  ! freezing point. This is not very physical.
  !
  ! Positive flux upwards.
 
  
  SUBROUTINE oce_ice_heatflx (p_patch, p_os,ice,Tfw,zHeatOceI)
    TYPE(t_patch)            , INTENT(IN), TARGET    :: p_patch
    TYPE(t_hydro_ocean_state), INTENT(IN)  :: p_os
    TYPE(t_sea_ice)          , INTENT(IN)  :: ice
    REAL(wp)                 , INTENT(IN)  :: Tfw(:,:,:)      ! freezing temperature
    REAL(wp)                 , INTENT(OUT) :: zHeatOceI(:,:,:)

    ! Local
    INTEGER :: jb, k, jc, i_startidx_c, i_endidx_c
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: u_star

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice_shared_sr:oce_ice_heatflx'
    
    all_cells => p_patch%cells%all 
    zHeatOceI = 0.0_wp

    ! calculate heat flux from ocean to ice  (zHeatOceI) 
    SELECT CASE ( i_Qio_type )
    CASE (1)
!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, k, jc) SCHEDULE(dynamic)
      DO jb = 1,p_patch%nblks_c
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c,i_endidx_c
          DO k=1,ice%kice
            IF (ice%hi(jc,k,jb) > 0._wp) THEN
              ! ALL energy of warm water over the whole grid area is used for melting ice - divide by concentration
              ! SLO/EO 2014-04-11 - this is wrong, must be accompanied by correction elsewhere,
              !                     since open part of water is still losing heat
              !                   - old formulation (without concentration, below) is reactivated
              ! zHeatOceI(jc,k,jb) = ( p_os%p_prog(nold(1))%tracer(jc,1,jb,1) - Tfw(jc,k,jb) )  &
              !   &                 * ice%zUnderIce(jc,jb) * clw*rho_ref/(dtime*ice%conc(jc,k,jb))
              ! Old (2014-02) formulation: part of warm water below ice covered area used only
                zHeatOceI(jc,k,jb) = ( p_os%p_prog(nold(1))%tracer(jc,1,jb,1) - Tfw(jc,k,jb) )  &
                    &                 * ice%zUnderIce(jc,jb) * clw*rho_ref/dtime
            ENDIF
          ENDDO
        ENDDO
      END DO
!ICON_OMP_END_PARALLEL_DO

   CASE(2)

!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, k, jc) SCHEDULE(dynamic)
    DO jb = 1,p_patch%nblks_c
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c,i_endidx_c
        DO k=1,ice%kice
          IF (ice%hi(jc,k,jb) > 0._wp) THEN
            ! energy of warm water over the ice covered part of grid area is used for melting ice - no division
              u_star = SQRT(Cd_io*( (p_os%p_diag%u(jc,1,jb)-ice%u(jc,jb))**2 + &
                &         (p_os%p_diag%v(jc,1,jb)-ice%v(jc,jb))**2 ))
              zHeatOceI(jc,k,jb) = ( p_os%p_prog(nold(1))%tracer(jc,1,jb,1) - Tfw(jc,k,jb) )  &
                &                         *rho_ref*clw*Ch_io*u_star
          ENDIF
        ENDDO
      ENDDO
    END DO
!ICON_OMP_END_PARALLEL_DO

    CASE (3)
!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, k, jc) SCHEDULE(dynamic)
      DO jb = 1,p_patch%nblks_c
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c,i_endidx_c
          DO k=1,ice%kice
            IF (ice%hi(jc,k,jb) > 0._wp) THEN
              ! ALL energy of warm water over the whole grid area is used for melting ice - divide by concentration
              ! SLO/EO 2014-04-11 - this is wrong, must be accompanied by correction elsewhere,
              !                     since open part of water is still losing heat
              !                   - old formulation (without concentration, below) is reactivated
              zHeatOceI(jc,k,jb) = ( p_os%p_prog(nold(1))%tracer(jc,1,jb,1) - Tfw(jc,k,jb) )  &
                &                 * ice%zUnderIce(jc,jb) * clw*rho_ref/(dtime*ice%conc(jc,k,jb))
            ENDIF
          ENDDO
        ENDDO
      END DO
!ICON_OMP_END_PARALLEL_DO

    
    CASE DEFAULT
      CALL finish(TRIM(routine), 'Invalid i_Qio_type')
    END SELECT
    
    CALL dbg_print('o-i-heat: Tfw'       ,Tfw          ,str_module,4, in_subset=p_patch%cells%owned)
    CALL dbg_print('o-i-heat: zUnderIce' ,ice%zUnderIce,str_module,4, in_subset=p_patch%cells%owned)
    CALL dbg_print('o-i-heat: zHeatOceI' ,zHeatOceI    ,str_module,3, in_subset=p_patch%cells%owned)
    
  END SUBROUTINE oce_ice_heatflx




!!$ debugging stuff and misc.

  !---
  !
  ! Do some stuff to check some output. This is not even close to a final version
  ! 

  SUBROUTINE print_maxmin_si(A,ice,p_patch,description)
    REAL(wp), INTENT(IN)         :: A(:,:)
    TYPE(t_sea_ice), INTENT(IN)  :: ice
    TYPE(t_patch), INTENT(IN), TARGET :: p_patch
    CHARACTER(len=*), INTENT(IN) :: description
    
    !local
    REAL(wp),ALLOCATABLE :: values(:)
    INTEGER :: ctr, jb, k, jc, i_startidx_c, i_endidx_c
    
    TYPE(t_subset_range), POINTER :: all_cells
    
    all_cells => p_patch%cells%all
    ALLOCATE(values(p_patch%nblks_c * nproma))
    
    ctr = 0
    values(:) = 1.0_wp
    DO jb = 1,p_patch%nblks_c
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c) 
      DO jc = i_startidx_c,i_endidx_c
        DO k=1,ice%kice
          IF (ice%hi(jc,k,jb) > 0._wp) THEN
            ctr = ctr+1
            values(ctr) = A(jc,jb)
          END IF
        END DO
      END DO
    END DO
    
101 FORMAT(a10,a25,':     ', 3g26.18) 
102 FORMAT(a10,a25,':     ', 4i4)
    
    WRITE(nerr,101) ' MAX/MIN ',description, &
      &              MAXVAL(A(:,:),MASK=(ice%hi(:,1,:)>0._wp)),     &
      &              MINVAL(A(:,:),MASK=(ice%hi(:,1,:)>0._wp)),     &
!!$      &              MAXVAL(values(1:ctr-1)), &
!!$      &              MINVAL(values(1:ctr-1)), &
      &              SUM(values(1:ctr-1))/SIZE(values(1:ctr-1))
    
    WRITE(nerr,102) ' LOC ',description, &
      &              MAXLOC(A(:,:),MASK=(ice%hi(:,1,:)>0._wp)),     &
      &              MINLOC(A(:,:),MASK=(ice%hi(:,1,:)>0._wp))
!!$      &              MAXLOC(values(1:ctr-1)), &
!!$      &              MINLOC(values(1:ctr-1)) 
      
  END SUBROUTINE print_maxmin_si
    

  !--- check output in a single grid cell
  
  SUBROUTINE print_cells (A,description)
    REAL(wp), INTENT(IN)         :: A(:,:)
!    TYPE(t_sea_ice), INTENT(IN)  :: ice
!    TYPE(t_patch), INTENT(IN)    :: p_patch
    CHARACTER(len=*), INTENT(IN) :: description

    !local
    !REAL(wp),ALLOCATABLE :: values(:)
!    INTEGER :: ctr, jb, k, jc, i_startidx_c, i_endidx_c



103 FORMAT(a18,a25,':     ', 3g26.18)
   WRITE(nerr,103) ' at single cell(s) ',description, &
     &            A(26,295), &
     &            A(21,289), &
     &            A(23,38)
  END SUBROUTINE print_cells


END MODULE mo_sea_ice_shared_sr
