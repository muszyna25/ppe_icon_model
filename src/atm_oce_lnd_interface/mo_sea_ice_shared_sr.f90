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
MODULE mo_sea_ice_shared_sr
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2007
  !
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_dynamics_config,     ONLY: nold
  USE mo_model_domain,        ONLY: t_patch
  USE mo_exception,           ONLY: finish, message
!  USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell, sea_boundary 
!  USE mo_loopindices,         ONLY: get_indices_c
!  USE mo_math_utilities,      ONLY: t_cartesian_coordinates
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref,ki,ks,Tf,albi,albim,albsm,albs,&
    &                               mu,mus,ci, alf, I_0, alv, albedoW, clw,            &
    &                               cpd, zemiss_def,rd, stbo,tmelt   
  USE mo_math_constants,      ONLY: rad2deg
!  USE mo_ocean_nml,           ONLY: no_tracer, init_oce_prog, iforc_oce, &
!    &                               FORCING_FROM_FILE_FLUX
  USE mo_oce_state,           ONLY: t_hydro_ocean_state, v_base, ocean_var_list
!  USE mo_oce_index,           ONLY: print_mxmn, ipl_src
!  USE mo_var_list,            ONLY: add_var
!  USE mo_master_control,      ONLY: is_restart_run
!  USE mo_cf_convention
!  USE mo_grib2
!  USE mo_cdi_constants
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes,&
    &                               t_atmos_for_ocean
  USE mo_io_units,            ONLY: nerr
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range 
  USE mo_loopindices,         ONLY: get_indices_c


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: oce_ice_heatflx
  PUBLIC :: print_maxmin_si
  PUBLIC :: print_cells

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
 
  
  SUBROUTINE oce_ice_heatflx (p_os,ice,Tfw,zHeatOceI)
    TYPE(t_hydro_ocean_state), INTENT(IN) :: p_os
    TYPE(t_sea_ice),           INTENT(IN) :: ice
    REAL(wp),                  INTENT(IN):: Tfw(:,:,:) ! freezing temperature
    REAL(wp),                  INTENT(OUT):: zHeatOceI(:,:,:)

    ! Local
    INTEGER :: k ! counter for ice thickness categories

    
    ! calculate heat flux from ocean to ice  (zHeatOceI) 
    DO k=1,ice%kice
      WHERE (ice%isice(:,k,:)) 
        zHeatOceI(:,k,:) = ( p_os%p_prog(nold(1))%tracer(:,1,:,1) - Tfw(:,k,:) ) &
          &                 * ice%zUnderIce(:,:) * clw*rho_ref/dtime
      ENDWHERE
    END DO
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
    DO jb = 1,p_patch%nblks_c
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c) 
      DO jc = i_startidx_c,i_endidx_c
!        IF (ice%isice(jc,1,jb)) THEN
        DO k=1,ice%kice
          IF (ice%isice(jc,k,jb)) THEN
            ctr = ctr+1
            values(ctr) = A(jc,jb)
          END IF
        END DO
      END DO
    END DO
    
991 FORMAT(a10,a25,':     ', 2g26.18) 
983 FORMAT(a10,a25,':     ', 4i4)
    WRITE(nerr,991) ' MAX/MIN ',description, &
      &              MAXVAL(A(:,:),MASK=ice%isice(:,1,:)),     &
      &              MINVAL(A(:,:),MASK=ice%isice(:,1,:))
    WRITE(nerr,983) ' LOC ',description, &
      &              MAXLOC(A(:,:),MASK=ice%isice(:,1,:)),     &
      &              MINLOC(A(:,:),MASK=ice%isice(:,1,:))
  

  END SUBROUTINE print_maxmin_si


  !--- check output in a single grid cell
  
  SUBROUTINE print_cells (A,ice,p_patch,description)
    REAL(wp), INTENT(IN)         :: A(:,:)
    TYPE(t_sea_ice), INTENT(IN)  :: ice
    TYPE(t_patch), INTENT(IN)    :: p_patch
    CHARACTER(len=*), INTENT(IN) :: description

    !local
    !REAL(wp),ALLOCATABLE :: values(:)
    INTEGER :: ctr, jb, k, jc, i_startidx_c, i_endidx_c



991 FORMAT(a15,a25,':  ', 2g26.18) 
   WRITE(nerr,991) ' at single cell(s) ',description, &
     &            A(7,17), &
     &            A(3,100)
  END SUBROUTINE print_cells


END MODULE mo_sea_ice_shared_sr
