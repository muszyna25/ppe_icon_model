!>
!! @brief Interface between ECHAM physics and the ocean, through a coupler
!!
!! @author Marco Giorgetta (MPI-M)
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

MODULE mo_interface_echam_ocean

  USE mo_kind                ,ONLY: wp
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_echam_phy_memory    ,ONLY: prm_field
                                
  USE mo_parallel_config     ,ONLY: nproma
  
  USE mo_run_config          ,ONLY: nlev
  USE mo_icoham_sfc_indices  ,ONLY: iwtr, iice

  USE mo_sync                ,ONLY: SYNC_C, sync_patch_array
#ifdef YAC_coupling
  USE finterface_description ,ONLY: yac_fput, yac_fget, yac_fget_nbr_fields, yac_fget_field_ids
#else
  USE mo_icon_cpl_exchg      ,ONLY: ICON_cpl_put, ICON_cpl_get
  USE mo_icon_cpl_def_field  ,ONLY: ICON_cpl_get_nbr_fields, ICON_cpl_get_field_ids
  USE mo_icon_cpl_restart    ,ONLY: icon_cpl_write_restart
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_ocean
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_interface_echam_ocean'

CONTAINS
  !>
  !! SUBROUTINE interface_icoham_echam -- the interface between
  !! ECHAM physics and the ocean, through a coupler
  !!
  !! This subroutine is called in the time loop of the ICONAM model.
  !! It takes the following as input:
  !! <ol>
  !! <li> prognostic and diagnostic variables of the dynamical core;
  !! <li> tendency of the prognostic varibles induced by adiabatic dynamics;
  !! <li> time step;
  !! <li> information about the dynamics grid;
  !! <li> interplation coefficients.
  !! </ol>
  !!
  !! The output includes tendencies of the prognostic variables caused by
  !! the parameterisations.
  !!
  !! Note that each call of this subroutine deals with a single grid level
  !! rather than the entire grid tree.

  SUBROUTINE interface_echam_ocean( jg      ,&! in
    &                               p_patch ) ! in

    ! Arguments

    INTEGER,               INTENT(IN)    :: jg            !< grid level/domain index
    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch

    ! Local variables

    LOGICAL               :: write_coupler_restart
    INTEGER               :: nbr_fields
    INTEGER               :: nbr_hor_points ! = inner and halo points
    INTEGER               :: nbr_points     ! = nproma * nblks
    INTEGER               :: field_shape(3)
    INTEGER , ALLOCATABLE :: field_id(:)
    REAL(wp), ALLOCATABLE :: buffer(:,:)

    INTEGER               :: info, ierror !< return values form cpl_put/get calls

!!$    CHARACTER(*), PARAMETER :: method_name = "interface_echam_ocean"

    !-------------------------------------------------------------------------
    ! If running in atm-oce coupled mode, exchange information 
    !-------------------------------------------------------------------------

    ! Possible fields that contain information to be sent to the ocean include
    !
    ! 1. prm_field(jg)% u_stress_tile(:,:,iwtr)  and 
    !    prm_field(jg)% v_stress_tile(:,:,iwtr)  which are the wind stress components;
    !
    ! 2. prm_field(jg)% evap_tile(:,:,iwtr) evaporation rate
    !
    ! 3. prm_field(jg)%rsfl + prm_field(jg)%rsfc + prm_field(jg)%ssfl + prm_field(jg)%ssfc
    !    which gives the precipitation rate;
    !
    ! 4. prm_field(jg)% temp(:,nlev,:)  temperature at the lowest model level, or
    !    prm_field(jg)% temp_2m(:,:)    2-m temperature, not available yet, or
    !    prm_field(jg)% shflx_tile(:,:,iwtr) sensible heat flux
    !
    ! 5  prm_field(jg)% lhflx_tile(:,:,iwtr) latent heat flux
    ! 6. shortwave radiation flux at the surface
    !
    ! Possible fields to receive from the ocean include
    !
    ! 1. prm_field(jg)% tsfc_tile(:,:,iwtr)   SST
    ! 2. prm_field(jg)% ocu(:,:) and ocv(:,:) ocean surface current
    ! 

    nbr_hor_points = p_patch%n_patch_cells
    nbr_points     = nproma * p_patch%nblks_c

    ALLOCATE(buffer(nproma*p_patch%nblks_c,5))
    buffer(:,:) = 0.0_wp

    !
    !  see drivers/mo_atmo_model.f90:
    !
    !   field_id(1) represents "TAUX"   wind stress component
    !   field_id(2) represents "TAUY"   wind stress component
    !   field_id(3) represents "SFWFLX" surface fresh water flux
    !   field_id(4) represents "SFTEMP" surface temperature
    !   field_id(5) represents "THFLX"  total heat flux
    !   field_id(6) represents "ICEATM" ice temperatures and melt potential
    !
    !   field_id(7) represents "SST"    sea surface temperature
    !   field_id(9) represents "OCEANU" u component of ocean surface current
    !   field_id(9) represents "OCEANV" v component of ocean surface current
    !   field_id(10)represents "ICEOCE" ice thickness, concentration and temperatures
    !
    !
#ifdef YAC_Coupling
    CALL yac_fget_nbr_fields ( nbr_fields )
    ALLOCATE(field_id(nbr_fields))
    CALL yac_fget_field_ids ( nbr_fields, field_id )
#else
    CALL ICON_cpl_get_nbr_fields ( nbr_fields )
    ALLOCATE(field_id(nbr_fields))
    CALL ICON_cpl_get_field_ids ( nbr_fields, field_id )
#endif
    !
    !
    field_shape(1) = 1
    field_shape(2) = nbr_hor_points
    field_shape(3) = 1

    !
    ! Send fields away
    ! ----------------
    !
    write_coupler_restart = .FALSE.
    !
    ! TAUX
    !
    buffer(:,:) = 0.0_wp
    buffer(:,1) = RESHAPE ( prm_field(jg)%u_stress_tile(:,:,iwtr), (/ nbr_points /) )
    buffer(:,2) = RESHAPE ( prm_field(jg)%u_stress_tile(:,:,iice), (/ nbr_points /) )

#ifdef YAC_coupling
    CALL yac_fput ( field_id(1), nbr_hor_points, 1, 1, 1, buffer, ierror )
#else
    field_shape(3) = 2
    CALL ICON_cpl_put ( field_id(1), field_shape, buffer(1:nbr_hor_points,1:2), info, ierror )
#endif
    IF ( info == 2 ) write_coupler_restart = .TRUE.
    !
    ! TAUY
    !
    buffer(:,1) = RESHAPE ( prm_field(jg)%v_stress_tile(:,:,iwtr), (/ nbr_points /) )
    buffer(:,2) = RESHAPE ( prm_field(jg)%v_stress_tile(:,:,iice), (/ nbr_points /) )
#ifdef YAC_coupling
    CALL yac_fput ( field_id(2), nbr_hor_points, 1, 1, 1, buffer, ierror )
#else
    CALL ICON_cpl_put ( field_id(2), field_shape, buffer(1:nbr_hor_points,1:2), info, ierror )
#endif
    IF ( info == 2 ) write_coupler_restart = .TRUE.
    !
    ! SFWFLX Note: the evap_tile should be properly updated and added
    !
    !       write(0,*)  prm_field(jg)%rsfl(:,:)
    !       write(0,*)  prm_field(jg)%rsfc(:,:)
    !       write(0,*)  prm_field(jg)%ssfl(:,:)
    !       write(0,*)  prm_field(jg)%ssfc(:,:)
    !       write(0,*)  prm_field(jg)%evap_tile(:,:,iwtr)

    buffer(:,1) = RESHAPE ( prm_field(jg)%rsfl(:,:), (/ nbr_points /) ) + &
      &        RESHAPE ( prm_field(jg)%rsfc(:,:), (/ nbr_points /) ) ! total rain
    buffer(:,2) = RESHAPE ( prm_field(jg)%ssfl(:,:), (/ nbr_points /) ) + &
      &        RESHAPE ( prm_field(jg)%ssfc(:,:), (/ nbr_points /) ) ! total snow
    buffer(:,3) = RESHAPE ( prm_field(jg)%evap_tile(:,:,iwtr), (/ nbr_points /) )

#ifdef YAC_coupling
    CALL yac_fput ( field_id(3), nbr_hor_points, 3, 1, 1, buffer, ierror )
#else
    field_shape(3) = 3
    CALL ICON_cpl_put ( field_id(3), field_shape, buffer(1:nbr_hor_points,1:3), info, ierror )
#endif
    IF ( info == 2 ) write_coupler_restart = .TRUE.
    !
    ! SFTEMP
    !
    buffer(:,1) =  RESHAPE ( prm_field(jg)%temp(:,nlev,:), (/ nbr_points /) )
#ifdef YAC_coupling
    CALL yac_fput ( field_id(4), nbr_hor_points, 1, 1, 1, buffer, ierror )
#else
    field_shape(3) = 1
    CALL ICON_cpl_put ( field_id(4), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
    IF ( info == 2 ) write_coupler_restart = .TRUE.
    !
    ! THFLX, total heat flux
    !
    buffer(:,1) =  RESHAPE ( prm_field(jg)%swflxsfc_tile(:,:,iwtr), (/ nbr_points /) ) !net shortwave flux for ocean
    buffer(:,2) =  RESHAPE ( prm_field(jg)%lwflxsfc_tile(:,:,iwtr), (/ nbr_points /) ) !net longwave flux
    buffer(:,3) =  RESHAPE ( prm_field(jg)%shflx_tile(:,:,iwtr),    (/ nbr_points /) ) !sensible heat flux
    buffer(:,4) =  RESHAPE ( prm_field(jg)%lhflx_tile(:,:,iwtr),    (/ nbr_points /) ) !latent heat flux for ocean
#ifdef YAC_coupling
    CALL yac_fput ( field_id(5), nbr_hor_points, 4, 1, 1, buffer, ierror )
#else
    field_shape(3) = 4
    CALL ICON_cpl_put ( field_id(5), field_shape, buffer(1:nbr_hor_points,1:4), info, ierror )
#endif
    !
    ! ICEATM, Ice state determined by atmosphere
    !
    buffer(:,1) =  RESHAPE ( prm_field(jg)%Qtop(:,1,:), (/ nbr_points /) ) !Melt-potential for ice - top
    buffer(:,2) =  RESHAPE ( prm_field(jg)%Qbot(:,1,:), (/ nbr_points /) ) !Melt-potential for ice - bottom
    buffer(:,3) =  RESHAPE ( prm_field(jg)%T1  (:,1,:), (/ nbr_points /) ) !Temperature of upper ice layer
    buffer(:,4) =  RESHAPE ( prm_field(jg)%T2  (:,1,:), (/ nbr_points /) ) !Temperature of lower ice layer
#ifdef YAC_coupling
    CALL yac_fput ( field_id(6), nbr_hor_points, 4, 1, 1, buffer, ierror )
#else
    field_shape(3) = 4
    CALL ICON_cpl_put ( field_id(6), field_shape, buffer(1:nbr_hor_points,1:4), info, ierror )
#endif
    IF ( info == 2 ) write_coupler_restart = .TRUE.
#ifdef YAC_coupling
    TODO
#else
    IF ( write_coupler_restart ) CALL icon_cpl_write_restart ( 6, field_id(1:6), ierror )
#endif
    !
    ! Receive fields, only assign values if something was received ( info > 0 )
    ! -------------------------------------------------------------------------
    !
    ! SST
    !
#ifdef YAC_coupling
    CALL yac_fget ( field_id(7), nbr_hor_points, 1, 1, 1, buffer, info, ierror )
#else
    field_shape(3) = 1
    CALL ICON_cpl_get ( field_id(7), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
    IF ( info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
      prm_field(jg)%tsfc_tile(:,:,iwtr) = RESHAPE (buffer(:,1), (/ nproma, p_patch%nblks_c /) )
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%tsfc_tile(:,:,iwtr))
    END IF
    !
    ! OCEANU
    !
#ifdef YAC_coupling
    CALL yac_fget ( field_id(8), nbr_hor_points, 1, 1, 1, buffer, info, ierror )
#else
    CALL ICON_cpl_get ( field_id(8), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
    IF ( info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
      prm_field(jg)%ocu(:,:) = RESHAPE (buffer(:,1), (/ nproma,  p_patch%nblks_c /) )
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocu(:,:))
    END IF
    !
    ! OCEANV
    !
#ifdef YAC_coupling
    CALL yac_fget ( field_id(9), nbr_hor_points, 1, 1, 1, buffer, info, ierror )
#else
    CALL ICON_cpl_get ( field_id(9), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
    IF ( info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
      prm_field(jg)%ocv(:,:) = RESHAPE (buffer(:,1), (/ nproma, p_patch%nblks_c /) )
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocv(:,:))
    END IF
    !
    ! ICEOCE
    !
#ifdef YAC_coupling
    CALL yac_fget ( field_id(7), nbr_hor_points, 4, 1, 1, buffer, info, ierror )
#else
    field_shape(3) = 5
    CALL ICON_cpl_get ( field_id(10), field_shape, buffer(1:nbr_hor_points,1:5), info, ierror )
#endif
    IF ( info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:5) = 0.0_wp
      prm_field(jg)%hi  (:,1,:) = RESHAPE (buffer(:,1), (/ nproma, p_patch%nblks_c /) )
      prm_field(jg)%hs  (:,1,:) = RESHAPE (buffer(:,2), (/ nproma, p_patch%nblks_c /) )
      prm_field(jg)%conc(:,1,:) = RESHAPE (buffer(:,3), (/ nproma, p_patch%nblks_c /) )
      prm_field(jg)%T1  (:,1,:) = RESHAPE (buffer(:,4), (/ nproma, p_patch%nblks_c /) )
      prm_field(jg)%T2  (:,1,:) = RESHAPE (buffer(:,5), (/ nproma, p_patch%nblks_c /) )
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%hi  (:,1,:))
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%hs  (:,1,:))
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%seaice(:,:))
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%T1  (:,1,:))
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%T2  (:,1,:))
      prm_field(jg)%seaice(:,:) = prm_field(jg)%conc(:,1,:)
    END IF

    DEALLOCATE(buffer)
    DEALLOCATE(field_id)

  END SUBROUTINE interface_echam_ocean

END MODULE mo_interface_echam_ocean
