!>
!! This module is the interface between nwp_nh_interface to the radiation schemes
!! (RRTM or Ritter-Geleyn).
!!
!! @author Leonidas Linardakis, MPI-M, 2012-08-21
!!
!! @par Revision History
!! Initial release by Leonidas Linardakis, MPI-M, 2012-08-21
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_rrtm_data_interface

  USE mo_exception,           ONLY: message,  finish !message_tex
  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch, p_patch_local_parent
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_lnd_diag
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_nonhydro_types,      ONLY: t_nh_diag
  
  USE mo_parallel_config,     ONLY: parallel_radiation_mode
  USE mo_mpi,                 ONLY: my_process_is_mpi_seq
  
  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER:: version = '$Id$'
    
    
  !--------------------------------------------------------------------------
  TYPE t_rrtm_data

    ! pointers to grids/states from input data
    ! not used currently
    TYPE(t_patch),        POINTER :: pt_patch     !<grid/patch info.
    TYPE(t_external_data),POINTER :: ext_data
    TYPE(t_lnd_diag),     POINTER :: lnd_diag     !< diag vars for sfc
    TYPE(t_nh_diag),      POINTER :: pt_diag      !<the diagnostic variables
    TYPE(t_nwp_phy_diag), POINTER :: prm_diag
    TYPE(t_lnd_prog),     POINTER :: lnd_prog_now

    !------------------------------------------
    ! grid parameters
    INTEGER   :: no_of_cells
    INTEGER   :: full_levels, half_levels
    ! blocking parameters
    INTEGER   :: nproma
    INTEGER   :: no_of_blocks, end_index

    ! communication patterns
    INTEGER :: send_out_data, recv_input_data

    !------------------------------------------
    ! INPUT
    REAL(wp)              :: p_sim_time
    
    INTEGER, POINTER  :: convection_type(:,:)     ! always zero
    
    REAL(wp), POINTER ::  fr_land_smt(:,:)   !< fraction land in a grid element        [ ]
                                                         !  = smoothed fr_land
    REAL(wp), POINTER ::  fr_glac_smt(:,:)   !< fraction land glacier in a grid element [ ]
                                                         ! = smoothed fr_glac
    REAL(wp), POINTER ::  cosmu0(:,:)        ! cosine of solar zenith angle
        
    REAL(wp), POINTER ::  albedo_vis_dir(:,:) ! surface albedo for visible and near IR range, direct
    REAL(wp), POINTER ::  albedo_vis_dif(:,:) !< in surface albedo for visible range, diffuse
    REAL(wp), POINTER ::  albedo_nir_dif(:,:) !< in surface albedo for visible range, diffuse
    
    REAL(wp), POINTER ::  emis_rad(:,:)      ! lw sfc emissivity
    REAL(wp), POINTER ::  tsfctrad(:,:)      ! surface temperature at trad [K]
    
    REAL(wp), POINTER ::  pres_ifc(:,:,:)    ! pressure at interfaces (nproma,nlevp1,nblks_c)  [Pa]
    REAL(wp), POINTER ::  pres(:,:,:)        ! pressure (nproma,nlev,nblks_c)                  [Pa]
    REAL(wp), POINTER ::  temp(:,:,:)        ! temperature (nproma,nlev,nblks_c)                 [K]
    
    REAL(wp), POINTER ::  qm_vapor(:,:,:)    !< water vapor mass mix ratio at t-dt
    REAL(wp), POINTER ::  qm_liquid(:,:,:)   !< cloud water mass mix ratio at t-dt
    REAL(wp), POINTER ::  qm_ice(:,:,:)      !< cloud ice mass mixing ratio at t-dt
    
    REAL(wp), POINTER ::  qm_o3(:,:,:)       ! in o3 mass mixing ratio
    REAL(wp), POINTER ::  acdnc(:,:,:)       ! cloud droplet number concentration [1/m**3]
    REAL(wp), POINTER ::  cld_frc(:,:,:)     !< cloud fraction [m2/m2]
    
    REAL(wp), POINTER ::  zaeq1(:,:,:)       !< in aerosol continental
    REAL(wp), POINTER ::  zaeq2(:,:,:)       !< in aerosol maritime
    REAL(wp), POINTER ::  zaeq3(:,:,:)       !< in aerosol urban
    REAL(wp), POINTER ::  zaeq4(:,:,:)       !< in aerosol volcano ashes
    REAL(wp), POINTER ::  zaeq5(:,:,:)       !< in aerosol stratospheric background

    !------------------------------------------
    ! OUTPUT    
    REAL(wp), POINTER ::  aclcov (:,:)       ! cloud cover in a column [m2/m2], not used at the moment
    REAL(wp), POINTER ::  lwflxclr(:,:,:)    ! terrestrial flux, clear sky, net down, longwave clear-sky net flux [W/m2]
    REAL(wp), POINTER ::  trsolclr(:,:,:)    ! sol. transmissivity, clear sky, net down, shortwave clear-sky net tranmissivity []
    REAL(wp), POINTER ::  lwflxall(:,:,:)    ! terrestrial flux, all sky, net down, longwave net flux           [W/m2]
    REAL(wp), POINTER ::  trsolall(:,:,:)    ! solar transmissivity, all sky, net dow, shortwave net tranmissivity []
    
  END TYPE t_rrtm_data
  !--------------------------------------------------------------------------

  TYPE(t_rrtm_data) :: rrtm_model_data
    

CONTAINS



  !-----------------------------------------
  !>
  SUBROUTINE init_rrtm_model_data(no_of_cells, full_levels, half_levels, nproma)
    INTEGER, INTENT(in) :: no_of_cells, full_levels, half_levels, nproma

    IF (parallel_radiation_mode <= 0) RETURN
    IF (my_process_is_mpi_seq()) THEN
      parallel_radiation_mode = 0
      RETURN
    ENDIF
    
    ! fill the rrtm_data parameters
    rrtm_model_data%no_of_cells = no_of_cells
    rrtm_model_data%full_levels = full_levels 
    rrtm_model_data%half_levels = half_levels
    
    rrtm_model_data%nproma       = nproma
    rrtm_model_data%no_of_blocks = ( (no_of_cells - 1) / nproma ) + 1
    rrtm_model_data%end_index    = no_of_cells - ( (rrtm_model_data%no_of_blocks - 1) * nproma )
    
    CALL allocate_rrtm_model_data(rrtm_model_data)
  
  END SUBROUTINE init_rrtm_model_data
  !-----------------------------------------

    
  !-----------------------------------------
  !>
  SUBROUTINE allocate_rrtm_model_data(rrtm_data)

    TYPE(t_rrtm_data), INTENT(inout) :: rrtm_data
    
    INTEGER :: no_of_blocks, nproma, full_levels, half_levels
    INTEGER :: return_status

    nproma       = rrtm_data%nproma
    no_of_blocks = rrtm_data%no_of_blocks
    full_levels  = rrtm_data%full_levels
    half_levels  = rrtm_data%half_levels
    
    !------------------------------------------
    ALLOCATE( &
      & rrtm_data%convection_type(nproma,            no_of_blocks),   &
      & rrtm_data%fr_land_smt    (nproma,            no_of_blocks),   &
      & rrtm_data%fr_glac_smt    (nproma,            no_of_blocks),   &
      & rrtm_data%cosmu0         (nproma,            no_of_blocks),   &
      & rrtm_data%albedo_vis_dir (nproma,            no_of_blocks),   &
      & rrtm_data%albedo_vis_dif (nproma,            no_of_blocks),   &
      & rrtm_data%albedo_nir_dif (nproma,            no_of_blocks),   &
      & rrtm_data%emis_rad       (nproma,            no_of_blocks),   &
      & rrtm_data%tsfctrad       (nproma,            no_of_blocks),   &
      & rrtm_data%pres_ifc       (nproma,half_levels,no_of_blocks),   &
      & rrtm_data%pres           (nproma,full_levels,no_of_blocks),   &
      & rrtm_data%temp           (nproma,full_levels,no_of_blocks),   &
      & rrtm_data%qm_vapor       (nproma,full_levels,no_of_blocks),   &
      & rrtm_data%qm_liquid      (nproma,full_levels,no_of_blocks),   &
      & rrtm_data%qm_ice         (nproma,full_levels,no_of_blocks),   &
      & rrtm_data%qm_o3          (nproma,full_levels,no_of_blocks),   &
      & rrtm_data%acdnc          (nproma,full_levels,no_of_blocks),   &
      & rrtm_data%cld_frc        (nproma,full_levels,no_of_blocks),   &
      & rrtm_data%zaeq1          (nproma,full_levels,no_of_blocks),   &
      & rrtm_data%zaeq2          (nproma,full_levels,no_of_blocks),   &
      & rrtm_data%zaeq3          (nproma,full_levels,no_of_blocks),   &
      & rrtm_data%zaeq4          (nproma,full_levels,no_of_blocks),   &
      & rrtm_data%zaeq5          (nproma,full_levels,no_of_blocks),   &
      & rrtm_data%aclcov         (nproma,            no_of_blocks),   &
      & rrtm_data%lwflxclr       (nproma,half_levels,no_of_blocks),   &
      & rrtm_data%trsolclr       (nproma,half_levels,no_of_blocks),   &
      & rrtm_data%lwflxall       (nproma,half_levels,no_of_blocks),   &
      & rrtm_data%trsolall       (nproma,half_levels,no_of_blocks),   &
      & STAT=return_status)

    IF (return_status /= 0 ) THEN
      CALL finish ("allocate_rrtm_model_data",'failed')
    ENDIF

    
  END SUBROUTINE allocate_rrtm_model_data
  !-----------------------------------------

  !-----------------------------------------
  !>
  SUBROUTINE deallocate_rrtm_model_data(rrtm_data)

    TYPE(t_rrtm_data), INTENT(inout) :: rrtm_data
    
    
    !------------------------------------------
    DEALLOCATE( &
      & rrtm_data%convection_type,   &
      & rrtm_data%fr_land_smt    ,   &
      & rrtm_data%fr_glac_smt    ,   &
      & rrtm_data%cosmu0         ,   &
      & rrtm_data%albedo_vis_dir ,   &
      & rrtm_data%albedo_vis_dif ,   &
      & rrtm_data%albedo_nir_dif ,   &
      & rrtm_data%emis_rad       ,   &
      & rrtm_data%tsfctrad       ,   &
      & rrtm_data%pres_ifc       ,   &
      & rrtm_data%pres           ,   &
      & rrtm_data%temp           ,   &
      & rrtm_data%qm_vapor       ,   &
      & rrtm_data%qm_liquid      ,   &
      & rrtm_data%qm_ice         ,   &
      & rrtm_data%qm_o3          ,   &
      & rrtm_data%acdnc          ,   &
      & rrtm_data%cld_frc        ,   &
      & rrtm_data%zaeq1          ,   &
      & rrtm_data%zaeq2          ,   &
      & rrtm_data%zaeq3          ,   &
      & rrtm_data%zaeq4          ,   &
      & rrtm_data%zaeq5          ,   &
      & rrtm_data%aclcov         ,   &
      & rrtm_data%lwflxclr       ,   &
      & rrtm_data%trsolclr       ,   &
      & rrtm_data%lwflxall       ,   &
      & rrtm_data%trsolall)
    
  END SUBROUTINE deallocate_rrtm_model_data
  !-----------------------------------------


END MODULE mo_rrtm_data_interface
