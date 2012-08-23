!>
!! This module is the data communication interface between nwp_rrtm_interface and RRTM
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

  USE mo_exception,           ONLY: message, get_filename_noext, finish !message_tex
  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: find_next_free_unit, filename_max
  
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_lnd_diag
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_nonhydro_types,      ONLY: t_nh_diag
  
  USE mo_parallel_config,     ONLY: parallel_radiation_mode, radiation_division_file_name, nproma
  USE mo_mpi,                 ONLY: my_process_is_mpi_seq, my_process_is_mpi_workroot,         &
    & num_work_procs, get_my_mpi_work_id, get_my_mpi_work_comm_size, get_mpi_all_workroot_id,  &
    & get_my_mpi_work_communicator, p_bcast
  USE mo_icon_comm_lib,       ONLY: new_icon_comm_pattern, inverse_of_icon_comm_pattern,   &
    & print_grid_comm_pattern, new_icon_comm_variable, is_ready, until_sync, icon_comm_sync_all 


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_rrtm_data
  PUBLIC :: construct_rrtm_model_repart
  PUBLIC :: destruct_rrtm_model_repart
  PUBLIC :: recv_radiation_input
    
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

    ! communication patterns
    INTEGER, POINTER :: global_index(:), dynamics_owner(:)    
    INTEGER :: radiation_recv_comm_pattern, radiation_send_comm_pattern
    
    !------------------------------------------
    ! grid parameters
    INTEGER   :: no_of_cells
    INTEGER   :: full_levels, half_levels
    ! blocking parameters
    INTEGER   :: block_size
    INTEGER   :: no_of_blocks, end_index

    
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
    REAL(wp), POINTER ::  albedo_nir_dir(:,:) ! surface albedo for visible and near IR range, direct
    REAL(wp), POINTER ::  albedo_nir_dif(:,:) !< in surface albedo for visible range, diffuse
    REAL(wp), POINTER ::  albedo_vis_dif(:,:) !< in surface albedo for visible range, diffuse
    
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

  TYPE(t_rrtm_data), TARGET :: rrtm_model_data
    

CONTAINS



  !-----------------------------------------
  !>
  SUBROUTINE construct_rrtm_model_repart(patch)
    TYPE(t_patch), INTENT(in), TARGET :: patch

    INTEGER, POINTER :: radiation_owner(:)
    INTEGER :: my_radiation_cells
    
    INTEGER :: my_mpi_work_id, my_mpi_work_comm_size, workroot_mpi_id, my_mpi_work_communicator    
    INTEGER :: return_status, file_id, j
    
    CHARACTER(*), PARAMETER :: method_name = "construct_rrtm_model_repart"
    
    IF (parallel_radiation_mode <= 0) RETURN
    IF (my_process_is_mpi_seq()) THEN
      parallel_radiation_mode = 0
      RETURN
    ENDIF

    my_mpi_work_id        = get_my_mpi_work_id()
    my_mpi_work_comm_size = get_my_mpi_work_comm_size()
    workroot_mpi_id       = get_mpi_all_workroot_id()
    my_mpi_work_communicator = get_my_mpi_work_communicator()

    ! NOTE: this global array should be eventully replaced by local arrays
    !       this is a first implementation
    ALLOCATE(radiation_owner(patch%n_patch_cells_g), STAT=return_status)
    IF (return_status /= 0 ) &
      CALL finish (method_name,'ALLOCATE(radiation_owner) failed')
    
    ! read the radiation re-distribution from file
    IF (my_mpi_work_id == workroot_mpi_id) THEN

      ! read all the radiation owners and sent the info to each procs
        
      IF (radiation_division_file_name == "") THEN
        radiation_division_file_name = &
          & TRIM(get_filename_noext(patch%grid_filename))//'.radiation_cell_domain_ids'
      ENDIF
        
      WRITE(0,*) "Read radiation redistribution from file: ", TRIM(radiation_division_file_name), &
        & " total cells=", patch%n_patch_cells_g
      file_id = find_next_free_unit(10,900)
      OPEN(file_id,FILE=TRIM(radiation_division_file_name), STATUS='OLD', IOSTAT=return_status)
      IF(return_status /= 0) CALL finish(method_name,&
        & 'Unable to open input file: '//TRIM(radiation_division_file_name))
        
      DO j = 1, patch%n_patch_cells_g
        READ(file_id,*,IOSTAT=return_status) radiation_owner(j)
        IF(return_status /= 0) CALL finish(method_name,&
          & 'Error reading: '//TRIM(radiation_division_file_name))
      ENDDO
      
      CLOSE(file_id)

      ! some checks
      IF ( MINVAL(radiation_owner(:)) < 0 .OR. &
         & MAXVAL(radiation_owner(:)) >= my_mpi_work_comm_size) THEN
        WRITE(0,*) "mpi_work_comm_size=",my_mpi_work_comm_size, &
        & " MINAVAL=", MINVAL(radiation_owner(:)), &
          " MAXVAL=",  MAXVAL(radiation_owner(:))
        CALL finish(method_name,'Invalid redistribution')
      ENDIF
      
    ENDIF
      
    CALL p_bcast(radiation_owner, workroot_mpi_id, comm=my_mpi_work_communicator)

    !-------------------------------
    ! get my radiation cells
    my_radiation_cells = 0
    DO j=1, patch%n_patch_cells_g
      IF (radiation_owner(j) == my_mpi_work_id) THEN
        my_radiation_cells = my_radiation_cells + 1
        ! keep the global index in the first part of the radiation_owner
        radiation_owner(my_radiation_cells) = j
      ENDIF
    ENDDO
    
    ALLOCATE(rrtm_model_data%global_index(my_radiation_cells), STAT=return_status)
    IF (return_status /= 0 ) &
      CALL finish (method_name,'ALLOCATE(global_index) failed')
    rrtm_model_data%global_index(1:my_radiation_cells) = radiation_owner(1:my_radiation_cells)
    DEALLOCATE(radiation_owner)

    !-------------------------------
    ! get the dymamics owners of my radiation cells
    ALLOCATE(rrtm_model_data%dynamics_owner(my_radiation_cells), STAT=return_status)
    DO j=1, my_radiation_cells
      rrtm_model_data%dynamics_owner(j) = &
      & patch%cells%owner_g(rrtm_model_data%global_index(j))
    ENDDO

    ! create the receive communicator from the dynamics
    rrtm_model_data%radiation_recv_comm_pattern =            &
      & new_icon_comm_pattern(                               &
      & total_no_of_points = my_radiation_cells,             &
      & receive_from_owner = rrtm_model_data%dynamics_owner, &
      & my_global_index = rrtm_model_data%global_index,      &
      & owners_local_index = patch%cells%loc_index,          &
      & allow_send_to_myself = .true. ,                      &
      & name = "radiation_rcv_from_dynamics" )

    ! create the inverse communicator, from radiation to dynamics
    rrtm_model_data%radiation_send_comm_pattern =            &
      & inverse_of_icon_comm_pattern(                        &
      &  in_comm_pattern_id = rrtm_model_data%radiation_recv_comm_pattern, &
      &  name = "radiation_send_to_dynamics" )

    CALL print_grid_comm_pattern(rrtm_model_data%radiation_recv_comm_pattern)
    CALL print_grid_comm_pattern(rrtm_model_data%radiation_send_comm_pattern)

    ! now allocate the data for the radiation interface
    CALL init_rrtm_data( &
      & rrtm_data   = rrtm_model_data   , &
      & no_of_cells = my_radiation_cells, &
      & full_levels = patch%nlev,         &
      & half_levels = patch%nlevp1,       &
      & block_size  = nproma)

  END SUBROUTINE construct_rrtm_model_repart
  !-----------------------------------------

  !-----------------------------------------
  !>
  SUBROUTINE init_rrtm_data(rrtm_data, no_of_cells, full_levels, half_levels, block_size)
    INTEGER, INTENT(in) :: no_of_cells, full_levels, half_levels, block_size
    TYPE(t_rrtm_data), INTENT(inout) :: rrtm_data
    
    ! fill the rrtm_data parameters
    rrtm_data%no_of_cells = no_of_cells
    rrtm_data%full_levels = full_levels
    rrtm_data%half_levels = half_levels
    
    rrtm_data%block_size       = block_size
    rrtm_data%no_of_blocks = ( (no_of_cells - 1) / block_size ) + 1
    rrtm_data%end_index    = no_of_cells - ( (rrtm_data%no_of_blocks - 1) * block_size )
    
    CALL allocate_rrtm_model_data(rrtm_data)

    rrtm_data%convection_type(:,:) = 0 ! this is always 0 and will not ber communicated
      
  END SUBROUTINE init_rrtm_data
  !-----------------------------------------

    
  !-----------------------------------------
  !>
  SUBROUTINE allocate_rrtm_model_data(rrtm_data)

    TYPE(t_rrtm_data), INTENT(inout) :: rrtm_data
    
    INTEGER :: no_of_blocks, block_size, full_levels, half_levels
    INTEGER :: return_status

    block_size   = rrtm_data%block_size
    no_of_blocks = rrtm_data%no_of_blocks
    full_levels  = rrtm_data%full_levels
    half_levels  = rrtm_data%half_levels
    
    !------------------------------------------
    ALLOCATE( &
      & rrtm_data%convection_type(block_size,            no_of_blocks),   &
      & rrtm_data%fr_land_smt    (block_size,            no_of_blocks),   &
      & rrtm_data%fr_glac_smt    (block_size,            no_of_blocks),   &
      & rrtm_data%cosmu0         (block_size,            no_of_blocks),   &
      & rrtm_data%albedo_vis_dir (block_size,            no_of_blocks),   &
      & rrtm_data%albedo_nir_dir (block_size,            no_of_blocks),   &
      & rrtm_data%albedo_vis_dif (block_size,            no_of_blocks),   &
      & rrtm_data%albedo_nir_dif (block_size,            no_of_blocks),   &
      & rrtm_data%emis_rad       (block_size,            no_of_blocks),   &
      & rrtm_data%tsfctrad       (block_size,            no_of_blocks),   &
      & rrtm_data%pres_ifc       (block_size,half_levels,no_of_blocks),   &
      & rrtm_data%pres           (block_size,full_levels,no_of_blocks),   &
      & rrtm_data%temp           (block_size,full_levels,no_of_blocks),   &
      & rrtm_data%qm_vapor       (block_size,full_levels,no_of_blocks),   &
      & rrtm_data%qm_liquid      (block_size,full_levels,no_of_blocks),   &
      & rrtm_data%qm_ice         (block_size,full_levels,no_of_blocks),   &
      & rrtm_data%qm_o3          (block_size,full_levels,no_of_blocks),   &
      & rrtm_data%acdnc          (block_size,full_levels,no_of_blocks),   &
      & rrtm_data%cld_frc        (block_size,full_levels,no_of_blocks),   &
      & rrtm_data%zaeq1          (block_size,full_levels,no_of_blocks),   &
      & rrtm_data%zaeq2          (block_size,full_levels,no_of_blocks),   &
      & rrtm_data%zaeq3          (block_size,full_levels,no_of_blocks),   &
      & rrtm_data%zaeq4          (block_size,full_levels,no_of_blocks),   &
      & rrtm_data%zaeq5          (block_size,full_levels,no_of_blocks),   &
      & rrtm_data%aclcov         (block_size,            no_of_blocks),   &
      & rrtm_data%lwflxclr       (block_size,half_levels,no_of_blocks),   &
      & rrtm_data%trsolclr       (block_size,half_levels,no_of_blocks),   &
      & rrtm_data%lwflxall       (block_size,half_levels,no_of_blocks),   &
      & rrtm_data%trsolall       (block_size,half_levels,no_of_blocks),   &
      & STAT=return_status)

    IF (return_status /= 0 ) THEN
      CALL finish ("allocate_rrtm_model_data",'failed')
    ENDIF

    
  END SUBROUTINE allocate_rrtm_model_data
  !-----------------------------------------

  !-----------------------------------------
  !>
  SUBROUTINE destruct_rrtm_model_repart()
    
    IF (parallel_radiation_mode <= 0) RETURN
    IF (my_process_is_mpi_seq()) THEN
      parallel_radiation_mode = 0
      RETURN
    ENDIF

    CALL deallocate_rrtm_model_data(rrtm_model_data)
  
  END SUBROUTINE destruct_rrtm_model_repart
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
      & rrtm_data%albedo_nir_dir ,   &
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

  !-----------------------------------------
  SUBROUTINE recv_radiation_input( &
    & zland      ,&!< in     land fraction
    & zglac      ,&!< in     land glacier fraction
    & cos_mu0    ,&!< in  cos of zenith angle mu0
    & alb_vis_dir,&!< in surface albedo for visible range, direct
    & alb_nir_dir,&!< in surface albedo for near IR range, direct
    & alb_vis_dif,&!< in surface albedo for visible range, diffuse
    & alb_nir_dif,&!< in surface albedo for near IR range, diffuse
    & emis_rad   ,&!< in longwave surface emissivity
    & tk_sfc     ,&!< in surface temperature
    & pp_hl      ,&!< in  pres at half levels at t-dt [Pa]
    & pp_fl      ,&!< in  pres at full levels at t-dt [Pa]
    & tk_fl      ,&!< in  temperature at full level at t-dt
    & qm_vap     ,&!< in  water vapor mass mix ratio at t-dt
    & qm_liq     ,&!< in cloud water mass mix ratio at t-dt
    & qm_ice     ,&!< in cloud ice mass mixing ratio at t-dt
    & qm_o3      ,&!< in o3 mass mixing ratio at t-dt
    & cdnc       ,&!< in  cloud droplet numb conc. [1/m**3]
    & cld_frc    ,&!< in  cloud fraction [m2/m2]
    & zaeq1      ,&!< in aerosol continental
    & zaeq2      ,&!< in aerosol maritime
    & zaeq3      ,&!< in aerosol urban
    & zaeq4      ,&!< in aerosol volcano ashes
    & zaeq5      ,&!< in aerosol stratospheric background
    & patch      ,&!< in
    & rrtm_data)   !< pointer out
    
    REAL(wp), TARGET ::  zland(:,:)   !< fraction land in a grid element        [ ]
                                                        !  = smoothed fr_land
    REAL(wp), TARGET ::  zglac(:,:)   !< fraction land glacier in a grid element [ ]
                                                        ! = smoothed fr_glac
    REAL(wp), TARGET ::  cos_mu0(:,:)        ! cosine of solar zenith angle

    REAL(wp), TARGET ::  alb_vis_dir(:,:) ! surface albedo for visible and near IR range, direct
    REAL(wp), TARGET ::  alb_nir_dir(:,:) ! surface albedo for visible and near IR range, direct
    REAL(wp), TARGET ::  alb_nir_dif(:,:) !< in surface albedo for visible range, diffuse
    REAL(wp), TARGET ::  alb_vis_dif(:,:) !< in surface albedo for visible range, diffuse

    REAL(wp), TARGET ::  emis_rad(:,:)      ! lw sfc emissivity
    REAL(wp), TARGET ::  tk_sfc(:,:)      ! surface temperature at trad [K]

    REAL(wp), TARGET ::  pp_hl(:,:,:)    ! pressure at interfaces (nproma,nlevp1,nblks_c)  [Pa]
    REAL(wp), TARGET ::  pp_fl(:,:,:)        ! pressure (nproma,nlev,nblks_c)                  [Pa]
    REAL(wp), TARGET ::  tk_fl(:,:,:)        ! temperature (nproma,nlev,nblks_c)                 [K]

    REAL(wp), TARGET ::  qm_vap(:,:,:)    !< water vapor mass mix ratio at t-dt
    REAL(wp), TARGET ::  qm_liq(:,:,:)   !< cloud water mass mix ratio at t-dt
    REAL(wp), TARGET ::  qm_ice(:,:,:)      !< cloud ice mass mixing ratio at t-dt

    REAL(wp), TARGET ::  qm_o3(:,:,:)       ! in o3 mass mixing ratio
    REAL(wp), TARGET ::  cdnc(:,:,:)       ! cloud droplet number concentration [1/m**3]
    REAL(wp), TARGET ::  cld_frc(:,:,:)     !< cloud fraction [m2/m2]

    REAL(wp), TARGET ::  zaeq1(:,:,:)       !< in aerosol continental
    REAL(wp), TARGET ::  zaeq2(:,:,:)       !< in aerosol maritime
    REAL(wp), TARGET ::  zaeq3(:,:,:)       !< in aerosol urban
    REAL(wp), TARGET ::  zaeq4(:,:,:)       !< in aerosol volcano ashes
    REAL(wp), TARGET ::  zaeq5(:,:,:)       !< in aerosol stratospheric background

    TYPE(t_patch) :: patch
    TYPE(t_rrtm_data), POINTER :: rrtm_data

    INTEGER :: recv_comm_pattern
    INTEGER :: recv_fr_land_smt
    INTEGER :: recv_tmp

    recv_comm_pattern = rrtm_model_data%radiation_recv_comm_pattern
! -----------------------------------
    recv_fr_land_smt = new_icon_comm_variable ( &
      & recv_var = rrtm_model_data%fr_land_smt, &
      & send_var = zland,                       &
      & comm_pattern_index = recv_comm_pattern, &
      & status   = is_ready,                    &
      & scope    = until_sync,                  &
      & name     = "fr_land" )

      
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%fr_glac_smt   , &
      &  send_var = zglac,                          &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
! 
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%cosmu0        , &
      &  send_var = cos_mu0,                        &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
! 
    
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%albedo_vis_dir, &
      &  send_var = alb_vis_dir,                    &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
          
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%albedo_nir_dir, &
      &  send_var = alb_nir_dir,                    &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
      
         
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%albedo_vis_dif, &
      &  send_var = alb_vis_dif,                    &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
          
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%albedo_nir_dif, &
      &  send_var = alb_nir_dif,                    &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
! 
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%emis_rad      , &
      &  send_var = emis_rad,                       &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
              
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%tsfctrad      , &
      &  send_var = tk_sfc,                         &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
! 
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%pres_ifc      , &
      &  send_var = pp_hl,                          &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
          
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%pres          , &
      &  send_var = pp_fl,                          &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
          
     recv_tmp     = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%temp          , &
      &  send_var = tk_fl,                          &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
! 
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%qm_vapor      , &
      &  send_var = qm_vap,                         &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
    
          
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%qm_liquid     , &
      &  send_var = qm_liq,                        &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
          
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%qm_ice        , &
      &  send_var = qm_ice,                         &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
! 
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%qm_o3         , &
      &  send_var = qm_o3,                          &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
          
    recv_tmp      = new_icon_comm_variable ( & 
      &  recv_var = rrtm_model_data%acdnc         , &
      &  send_var = cdnc,                           &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
          
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%cld_frc       , &
      &  send_var = cld_frc,                        &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
! 
    
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%zaeq1         , &
      &  send_var = zaeq1,                          &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
          
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%zaeq2         , &
      &  send_var = zaeq2,                          &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
          
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%zaeq3         , &
      &  send_var = zaeq3,                          &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
          
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%zaeq4         , &
      &  send_var = zaeq4,                          &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )
          
    recv_tmp      = new_icon_comm_variable ( &
      &  recv_var = rrtm_model_data%zaeq5         , &
      &  send_var = zaeq5,                          &
      &  comm_pattern_index = recv_comm_pattern,    &
      &  status   = is_ready,                       &
      &  scope    = until_sync,                     &
      &  name     = "tmp" )

    CALL icon_comm_sync_all
  

  END SUBROUTINE recv_radiation_input

END MODULE mo_rrtm_data_interface
