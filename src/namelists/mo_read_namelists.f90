!>
!! @brief
!!  Read namelists, make sanity checks specific to each namelist and make
!!  a cross check once all namelists of a component are available.
!!
!! @author
!!  Leonidas Linardakis (MPI-M)
!!  Hui Wan             (MPI-M)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_read_namelists

  USE mo_mpi                 ,ONLY: my_process_is_stdio
  USE mo_namelist            ,ONLY: open_nml_output, close_nml_output
  USE mo_nml_annotate        ,ONLY: log_nml_settings

  USE mo_time_nml            ,ONLY: read_time_namelist

  USE mo_parallel_nml        ,ONLY: read_parallel_namelist
  USE mo_run_nml             ,ONLY: read_run_namelist
  USE mo_io_nml              ,ONLY: read_io_namelist
  USE mo_gribout_nml         ,ONLY: read_gribout_namelist
  USE mo_dbg_nml             ,ONLY: read_dbg_namelist


  USE mo_grid_nml            ,ONLY: read_grid_namelist
  USE mo_grid_config         ,ONLY: init_grid_configuration
  USE mo_gridref_nml         ,ONLY: read_gridref_namelist
  USE mo_dynamics_nml        ,ONLY: read_dynamics_namelist
  USE mo_interpol_nml        ,ONLY: read_interpol_namelist
  USE mo_sleve_nml           ,ONLY: read_sleve_namelist
  USE mo_ha_dyn_nml          ,ONLY: read_ha_dyn_namelist
  USE mo_nonhydrostatic_nml  ,ONLY: read_nonhydrostatic_namelist
  USE mo_diffusion_nml       ,ONLY: read_diffusion_namelist

  USE mo_advection_nml       ,ONLY: read_transport_namelist

  USE mo_echam_phy_nml       ,ONLY: process_echam_phy_nml
  USE mo_echam_cld_nml       ,ONLY: process_echam_cld_nml
  USE mo_echam_cnv_nml       ,ONLY: process_echam_cnv_nml
  USE mo_echam_gwd_nml       ,ONLY: process_echam_gwd_nml
  USE mo_echam_rad_nml       ,ONLY: process_echam_rad_nml
  USE mo_echam_sso_nml       ,ONLY: process_echam_sso_nml
  USE mo_echam_vdf_nml       ,ONLY: process_echam_vdf_nml
  
  USE mo_nwp_phy_nml         ,ONLY: read_nwp_phy_namelist
  USE mo_nwp_tuning_nml      ,ONLY: read_nwp_tuning_namelist
  USE mo_ensemble_pert_nml   ,ONLY: read_ensemble_pert_namelist
  USE mo_radiation_nml       ,ONLY: read_radiation_namelist
  USE mo_psrad_radiation     ,ONLY: setup_psrad_radiation
  USE mo_ccycle_nml          ,ONLY: read_ccycle_nml
  USE mo_ccycle_config       ,ONLY: init_ccycle_config
  USE mo_synsat_nml          ,ONLY: read_synsat_namelist
  USE mo_turbdiff_nml        ,ONLY: read_turbdiff_namelist
  USE mo_lnd_nwp_nml         ,ONLY: read_nwp_lnd_namelist
  USE mo_art_nml             ,ONLY: read_art_namelist

  USE mo_initicon_nml        ,ONLY: read_initicon_namelist
  USE mo_ha_testcases        ,ONLY: read_ha_testcase_namelist
  USE mo_nh_testcases_nml    ,ONLY: read_nh_testcase_namelist
  USE mo_meteogram_nml       ,ONLY: read_meteogram_namelist

  USE mo_coupling_nml        ,ONLY: read_coupling_namelist
  USE mo_extpar_nml          ,ONLY: read_extpar_namelist

  USE mo_sea_ice_nml         ,ONLY: read_sea_ice_namelist

  USE mo_name_list_output_init ,ONLY: read_name_list_output_namelists
  USE mo_les_nml             ,ONLY: read_les_namelist
  USE mo_ls_forcing_nml      ,ONLY: read_ls_forcing_namelist
  USE mo_limarea_nml         ,ONLY: read_limarea_namelist

  USE mo_run_config          ,ONLY: iforcing
  USE mo_impl_constants      ,ONLY: IECHAM, ILDF_ECHAM, INWP
  USE mo_assimilation_nml    ,ONLY: read_assimilation_namelist
  USE mo_nudging_nml         ,ONLY: read_nudging_namelist

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_atmo_namelists

CONTAINS

  !---------------------------------------------------------------------
  !>
  !! Read namelists for atmospheric models
  !!
  SUBROUTINE read_atmo_namelists(atm_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: atm_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    !-----------------------------------------------------------------
    ! Create a new file in which all the namelist variables and their
    ! actual values used in the model run will be stored.
    !-----------------------------------------------------------------

    IF(my_process_is_stdio()) CALL open_nml_output('NAMELIST_ICON_output_atm')

    !-----------------------------------------------------------------
    ! Read namelists that are shared by all components of the model.
    ! This means that the same namelists with the same values are
    ! read by all components of a coupled system.
    !-----------------------------------------------------------------

    CALL read_time_namelist           (TRIM(shr_namelist_filename))

    !-----------------------------------------------------------------
    ! Read namelist that are specific to the atm model.
    ! In case of a coupled simulation, the ocean model may also
    ! read some of these namelists, but probably from a different
    ! ASCII file containing different values.
    !-----------------------------------------------------------------

    ! General
    !
    CALL read_parallel_namelist       (TRIM(atm_namelist_filename))
    CALL read_run_namelist            (TRIM(atm_namelist_filename))
    CALL read_io_namelist             (TRIM(atm_namelist_filename))
    CALL read_meteogram_namelist      (TRIM(atm_namelist_filename))
    CALL read_name_list_output_namelists (TRIM(atm_namelist_filename))
    CALL read_dbg_namelist            (TRIM(atm_namelist_filename))
    CALL read_synsat_namelist         (TRIM(atm_namelist_filename))
    
    ! Grid
    !
    CALL read_grid_namelist           (TRIM(atm_namelist_filename))
    CALL read_gridref_namelist        (TRIM(atm_namelist_filename))
    CALL read_interpol_namelist       (TRIM(atm_namelist_filename))
    CALL read_sleve_namelist          (TRIM(atm_namelist_filename))
    !
    CALL init_grid_configuration()    ! so that the number of grids is known
    !                                 ! and arrays can be allocated

    ! Dynamics
    !
    CALL read_dynamics_namelist       (TRIM(atm_namelist_filename))
    CALL read_ha_dyn_namelist         (TRIM(atm_namelist_filename))
    CALL read_nonhydrostatic_namelist (TRIM(atm_namelist_filename))
    CALL read_diffusion_namelist      (TRIM(atm_namelist_filename))

    ! Transport
    !
    CALL read_transport_namelist      (TRIM(atm_namelist_filename))

    ! Physics
    !
    SELECT CASE (iforcing)
    CASE (IECHAM, ILDF_ECHAM)
       !
       ! ECHAM physics ...
       CALL process_echam_phy_nml        (TRIM(atm_namelist_filename))
       !
       ! ... and the employed parameterizations
       CALL process_echam_cld_nml        (TRIM(atm_namelist_filename))
       CALL process_echam_cnv_nml        (TRIM(atm_namelist_filename))
       CALL process_echam_gwd_nml        (TRIM(atm_namelist_filename))
       CALL process_echam_rad_nml        (TRIM(atm_namelist_filename))
       CALL process_echam_sso_nml        (TRIM(atm_namelist_filename))
       CALL process_echam_vdf_nml        (TRIM(atm_namelist_filename))
       !
       CALL read_sea_ice_namelist        (TRIM(atm_namelist_filename))
       CALL read_art_namelist            (TRIM(atm_namelist_filename))
       ! setup_psrad_radiation depends on cloud_config
       CALL setup_psrad_radiation        (TRIM(atm_namelist_filename))
       ! carbon cycle
       CALL init_ccycle_config
       CALL read_ccycle_nml              (TRIM(atm_namelist_filename))
       !
    CASE (INWP)
       !
       CALL read_nwp_phy_namelist        (TRIM(atm_namelist_filename))
       CALL read_nwp_tuning_namelist     (TRIM(atm_namelist_filename))
       CALL read_ensemble_pert_namelist  (TRIM(atm_namelist_filename))
       CALL read_radiation_namelist      (TRIM(atm_namelist_filename))
       CALL read_turbdiff_namelist       (TRIM(atm_namelist_filename))
       CALL read_nwp_lnd_namelist        (TRIM(atm_namelist_filename))
       CALL read_art_namelist            (TRIM(atm_namelist_filename))
       CALL read_les_namelist            (TRIM(atm_namelist_filename))
       CALL read_ls_forcing_namelist     (TRIM(atm_namelist_filename))
       !
    END SELECT

    ! Initial conditions
    !
    CALL read_initicon_namelist       (TRIM(atm_namelist_filename))
    CALL read_ha_testcase_namelist    (TRIM(atm_namelist_filename))
    CALL read_nh_testcase_namelist    (TRIM(atm_namelist_filename))

    ! Boundary conditions
    !
    CALL read_extpar_namelist         (TRIM(atm_namelist_filename))
    CALL read_limarea_namelist        (TRIM(atm_namelist_filename))
    CALL read_nudging_namelist        (TRIM(atm_namelist_filename))

    ! GRIB output
    !
    CALL read_gribout_namelist        (TRIM(atm_namelist_filename))

    ! Coupling
    !
    CALL read_coupling_namelist       (TRIM(atm_namelist_filename))

    ! Assimilation
    CALL read_assimilation_namelist   (TRIM(atm_namelist_filename))
    !-----------------------------------------------------------------
    ! Close the file in which all the namelist variables and their
    ! actual values were stored.
    !-----------------------------------------------------------------

    IF (my_process_is_stdio()) CALL close_nml_output

    ! write an annotate table of all namelist settings to a text file
    IF (my_process_is_stdio()) CALL log_nml_settings("nml.atmo.log")

  END SUBROUTINE read_atmo_namelists
  !-------------------------------------------------------------------------


END MODULE mo_read_namelists
