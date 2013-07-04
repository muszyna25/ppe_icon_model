!>
!! @brief
!!  Read namelists, make sanity checks specific to each namelist and make
!!  a cross check once all namelists of a component are available.
!!
!! @author
!!  Leonidas Linardakis (MPI-M)
!!  Hui Wan             (MPI-M)
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
  USE mo_gridref_nml         ,ONLY: read_gridref_namelist
  USE mo_dynamics_nml        ,ONLY: read_dynamics_namelist
  USE mo_interpol_nml        ,ONLY: read_interpol_namelist
  USE mo_sleve_nml           ,ONLY: read_sleve_namelist
  USE mo_nh_pzlev_nml        ,ONLY: read_nh_pzlev_namelist
  USE mo_ha_dyn_nml          ,ONLY: read_ha_dyn_namelist
  USE mo_nonhydrostatic_nml  ,ONLY: read_nonhydrostatic_namelist
  USE mo_diffusion_nml       ,ONLY: read_diffusion_namelist

  USE mo_advection_nml       ,ONLY: read_transport_namelist

  USE mo_echam_phy_nml       ,ONLY: read_echam_phy_namelist
  USE mo_nwp_phy_nml         ,ONLY: read_nwp_phy_namelist
  USE mo_radiation_nml       ,ONLY: read_radiation_namelist
  USE mo_vdiff_nml           ,ONLY: read_vdiff_namelist
  USe mo_turbdiff_nml        ,ONLY: read_turbdiff_namelist
  USE mo_echam_conv_nml      ,ONLY: read_echam_conv_namelist
  USE mo_gw_hines_nml        ,ONLY: read_gw_hines_namelist
  USE mo_lnd_nwp_nml         ,ONLY: read_nwp_lnd_namelist
  USE mo_art_nml             ,ONLY: read_art_namelist

  USE mo_initicon_nml        ,ONLY: read_initicon_namelist
  USE mo_ha_testcases        ,ONLY: read_ha_testcase_namelist
  USE mo_nh_testcases        ,ONLY: read_nh_testcase_namelist
  USE mo_meteogram_nml       ,ONLY: read_meteogram_namelist

  USE mo_coupling_nml        ,ONLY: read_coupling_namelist
  USE mo_extpar_nml          ,ONLY: read_extpar_namelist

  USE mo_sea_ice_nml         ,ONLY: read_sea_ice_namelist

  USE mo_name_list_output_init ,ONLY: read_name_list_output_namelists
  USE mo_les_nml             ,ONLY: read_les_namelist
  USE mo_ls_forcing_nml      ,ONLY: read_ls_forcing_namelist
  USE mo_limarea_nml         ,ONLY: read_limarea_namelist
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_atmo_namelists, read_cpl_dummy_namelists

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
    CALL read_nh_pzlev_namelist       (TRIM(atm_namelist_filename))
    CALL read_meteogram_namelist      (TRIM(atm_namelist_filename))
    CALL read_name_list_output_namelists (TRIM(atm_namelist_filename))
    CALL read_dbg_namelist            (TRIM(atm_namelist_filename))

    ! Grid
    !
    CALL read_grid_namelist           (TRIM(atm_namelist_filename))
    CALL read_gridref_namelist        (TRIM(atm_namelist_filename))
    CALL read_interpol_namelist       (TRIM(atm_namelist_filename))
    CALL read_sleve_namelist          (TRIM(atm_namelist_filename))
    !
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
    CALL read_echam_phy_namelist      (TRIM(atm_namelist_filename))
    CALL read_nwp_phy_namelist        (TRIM(atm_namelist_filename))
    CALL read_radiation_namelist      (TRIM(atm_namelist_filename))
    CALL read_vdiff_namelist          (TRIM(atm_namelist_filename))
    CALL read_turbdiff_namelist       (TRIM(atm_namelist_filename))
    CALL read_echam_conv_namelist     (TRIM(atm_namelist_filename))
    CALL read_gw_hines_namelist       (TRIM(atm_namelist_filename))
    CALL read_nwp_lnd_namelist        (TRIM(atm_namelist_filename))
    CALL read_sea_ice_namelist        (TRIM(atm_namelist_filename))
    CALL read_art_namelist            (TRIM(atm_namelist_filename))
    CALL read_les_namelist            (TRIM(atm_namelist_filename))
    CALL read_ls_forcing_namelist     (TRIM(atm_namelist_filename))

    ! Initial conditions
    !
    CALL read_initicon_namelist       (TRIM(atm_namelist_filename))
    CALL read_ha_testcase_namelist    (TRIM(atm_namelist_filename))
    CALL read_nh_testcase_namelist    (TRIM(atm_namelist_filename))

    ! Boundary conditions
    !
    CALL read_extpar_namelist         (TRIM(atm_namelist_filename))
    CALL read_limarea_namelist        (TRIM(atm_namelist_filename))

    !
    ! GRIB output
    CALL read_gribout_namelist        (TRIM(atm_namelist_filename))

    ! Coupling
    !
    CALL read_coupling_namelist       (TRIM(atm_namelist_filename))
    !-----------------------------------------------------------------
    ! Close the file in which all the namelist variables and their
    ! actual values were stored.
    !-----------------------------------------------------------------

    IF (my_process_is_stdio()) CALL close_nml_output

    ! write an annotate table of all namelist settings to a text file
    IF (my_process_is_stdio()) CALL log_nml_settings("nml.atmo.log")

  END SUBROUTINE read_atmo_namelists
  !-------------------------------------------------------------------------


  !---------------------------------------------------------------------
  !>
  !! Read namelists for coupling dummy models
  !!
  SUBROUTINE read_cpl_dummy_namelists(cpl_dummy_namelist,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: cpl_dummy_namelist
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
    CALL read_parallel_namelist       (TRIM(cpl_dummy_namelist))
    CALL read_run_namelist            (TRIM(cpl_dummy_namelist))
    CALL read_io_namelist             (TRIM(cpl_dummy_namelist))
    CALL read_dbg_namelist            (TRIM(cpl_dummy_namelist))

    ! Grid
    !
    CALL read_grid_namelist           (TRIM(cpl_dummy_namelist))
    CALL read_gridref_namelist        (TRIM(cpl_dummy_namelist))
    CALL read_interpol_namelist       (TRIM(cpl_dummy_namelist))
    CALL read_sleve_namelist          (TRIM(cpl_dummy_namelist))
    !
    ! Dynamics
    !
    CALL read_dynamics_namelist       (TRIM(cpl_dummy_namelist))
    CALL read_ha_dyn_namelist         (TRIM(cpl_dummy_namelist))
    CALL read_nonhydrostatic_namelist (TRIM(cpl_dummy_namelist))
    CALL read_diffusion_namelist      (TRIM(cpl_dummy_namelist))

    ! Transport
    !
    CALL read_transport_namelist      (TRIM(cpl_dummy_namelist))

    ! Physics
    !
    CALL read_echam_phy_namelist      (TRIM(cpl_dummy_namelist))
    CALL read_nwp_phy_namelist        (TRIM(cpl_dummy_namelist))
    CALL read_radiation_namelist      (TRIM(cpl_dummy_namelist))
    CALL read_vdiff_namelist          (TRIM(cpl_dummy_namelist))
    CALL read_echam_conv_namelist     (TRIM(cpl_dummy_namelist))
    CALL read_gw_hines_namelist       (TRIM(cpl_dummy_namelist))
    CALL read_nwp_lnd_namelist        (TRIM(cpl_dummy_namelist))
    CALL read_sea_ice_namelist        (TRIM(cpl_dummy_namelist))

    ! Initial conditions
    !
    CALL read_initicon_namelist       (TRIM(cpl_dummy_namelist))
    CALL read_ha_testcase_namelist    (TRIM(cpl_dummy_namelist))
    CALL read_nh_testcase_namelist    (TRIM(cpl_dummy_namelist))
    ! Boundary conditions
    !
    CALL read_extpar_namelist         (TRIM(cpl_dummy_namelist))

    ! Coupling
    !
    CALL read_coupling_namelist       (TRIM(cpl_dummy_namelist))

    !-----------------------------------------------------------------
    ! Close the file in which all the namelist variables and their
    ! actual values were stored.
    !-----------------------------------------------------------------

    IF (my_process_is_stdio()) CALL close_nml_output

    ! write an annotate table of all namelist settings to a text file
    IF (my_process_is_stdio()) CALL log_nml_settings("nml.cpl.log")

  END SUBROUTINE read_cpl_dummy_namelists
  !-------------------------------------------------------------------------

END MODULE mo_read_namelists
