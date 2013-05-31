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
MODULE mo_ocean_read_namelists

  USE mo_mpi                 ,ONLY: my_process_is_stdio
  USE mo_namelist            ,ONLY: open_nml_output, close_nml_output

  USE mo_time_nml            ,ONLY: read_time_namelist

  USE mo_parallel_nml        ,ONLY: read_parallel_namelist
  USE mo_run_nml             ,ONLY: read_run_namelist
  USE mo_io_nml              ,ONLY: read_io_namelist
  USE mo_gribout_nml         ,ONLY: read_gribout_namelist
  USE mo_dbg_nml             ,ONLY: read_dbg_namelist


  USE mo_grid_nml            ,ONLY: read_grid_namelist
  USE mo_gridref_nml         ,ONLY: read_gridref_namelist
  USE mo_dynamics_nml        ,ONLY: read_dynamics_namelist
  USE mo_extpar_nml          ,ONLY: read_extpar_namelist

  USE mo_ocean_nml           ,ONLY: setup_ocean_nml

  USE mo_sea_ice_nml         ,ONLY: read_sea_ice_namelist

  USE mo_name_list_output    ,ONLY: read_name_list_output_namelists
  USE mo_les_nml             ,ONLY: read_les_namelist
#ifndef __ICON_OCEAN__
  USE mo_coupling_nml        ,ONLY: read_coupling_namelist
#endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_ocean_namelists

CONTAINS


  !---------------------------------------------------------------------
  !>
  !! Read namelists for ocean models
  !!
  SUBROUTINE read_ocean_namelists(oce_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: oce_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    !-----------------------------------------------------------------
    ! Create a new file in which all the namelist variables and their
    ! actual values used in the model run will be stored.
    !-----------------------------------------------------------------

    IF(my_process_is_stdio()) CALL open_nml_output('NAMELIST_ICON_output_oce')

    !-----------------------------------------------------------------
    ! Read namelists that are shared by all components of the model.
    ! This means that the same namelists with the same values are
    ! read by all components of a coupled system.
    !-----------------------------------------------------------------

    CALL read_time_namelist           (TRIM(shr_namelist_filename))

    !-----------------------------------------------------------------
    ! Read namelists that are specific to the oce model.
    ! In case of a coupled simulation, the atmosphere model may also
    ! read some of these namelists, but probably from a different
    ! ASCII file containing different values.
    !-----------------------------------------------------------------

    ! General
    !
    CALL read_parallel_namelist       (TRIM(oce_namelist_filename))
    CALL read_run_namelist            (TRIM(oce_namelist_filename))
    CALL read_io_namelist             (TRIM(oce_namelist_filename))
    CALL read_name_list_output_namelists (TRIM(oce_namelist_filename))
    CALL read_dbg_namelist            (TRIM(oce_namelist_filename))

    ! Grid
    !
    CALL read_grid_namelist           (TRIM(oce_namelist_filename))
!    CALL read_interpol_namelist       (TRIM(oce_namelist_filename))

    ! Dynamics, transport and physics
    ! (still using the old setup)
    !
    CALL read_dynamics_namelist       (TRIM(oce_namelist_filename))

    ! Boundary conditions
    !
    CALL read_extpar_namelist         (TRIM(oce_namelist_filename))

    !
    ! GRIB2 output
    CALL read_gribout_namelist        (TRIM(oce_namelist_filename))

    ! Coupling
    !
#ifndef __ICON_OCEAN__
    CALL read_coupling_namelist       (TRIM(oce_namelist_filename))
#endif

    ! More namelists from the old setup
    !
    CALL setup_ocean_nml              (TRIM(oce_namelist_filename))

    ! Sea ice namelist
    !
    CALL read_sea_ice_namelist        (TRIM(oce_namelist_filename))

    !-----------------------------------------------------------------
    ! Close the file in which all the namelist variables and their
    ! actual values were stored.
    !-----------------------------------------------------------------

    IF (my_process_is_stdio()) CALL close_nml_output

  END SUBROUTINE read_ocean_namelists
  !-------------------------------------------------------------------------


END MODULE mo_ocean_read_namelists
