!>
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
!!
MODULE mo_extpar_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_impl_constants,      ONLY: max_dom

  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist         , &
                                  & open_and_restore_namelist, close_tmpfile

  USE mo_extpar_config,       ONLY: config_itopo                    => itopo             ,           &
                                  & config_fac_smooth_topo          => fac_smooth_topo   ,           &
                                  & config_n_iter_smooth_topo       => n_iter_smooth_topo,           &
                                  & config_hgtdiff_max_smooth_topo  => hgtdiff_max_smooth_topo,      &
                                  & config_l_emiss                  => l_emiss,                      &
                                  & config_heightdiff_threshold     => heightdiff_threshold,         &
                                  & config_extpar_filename          => extpar_filename,              &
                                  & config_extpar_varnames_map_file => extpar_varnames_map_file,     &
                                  & config_lrevert_sea_height       => lrevert_sea_height
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC read_extpar_namelist

  !------------------------------------------------------------------------
  ! Namelist variables
  !------------------------------------------------------------------------
  INTEGER  :: itopo   ! 0: topography specified by analytical functions,
                      ! 1: topography read from netcdf files

  REAL(wp) :: fac_smooth_topo
  INTEGER  :: n_iter_smooth_topo(max_dom)
  REAL(wp) :: hgtdiff_max_smooth_topo(max_dom)
  LOGICAL  :: l_emiss ! if true: read external emissivity map
  REAL(wp) :: heightdiff_threshold(max_dom)
  LOGICAL  :: lrevert_sea_height  ! if true: bring sea points back to original height
  CHARACTER(LEN=filename_max) :: extpar_filename

  ! external parameter: dictionary which maps internal variable names
  ! onto GRIB2 shortnames or NetCDF var names.
  CHARACTER(LEN=filename_max) :: extpar_varnames_map_file

  NAMELIST /extpar_nml/ itopo, fac_smooth_topo,n_iter_smooth_topo,l_emiss, &
                        heightdiff_threshold, extpar_filename,             &
                        extpar_varnames_map_file, hgtdiff_max_smooth_topo, &
                        lrevert_sea_height

CONTAINS
  !>
  !!
  SUBROUTINE read_extpar_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit
    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_extpar_nml:read_extpar_namelist'

    !------------------------------------------------------------
    ! Default settings
    !------------------------------------------------------------
    itopo                   = 0
    fac_smooth_topo         = 0.015625_wp
    n_iter_smooth_topo(:)   = 0
    hgtdiff_max_smooth_topo(:) = 0._wp
    l_emiss                 = .TRUE.
    heightdiff_threshold(:) = 3000._wp
    lrevert_sea_height      = .FALSE.
    extpar_filename         = "<path>extpar_<gridfile>"
    extpar_varnames_map_file = " "

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('extpar_nml')
      READ(funit,NML=extpar_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('extpar_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, extpar_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, extpar_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, extpar_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! Sanity check
    !----------------------------------------------------
    SELECT CASE (itopo)
    CASE (0,1,2) !OK
    CASE default
      CALL finish(TRIM(routine),'Wrong value for itopo. Must be 0 or 1.')
    END SELECT

    !----------------------------------------------------
    ! Fill the configuration state
    !----------------------------------------------------
    config_itopo              = itopo 
    config_fac_smooth_topo    = fac_smooth_topo 
    config_n_iter_smooth_topo = n_iter_smooth_topo
    config_hgtdiff_max_smooth_topo = hgtdiff_max_smooth_topo
    config_l_emiss            = l_emiss
    config_heightdiff_threshold = heightdiff_threshold
    config_lrevert_sea_height = lrevert_sea_height
    config_extpar_filename    = extpar_filename
    config_extpar_varnames_map_file = extpar_varnames_map_file

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=extpar_nml)
      CALL store_and_close_namelist(funit, 'extpar_nml')
    ENDIF
    !write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=extpar_nml)

  END SUBROUTINE read_extpar_namelist

END MODULE mo_extpar_nml
