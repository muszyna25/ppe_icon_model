!>
!! Contains the setup of the variables for io.
!!
!! @par Revision History
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
!!
MODULE mo_io_nml
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
!
  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_char_length, max_ntracer, max_dom, &
    &                              PRES_MSL_METHOD_GME, RH_METHOD_WMO
  USE mo_io_units,           ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio, p_n_work
  USE mo_master_control,     ONLY: is_restart_run
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings
  USE mo_io_config,          ONLY: config_lkeep_in_sync           => lkeep_in_sync          , &
                                 & config_dt_diag                 => dt_diag                , &
                                 & config_dt_checkpoint           => dt_checkpoint          , &
                                 & config_inextra_2d              => inextra_2d             , &
                                 & config_inextra_3d              => inextra_3d             , &
                                 & config_lflux_avg               => lflux_avg              , &
                                 & config_itype_pres_msl          => itype_pres_msl         , &
                                 & config_output_nml_dict         => output_nml_dict        , &
                                 & config_netcdf_dict             => netcdf_dict            , &
                                 & config_lzaxis_reference        => lzaxis_reference       , &
                                 & config_itype_rh                => itype_rh

  USE mo_exception,        ONLY: message, message_text, finish
  USE mo_parallel_config,  ONLY: nproma

  IMPLICIT NONE
  PUBLIC :: read_io_namelist
  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  
CONTAINS
  !>
  !! Read Namelist for I/O.
  !!
  !! This subroutine
  !! - reads the Namelist for I/O
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-06-07)
  !!
  SUBROUTINE read_io_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN)   :: filename
    INTEGER                        :: istat, funit
    INTEGER                        :: iunit

    !-------------------------------------------------------------------------
    ! Namelist variables
    !-------------------------------------------------------------------------

    LOGICAL :: lkeep_in_sync              ! if .true., sync stream after each timestep
    REAL(wp):: dt_diag                    ! diagnostic output timestep [seconds]
    REAL(wp):: dt_checkpoint              ! timestep [seconds] for triggering new restart file

    INTEGER :: inextra_2d                 ! number of extra output fields for debugging
    INTEGER :: inextra_3d                 ! number of extra output fields for debugging
    LOGICAL :: lflux_avg                  ! if .FALSE. the output fluxes are accumulated
                                          !  from the beginning of the run
                                          ! if .TRUE. the output fluxex are average values
                                          !  from the beginning of the run, except of
                                          !  TOT_PREC that would be accumulated
    INTEGER :: itype_pres_msl             ! Specifies method for computation of mean sea level pressure
                                          ! 1: GME-type extrapolation
                                          ! 2: stepwise analytical integration
                                          ! 3: IFS method
                                          ! 4: IFS method with consistency correction
    INTEGER :: itype_rh                   ! Specifies method for computation of relative humidity
                                          ! 1: WMO: water only (e_s=e_s_water)
                                          ! 2: IFS: mixed phases (e_s=a*e_s_water + b*e_s_ice)

    LOGICAL :: lzaxis_reference           ! use ZAXIS_REFERENCE instead of ZAXIS_HYBRID for atmospheric
                                          ! output fields


    CHARACTER(LEN=filename_max) :: &
      &        output_nml_dict,    &     !< maps variable names onto the internal ICON names.
      &        netcdf_dict               !< maps internal variable names onto names in output file (NetCDF only).

    NAMELIST/io_nml/ lkeep_in_sync, dt_diag, dt_checkpoint,  &
      &              inextra_2d, inextra_3d,                 &
      &              lflux_avg, itype_pres_msl, itype_rh,    &
      &              output_nml_dict, netcdf_dict,           &
      &              lzaxis_reference

    !-----------------------
    ! 1. default settings
    !-----------------------
    lkeep_in_sync           = .FALSE.

    dt_diag                 = 86400._wp    !  1 day
    dt_checkpoint           = 2592000._wp  ! 30 days

    inextra_2d              = 0     ! no extra output 2D fields
    inextra_3d              = 0     ! no extra output 3D fields
    lflux_avg               = .TRUE.
    itype_pres_msl          = PRES_MSL_METHOD_GME
    itype_rh                = RH_METHOD_WMO       ! WMO: water only
    output_nml_dict         = ' '
    netcdf_dict             = ' '

    lzaxis_reference        = .TRUE. ! use ZAXIS_REFERENCE (generalVertical)

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('io_nml')
      READ(funit,NML=io_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !-------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('io_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, io_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, io_nml)                                       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, io_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------

    config_lkeep_in_sync           = lkeep_in_sync
    config_dt_diag                 = dt_diag
    config_dt_checkpoint           = dt_checkpoint
    config_inextra_2d              = inextra_2d
    config_inextra_3d              = inextra_3d
    config_lflux_avg               = lflux_avg
    config_itype_pres_msl          = itype_pres_msl
    config_itype_rh                = itype_rh
    config_output_nml_dict         = output_nml_dict
    config_netcdf_dict             = netcdf_dict
    config_lzaxis_reference        = lzaxis_reference
    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=io_nml)
      CALL store_and_close_namelist(funit, 'io_nml')
    ENDIF
    !-----------------------------------------------------
    ! 6. write the contents of the namelist to an ASCII file
    !-----------------------------------------------------
    IF(my_process_is_stdio()) THEN
      WRITE(nnml_output,nml=io_nml)
    END IF

  END SUBROUTINE read_io_namelist

END MODULE mo_io_nml
