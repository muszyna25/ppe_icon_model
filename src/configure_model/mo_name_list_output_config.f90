!>
!! @brief configuration setup for output name lists
!!
!! Configuration setup for output name lists
!!
!! @author R.Johanni, F. Prill
!!
!! @note This is only a preliminary implementation of a configure state
!!       for output name lists; based on R. Johanni's name list
!!       handling which was originally implemented in
!!       "shared/mo_name_list_output"
!!
!! @par Revision History
!! Moved configure state from shared/mo_name_list_output:
!! F. Prill, DWD (2012-01-26)
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
MODULE mo_name_list_output_config

  USE mo_kind,                  ONLY: wp
  USE mo_io_units,              ONLY: filename_max
  USE mo_impl_constants,        ONLY: max_phys_dom
  USE mo_cdi_constants,         ONLY: FILETYPE_GRB, FILETYPE_GRB2

  IMPLICIT NONE
  PUBLIC

  PUBLIC :: is_grib_output

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  ! Flag whether name_list output is active, i.e. at least one /output_nml/ has been read

  LOGICAL :: name_list_output_active = .FALSE.

  ! Flag whether async name_list I/O is used, it is set in the main program:

  LOGICAL :: use_async_name_list_io = .FALSE.

  ! The following parameter decides whether physical or logical patches are output
  ! and thus whether the domain number in output name lists pertains to physical
  ! or logical patches.

  LOGICAL, PARAMETER :: l_output_phys_patch = .TRUE.

  INTEGER, PARAMETER :: &
    max_var_ml = 100, & ! maximum number of output model-level variables
    max_var_pl = 100, & ! maximum number of pressure-level variables
    max_var_hl = 100, & ! maximum number of height-level variables
    max_bounds = 100, & ! maximum number of output_bounds
    max_levels = 100, & ! maximum number of pressure/height levels
    vname_len  =  32    ! variable name length in I/O namelists

  TYPE t_output_name_list

    INTEGER  :: filetype            ! One of CDI's FILETYPE_XXX constants
    CHARACTER(LEN=8) :: namespace   ! 'DWD' - DWD short names (or 'MPIM', 'CMIP', 'ECMWF')
    CHARACTER(LEN=filename_max) :: map_file ! File containig mapping internal names -> names in NetCDF
    INTEGER  :: mode                ! 1 = forecast mode, 2 = climate mode
    INTEGER  :: dom(max_phys_dom)   ! domains for which this namelist is used, ending with -1
    INTEGER  :: output_time_unit    ! 1 = second, 2=minute, 3=hour, 4=day, 5=month, 6=year
    REAL(wp) :: output_bounds(3,max_bounds) ! post-processing times in units defined by output_time_unit: start, end, increment
    INTEGER  :: steps_per_file      ! Max number of output steps in one output file
    LOGICAL  :: include_last        ! Flag whether to include the last timestep in output
    LOGICAL  :: output_grid         ! Flag whether grid information is output (in NetCDF output)
    CHARACTER(LEN=filename_max) :: output_filename   ! output filename prefix
    LOGICAL  :: lwrite_ready        ! Flag. TRUE if a "ready file" (sentinel file) should be written at the end of each output stage
    CHARACTER(LEN=filename_max) :: ready_directory        ! output directory for ready files
    CHARACTER(LEN=vname_len)  :: ml_varlist(max_var_ml)   ! name of model level fields (translation to model by namespace)
    CHARACTER(LEN=vname_len)  :: pl_varlist(max_var_pl)   ! name of pressure level fields (translation to model by namespace)
    REAL(wp) :: p_levels(max_levels)                      ! pressure levels [hPa]
    CHARACTER(LEN=vname_len)  :: hl_varlist(max_var_hl)   ! name of height level fields
    REAL(wp) :: h_levels(max_levels)                      ! height levels
    INTEGER  :: remap               ! interpolate horizontally, 0: none, 1: to regular lat-lon grid, 2: to Gaussian grids, (3:...)
    LOGICAL  :: remap_internal      ! do interpolations online in the model or external (including triggering)
    REAL(wp) :: reg_lon_def(3)      ! if remap=1: start, increment, end longitude in degrees
    REAL(wp) :: reg_lat_def(3)      ! if remap=1: start, increment, end latitude in degrees
    INTEGER  :: gauss_tgrid_def     ! if remap=2: triangular truncation (e.g.63 for T63) for which the Gauss grid should be used
    REAL(wp) :: north_pole(2)       ! definition of north pole for rotated lon-lat grids.

    ! Internal members, not read from input
    INTEGER  :: lon_dim             ! Number of points in lon direction
    INTEGER  :: lat_dim             ! Number of points in lat direction
    INTEGER  :: cur_bounds_triple   ! current output_bounds triple in use
    REAL(wp) :: next_output_time    ! next output time (in seconds simulation time)
    INTEGER  :: n_output_steps
    TYPE(t_output_name_list), POINTER :: next ! Pointer to next output_name_list

  END TYPE t_output_name_list

  ! Pointer to a linked list of output name lists:
  TYPE(t_output_name_list), POINTER :: first_output_name_list => NULL()

  !------------------------------------------------------------------------------------------------
  ! Max number of time levels:
  INTEGER, PARAMETER :: max_time_levels = 5

CONTAINS
  
  !-------------------------------------------------------------------------------------------------
  !>
  !! @return .TRUE. if one of the output namelists has been specified with GRIB output.

  FUNCTION is_grib_output() RESULT(retval)

    LOGICAL                           :: retval
    TYPE(t_output_name_list), POINTER :: p_onl

    retval = .FALSE.
    p_onl => first_output_name_list
    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT
      retval = retval                           .OR.  &
        &      (p_onl%filetype == FILETYPE_GRB) .OR.  &
        &      (p_onl%filetype == FILETYPE_GRB2)
      p_onl => p_onl%next
    ENDDO

  END FUNCTION is_grib_output

END MODULE mo_name_list_output_config
