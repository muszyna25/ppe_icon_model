!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_extpar_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_dom, MAX_CHAR_LENGTH
  USE mo_io_units,           ONLY: filename_max
  USE mo_util_string,        ONLY: t_keyword_list, &
    &                              associate_keyword, with_keywords
  USE mo_exception,          ONLY: finish


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: itopo, fac_smooth_topo, n_iter_smooth_topo, l_emiss, heightdiff_threshold
  PUBLIC :: extpar_filename, generate_filename, generate_td_filename
  PUBLIC :: extpar_varnames_map_file

  !>
  !!----------------------------------------------------------------------------
  !! Derived type containing control variables specific to the nonhydrostatic 
  !! atm model

  !------------------------------------------------------------------------
  INTEGER  :: itopo       ! 0: topography specified by analytical functions,
                          ! 1: topography read from netcdf files

  REAL(wp) :: fac_smooth_topo
  INTEGER  :: n_iter_smooth_topo(max_dom)
  LOGICAL  :: l_emiss     ! if true: read external emissivity map 
  REAL(wp) :: heightdiff_threshold(max_dom)

  ! ExtPar input filename, may contain keywords, by default
  ! extpar_filename = "<path>extpar_<gridfile>"
  CHARACTER(LEN=filename_max) :: extpar_filename

  ! external parameter: dictionary which maps internal variable names
  ! onto GRIB2 shortnames or NetCDF var names.
  CHARACTER(LEN=filename_max) :: extpar_varnames_map_file
  
  !!----------------------------------------------------------------------------

CONTAINS

  FUNCTION generate_filename(extpar_filename, model_base_dir, grid_filename) &
    &  RESULT(result_str)
    CHARACTER(len=*), INTENT(IN)    :: extpar_filename, &
      &                                model_base_dir,  &
      &                                grid_filename
    CHARACTER(len=MAX_CHAR_LENGTH)  :: result_str
    TYPE (t_keyword_list), POINTER  :: keywords => NULL()

    CALL associate_keyword("<path>",     TRIM(model_base_dir), keywords)
    CALL associate_keyword("<gridfile>", TRIM(grid_filename),  keywords)
    ! replace keywords in "extpar_filename", which is by default
    ! extpar_filename = "<path>extpar_<gridfile>"
    result_str = TRIM(with_keywords(keywords, TRIM(extpar_filename)))

  END FUNCTION generate_filename
!-----------------------------------------------------------------------
  FUNCTION generate_td_filename(extpar_td_filename, model_base_dir, grid_filename, month, year, clim) &
    &  RESULT(result_str)
    CHARACTER(len=*), INTENT(IN)    :: extpar_td_filename, &
      &                                model_base_dir,  &
      &                                grid_filename
    INTEGER, INTENT(IN)             :: month
    INTEGER, INTENT(IN), OPTIONAL   :: year
    LOGICAL, INTENT(IN), OPTIONAL   :: clim
    CHARACTER(len=MAX_CHAR_LENGTH)  :: syear,smonth, result_str
    TYPE (t_keyword_list), POINTER :: keywords => NULL()
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &  routine = 'mo_td_ext_data:generate_td_filename:'

    IF (PRESENT (year)) THEN
     WRITE(syear, '(i4.4)') year
    ELSEIF ( PRESENT(clim) .AND. clim) THEN
     syear="CLIM"
    ELSE
          CALL finish(TRIM(ROUTINE),&
            & 'Missing year for a non climatological run')     
    END IF
    WRITE(smonth,'(i2.2)') month

    CALL associate_keyword("<path>",     TRIM(model_base_dir), keywords)
    CALL associate_keyword("<gridfile>", TRIM(grid_filename),  keywords)
    CALL associate_keyword("<year>", TRIM(syear),  keywords)
    CALL associate_keyword("<month>", TRIM(smonth),  keywords)
    ! replace keywords in "extpar_filename", which is by default
    ! extpar_td_filename = "<path>extpar_<year>_<month>_<gridfile>"
    ! if clim ist present and clim=.TRUE., <year> ist subst. by "CLIM"
    result_str = TRIM(with_keywords(keywords, TRIM(extpar_td_filename)))

  END FUNCTION generate_td_filename
END MODULE mo_extpar_config
