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
    &                              associate_keyword, with_keywords, &
    &                              int2string
  USE mo_exception,          ONLY: finish


  IMPLICIT NONE

  PRIVATE

  ! variables
  PUBLIC :: itopo
  PUBLIC :: fac_smooth_topo
  PUBLIC :: n_iter_smooth_topo
  PUBLIC :: hgtdiff_max_smooth_topo
  PUBLIC :: l_emiss
  PUBLIC :: read_nc_via_cdi
  PUBLIC :: heightdiff_threshold
  PUBLIC :: lrevert_sea_height
  PUBLIC :: itype_vegetation_cycle
  PUBLIC :: extpar_filename
  PUBLIC :: extpar_varnames_map_file
  PUBLIC :: i_lctype
  PUBLIC :: nclass_lu
  PUBLIC :: nmonths_ext

  ! subroutines/functions
  PUBLIC :: generate_filename
  PUBLIC :: generate_td_filename

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_extpar_config'


  !>
  !!----------------------------------------------------------------------------
  !! Derived type containing control variables specific to the nonhydrostatic 
  !! atm model

  !------------------------------------------------------------------------

  ! namelist variables

  INTEGER  :: itopo       ! 0: topography specified by analytical functions,
                          ! 1: topography read from netcdf files

  REAL(wp) :: fac_smooth_topo
  INTEGER  :: n_iter_smooth_topo(max_dom)
  REAL(wp) :: hgtdiff_max_smooth_topo(max_dom)
  LOGICAL  :: l_emiss     ! if true: read external emissivity map 
  LOGICAL  :: read_nc_via_cdi ! read netcdf input via cdi library (alternative: parallel netcdf)
  REAL(wp) :: heightdiff_threshold(max_dom)
  LOGICAL  :: lrevert_sea_height  ! if true: bring sea points back to original height
  INTEGER  :: itype_vegetation_cycle

  ! ExtPar input filename, may contain keywords, by default
  ! extpar_filename = "<path>extpar_<gridfile>"
  CHARACTER(LEN=filename_max) :: extpar_filename

  ! external parameter: dictionary which maps internal variable names
  ! onto GRIB2 shortnames or NetCDF var names.
  CHARACTER(LEN=filename_max) :: extpar_varnames_map_file

  ! information read from extpar file (no namelist parameter)
  !
  INTEGER ::  &           !< stores the landcover classification used for the external parameter data
    &  i_lctype(max_dom)  !< 1: Globcover2009, 2: GLC2000
                          !< defined in mo_ext_data_state:inquire_extpar_file

  INTEGER ::  &           !< number of landuse classes
    &  nclass_lu(max_dom)

  INTEGER ::  &           !< number of months in external data file
    &  nmonths_ext(max_dom)

  !!----------------------------------------------------------------------------

CONTAINS

  FUNCTION generate_filename(extpar_filename, model_base_dir, grid_filename, &
    &                        nroot, jlev, idom) &
    &  RESULT(result_str)
    CHARACTER(len=*), INTENT(IN)    :: extpar_filename, &
      &                                model_base_dir,  &
      &                                grid_filename
    INTEGER,          INTENT(IN)   :: nroot, jlev, idom
    CHARACTER(len=MAX_CHAR_LENGTH)  :: result_str
    TYPE (t_keyword_list), POINTER  :: keywords => NULL()

    CALL associate_keyword("<path>",     TRIM(model_base_dir), keywords)
    CALL associate_keyword("<gridfile>", TRIM(grid_filename),  keywords)
    CALL associate_keyword("<nroot>",  TRIM(int2string(nroot,"(i0)")),   keywords)
    CALL associate_keyword("<nroot0>", TRIM(int2string(nroot,"(i2.2)")), keywords)
    CALL associate_keyword("<jlev>",   TRIM(int2string(jlev, "(i2.2)")), keywords)
    CALL associate_keyword("<idom>",   TRIM(int2string(idom, "(i2.2)")), keywords)
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
    &  routine = modname//':generate_td_filename:'

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
