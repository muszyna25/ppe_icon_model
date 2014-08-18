!>
!! @author G. Zaengl
!!
!! @par Revision History
!! Moved configure state from namelists/mo_initicon_nml:
!! F. Prill, DWD (2012-01-31)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_initicon_config

  USE mo_kind,               ONLY: wp
  USE mo_util_string,        ONLY: t_keyword_list, associate_keyword, with_keywords, &
    &                              int2string
  USE mo_io_units,           ONLY: filename_max
  USE mo_impl_constants,     ONLY: max_dom, vname_len, max_var_ml, MAX_CHAR_LENGTH,  &
    &                              MODE_IFSANA, MODE_COMBINED, MODE_COSMODE,         &
    &                              MODE_DWDANA_INC, MODE_IAU
  USE mtime,                 ONLY: timedelta, newTimedelta, deallocateTimedelta,     &
    &                              max_timedelta_str_len
  USE mo_mtime_extensions,   ONLY: getPTStringFromMS

  IMPLICIT NONE

  PRIVATE 


  PUBLIC :: configure_initicon

  PUBLIC :: init_mode, nlev_in, nlevsoil_in, zpbl1, zpbl2
  PUBLIC :: dt_iau
  PUBLIC :: type_iau_wgt
  PUBLIC :: l_sst_in
  PUBLIC :: lread_ana     
  PUBLIC :: l_coarse2fine_mode
  PUBLIC :: lp2cintp_incr
  PUBLIC :: ifs2icon_filename
  PUBLIC :: dwdfg_filename
  PUBLIC :: dwdana_filename
  PUBLIC :: generate_filename
  PUBLIC :: ana_varlist
  PUBLIC :: filetype
  PUBLIC :: ana_varnames_map_file
  PUBLIC :: latbc_varnames_map_file
  PUBLIC :: init_mode_soil
  PUBLIC :: is_iau_active
  PUBLIC :: iau_wgt_dyn, iau_wgt_adv
  PUBLIC :: rho_incr_filter_wgt
  PUBLIC :: t_timeshift
  PUBLIC :: timeshift

  TYPE t_timeshift
    REAL(wp)                 :: dt_shift
    TYPE(timedelta)          :: mtime_shift
  END TYPE t_timeshift

  ! ----------------------------------------------------------------------------
  ! 1.0 Namelist variables for the init_icon preprocessing program
  ! ----------------------------------------------------------------------------
  !
  INTEGER  :: init_mode     ! initialization mode
  INTEGER  :: nlevsoil_in   ! number of soil levels of input data

  REAL(wp) :: zpbl1, zpbl2  ! AGL heights used for vertical gradient computation
  LOGICAL  :: l_sst_in      ! logical switch, if sea surface temperature is provided as input
  LOGICAL  :: lread_ana     ! If .TRUE., read analysis fields are read from analysis file
                            ! dwdana_filename. If .FALSE., ICON is soleyly started 
                            ! from first guess fields.   

  LOGICAL  :: l_coarse2fine_mode(max_dom)  ! If true, apply special corrections for interpolation from coarse
                                           ! to fine resolutions over mountainous terrain
  LOGICAL  :: lp2cintp_incr(max_dom) ! If true, perform parent-to-child interpolation of atmospheric data
                                     ! assimilation increments

  INTEGER  :: filetype      ! One of CDI's FILETYPE\_XXX constants. Possible values: 2 (=FILETYPE\_GRB2), 4 (=FILETYPE\_NC2)

  REAL(wp) :: dt_iau        ! Time interval during which incremental analysis update (IAU) is performed [s]. 
                            ! Only required for init_mode=MODE_DWDANA_INC, MODE_IAU

  TYPE(t_timeshift) :: &    ! Allows IAU runs to start earlier than the nominal simulation start date 
    &  timeshift            ! without showing up in the output metadata

  INTEGER  :: type_iau_wgt  ! Type of weighting function for IAU.
                            ! 1: Top-hat
                            ! 2: SIN2
                            ! Only required for init_mode=MODE_DWDANA_INC, MODE_IAU
  REAL(wp) :: rho_incr_filter_wgt  ! Vertical filtering weight for density increments 
                                   ! Only applicable for init_mode=MODE_DWDANA_INC, MODE_IAU

  CHARACTER(LEN=vname_len) :: ana_varlist(max_var_ml) ! list of mandatory analysis fields. 
                                                      ! This list can include a subset or the 
                                                      ! entire set of default analysis fields.

  ! IFS2ICON input filename, may contain keywords, by default
  ! ifs2icon_filename = "<path>ifs2icon_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: ifs2icon_filename

  ! DWD-FG input filename, may contain keywords, by default
  ! dwdfg_filename = "<path>dwdFG_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: dwdfg_filename

  ! DWD-ANA input filename, may contain keywords, by default
  ! dwdana_filename = "<path>dwdana_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: dwdana_filename

  ! analysis file: dictionary which maps internal variable names onto
  ! GRIB2 shortnames or NetCDF var names.
  CHARACTER(LEN=filename_max) :: ana_varnames_map_file      

  ! analysis file: dictionary which maps internal variable names onto
  ! GRIB2 shortnames or NetCDF var names used in lateral boundary nudging.
  CHARACTER(LEN=filename_max) :: latbc_varnames_map_file  
   
  ! ----------------------------------------------------------------------------
  ! Derived variables / variables based on input file contents
  ! ----------------------------------------------------------------------------

  INTEGER :: nlev_in   = 0  !< number of model levels of input data

  INTEGER :: init_mode_soil     !< initialization mode of soil model (coldstart, warmstart, warmstart+IAU)

  LOGICAL :: is_iau_active = .FALSE.  !< determines whether IAU is active at current time

  REAL(wp):: iau_wgt_dyn = 0._wp    !< IAU weight for dynamics fields 
  REAL(wp):: iau_wgt_adv = 0._wp    !< IAU weight for tracer fields

CONTAINS

  !>
  !! setup additional initicon control variables
  !!
  !! Setup of additional initicon control variables depending on the 
  !! initicon-NAMELIST and potentially other namelists. This routine is 
  !! called, after all namelists have been read and a synoptic consistency 
  !! check has been done.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2013-07-11)
  !!
  SUBROUTINE configure_initicon
  !
    CHARACTER(len=max_timedelta_str_len) :: PTshift
    TYPE(timedelta), POINTER             :: mtime_shift_local
    !-----------------------------------------------------------------------


    IF ( ANY((/MODE_IFSANA,MODE_COMBINED,MODE_COSMODE/) == init_mode) ) THEN
      init_mode_soil = 1   ! full coldstart is executed
                           ! i.e. w_so_ice and h_snow are re-diagnosed
    ELSE IF ( ANY((/MODE_DWDANA_INC, MODE_IAU/) == init_mode) ) THEN
      init_mode_soil = 3  ! warmstart (within assimilation cycle) with analysis increments for h_snow
    ELSE
      init_mode_soil = 2  ! warmstart with full fields for h_snow from snow analysis
    ENDIF


    !
    ! transform timeshift to mtime-format
    !
    CALL getPTStringFromMS(INT(timeshift%dt_shift * 1000._wp), PTshift)

    !******************************************************* 
    ! can be removed, as soon, as this issue is fixed in libmtime(timedeltaToString)
    IF (TRIM(PTshift)=="-P00.000S") THEN
      PTshift = "-PT00.000S"
    ELSE IF (TRIM(PTshift)=="+P00.000S") THEN
      PTshift = "+PT00.000S"
    ENDIF 
    !********************************************************

    mtime_shift_local => newTimedelta(TRIM(PTshift))
    timeshift%mtime_shift = mtime_shift_local

    ! cleanup
    CALL deallocateTimedelta(mtime_shift_local)

  END SUBROUTINE configure_initicon



  FUNCTION generate_filename(input_filename, model_base_dir, &
    &                        nroot, jlev, idom)  RESULT(result_str)
    CHARACTER(len=*), INTENT(IN)   :: input_filename, &
      &                               model_base_dir
    INTEGER,          INTENT(IN)   :: nroot, jlev, idom
    CHARACTER(len=MAX_CHAR_LENGTH) :: result_str
    TYPE (t_keyword_list), POINTER :: keywords => NULL()

    CALL associate_keyword("<path>",   TRIM(model_base_dir),             keywords)
    CALL associate_keyword("<nroot>",  TRIM(int2string(nroot,"(i1)")),   keywords)
    CALL associate_keyword("<nroot0>", TRIM(int2string(nroot,"(i2.2)")), keywords)
    CALL associate_keyword("<jlev>",   TRIM(int2string(jlev, "(i2.2)")), keywords)
    CALL associate_keyword("<idom>",   TRIM(int2string(idom, "(i2.2)")), keywords)
    ! replace keywords in "input_filename", which is by default
    ! ifs2icon_filename = "<path>ifs2icon_R<nroot>B<jlev>_DOM<idom>.nc"
    result_str = TRIM(with_keywords(keywords, TRIM(input_filename)))

  END FUNCTION generate_filename

END MODULE mo_initicon_config
