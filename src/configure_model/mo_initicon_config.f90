!>
!! @author G. Zaengl
!!
!! @par Revision History
!! Moved configure state from namelists/mo_initicon_nml:
!! F. Prill, DWD (2012-01-31)
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
MODULE mo_initicon_config

  USE mo_kind,               ONLY: wp
  USE mo_util_string,        ONLY: t_keyword_list, associate_keyword, with_keywords, &
    &                              int2string
  USE mo_io_units,           ONLY: filename_max
  USE mo_impl_constants,     ONLY: max_dom, vname_len, max_var_ml, MAX_CHAR_LENGTH,  &
    &                              MODE_IFSANA, MODE_COMBINED, MODE_COSMODE, MODE_DWDANA_INC

  IMPLICIT NONE

  PRIVATE 

  PUBLIC :: configure_initicon

  PUBLIC :: init_mode, nlev_in, nlevsoil_in, zpbl1, zpbl2
  PUBLIC :: dt_iau
  PUBLIC :: type_iau_wgt
  PUBLIC :: l_sst_in
  PUBLIC :: lread_ana     
  PUBLIC :: l_coarse2fine_mode
  PUBLIC :: ifs2icon_filename
  PUBLIC :: dwdfg_filename
  PUBLIC :: dwdana_filename
  PUBLIC :: generate_filename
  PUBLIC :: ana_varlist
  PUBLIC :: filetype
  PUBLIC :: ana_varnames_map_file
  PUBLIC :: init_mode_soil
  PUBLIC :: is_iau_active
  PUBLIC :: iau_wgt_dyn, iau_wgt_adv


  CHARACTER(len=*),PARAMETER,PRIVATE :: &
    &  version = '$Id$'

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

  INTEGER  :: filetype      ! One of CDI's FILETYPE\_XXX constants. Possible values: 2 (=FILETYPE\_GRB2), 4 (=FILETYPE\_NC2)

  REAL(wp) :: dt_iau        ! Time interval during which incremental analysis update (IAU) is performed [s]. 
                            ! Only required for init_mode=MODE_DWDANA_INC
  INTEGER  :: type_iau_wgt  ! Type of weighting function for IAU.
                            ! 1: Top-hat
                            ! 2: SIN2
                            ! Only required for init_mode=MODE_DWDANA_INC

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
    !-----------------------------------------------------------------------


    IF ( ANY((/MODE_IFSANA,MODE_COMBINED,MODE_COSMODE/) == init_mode) ) THEN
      init_mode_soil = 1   ! full coldstart is executed
                           ! i.e. w_so_ice and h_snow are re-diagnosed
    ELSE IF (init_mode == MODE_DWDANA_INC) THEN
      init_mode_soil = 3  ! warmstart (within assimilation cycle) with analysis increments for h_snow
    ELSE
      init_mode_soil = 2  ! warmstart with full fields for h_snow from snow analysis
    ENDIF

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
