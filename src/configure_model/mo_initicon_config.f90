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
  USE mo_util_string,        ONLY: t_keyword_list, MAX_STRING_LEN,   &
    &                              associate_keyword, with_keywords, &
    &                              int2string
  USE mo_io_units,           ONLY: filename_max
  USE mo_impl_constants,     ONLY: max_dom

  IMPLICIT NONE

  PUBLIC :: init_mode, nlev_in, nlevsoil_in, zpbl1, zpbl2
  PUBLIC :: l_hice_in, l_sst_in     
  PUBLIC :: l_coarse2fine_mode
  PUBLIC :: ifs2icon_filename
  PUBLIC :: dwdfg_filename
  PUBLIC :: dwdana_filename
  PUBLIC :: generate_filename
  PUBLIC :: filetype
  PUBLIC :: ana_varnames_map_file


  CHARACTER(len=*),PARAMETER,PRIVATE :: &
    &  version = '$Id: mo_initicon_config.f90 10934 2013-01-14 08:50:18Z dreinert $'

  ! ----------------------------------------------------------------------------
  ! 1.0 Namelist variables for the prep_icon preprocessing program
  ! ----------------------------------------------------------------------------
  !
  INTEGER  :: init_mode     ! initialization mode
  INTEGER  :: nlev_in       ! number of model levels of input data
  INTEGER  :: nlevsoil_in   ! number of soil levels of input data

  REAL(wp) :: zpbl1, zpbl2  ! AGL heights used for vertical gradient computation
  LOGICAL  :: l_hice_in     ! Logical switch, if sea-ice thickness field is provided as input
  LOGICAL  :: l_sst_in      ! logical switch, if sea surface temperature is provided as input

  LOGICAL  :: l_coarse2fine_mode(max_dom)  ! If true, apply special corrections for interpolation from coarse
                                           ! to fine resolutions over mountainous terrain

  INTEGER  :: filetype      ! One of CDI's FILETYPE\_XXX constants. Possible values: 2 (=FILETYPE\_GRB2), 4 (=FILETYPE\_NC2)

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

CONTAINS
  
  FUNCTION generate_filename(input_filename, model_base_dir, &
    &                        nroot, jlev, idom)  RESULT(result_str)
    CHARACTER(len=*), INTENT(IN)   :: input_filename, &
      &                               model_base_dir
    INTEGER,          INTENT(IN)   :: nroot, jlev, idom
    CHARACTER(len=MAX_STRING_LEN)  :: result_str
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
