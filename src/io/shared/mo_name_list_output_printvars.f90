!>
!! Auxiliary module: Print a list of all output variables.
!!
!! Initial implementation  06/2018 : F. Prill (DWD)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
#include "icon_contiguous_defines.h"
MODULE mo_name_list_output_printvars

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char

  USE mo_cdi,                               ONLY: streamOpenWrite, FILETYPE_GRB2, gridCreate,                &
    &                                             GRID_UNSTRUCTURED, TAXIS_ABSOLUTE,                         &
    &                                             vlistDestroy, streamClose, streamDefVlist,                 &
    &                                             streamWriteVarSlice, gridDestroy, cdi_max_name,            &
    &                                             streamOpenRead, streamInqVlist, vlistinqvarname,           &
    &                                             taxisDestroy, vlistCreate, taxisCreate, institutInq,       &
    &                                             vlistDefTaxis, vlistDefInstitut
  USE mo_mpi,                               ONLY: get_my_global_mpi_id
  USE mo_kind,                              ONLY: wp, dp
  USE mo_impl_constants,                    ONLY: ihs_ocean, SUCCESS
  USE mo_cf_convention,                     ONLY: t_cf_var
  USE mo_exception,                         ONLY: finish, message_text
  USE mo_linked_list,                       ONLY: t_var_list_intrinsic
  USE mo_var_metadata_types,                ONLY: t_var_metadata, VARNAME_LEN
  USE mo_gribout_config,                    ONLY: t_gribout_config
  USE mo_name_list_output_zaxes_types,      ONLY: t_verticalAxisList, t_verticalAxis
  USE mo_level_selection_types,             ONLY: t_level_selection
  USE mo_var_list,                          ONLY: get_var_container, &
    &                                             total_number_of_variables, &
    &                                             var_lists_apply
  USE mo_var_list_element,                  ONLY: t_var_list_element, level_type_ml
  USE mo_dictionary,                        ONLY: t_dictionary
  USE mo_util_sort,                         ONLY: quicksort
  USE mo_util_string,                       ONLY: remove_duplicates, toupper, tolower, int2string
  USE mo_util_file,                         ONLY: util_unlink
  USE mo_name_list_output_zaxes,            ONLY: setup_ml_axes_atmo,         &
    &                                             setup_zaxes_oce
  USE mo_name_list_output_types,            ONLY: var_list_search_out_patch_lev, &
    &                                             var_list_filter_output_patch_levtype, &
    &                                             var_filter_output
  USE mo_util_cdi,                          ONLY: create_cdi_variable
#ifndef __NO_JSBACH__
  USE mo_echam_phy_config,                  ONLY: echam_phy_config
  USE mo_jsb_vertical_axes,                 ONLY: setup_zaxes_jsbach
#endif


  IMPLICIT NONE

  PRIVATE

  ! subroutines
  PUBLIC :: print_var_list

  !------------------------------------------------------------------------------------------------

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output_printvars'

  INTEGER, PARAMETER :: max_str_len = cdi_max_name + 1 + 128 + VARNAME_LEN + 63

  !> type to hold aggregation of printable data for \a print_var_list
  !! output text is appended to out_vars, iout_var keeps track of used
  !! slots and otherwise mostly holds arguments to identify_grb2_shortname
  TYPE, EXTENDS(var_list_search_out_patch_lev) :: output_var_search_state
    INTEGER :: iout_var
    CHARACTER(LEN=max_str_len), ALLOCATABLE :: out_vars(:)
    TYPE(t_verticalAxisList) :: tmp_verticalAxisList
    TYPE(t_dictionary), POINTER :: out_varnames_dict
    TYPE(t_gribout_config), POINTER :: gribout_config
    INTEGER :: i_lctype
  END TYPE output_var_search_state


CONTAINS

  !------------------------------------------------------------------------------------------------
  !> @return variable name without time level or tile suffix
  !
  FUNCTION get_var_basename(field)
    CHARACTER(LEN=VARNAME_LEN) :: get_var_basename
    TYPE(t_var_list_element), INTENT(in) :: field
    ! local variable
    INTEGER          :: endidx
    CHARACTER(LEN=1) :: suffix_str

    ! first cut off the time level suffix:
    endidx = INDEX(field%info%name,'.TL')
    IF (endidx==0) THEN
      endidx = LEN_TRIM(field%info%name)
    ELSE
      endidx = endidx - 1
    END IF

    suffix_str = " "

    ! condense three-digit suffices "001, 002, 003, ..." into "*"
    ! for container variables
    IF (endidx > 3) THEN
      IF (field%info%lcontained .AND. &
          .NOT. is_number(field%info%name(endidx-3:endidx-3)) .AND. &
        &       is_number(field%info%name(endidx-2:endidx-2)) .AND. &
        &       is_number(field%info%name(endidx-1:endidx-1)) .AND. &
        &       is_number(field%info%name(endidx  :endidx  ))) THEN
        endidx = endidx - 3
        suffix_str = "*"
      END IF
    END IF

    ! condense two-digit suffices "_01, _02, _03, ..." into "*"
    IF (endidx > 2) THEN
      IF ((field%info%name(endidx-2:endidx-2) == "_")    .AND. &
        &  is_number(field%info%name(endidx-1:endidx-1)) .AND. &
        &  is_number(field%info%name(endidx  :endidx  ))) THEN
        endidx = endidx - 2
        suffix_str = "*"
      END IF
    END IF

    ! condense one-digit suffices "_1, _2, _3, ..." into "_*"
    IF ((endidx > 1) .AND. (suffix_str == " ")) THEN
      IF ((field%info%name(endidx-1:endidx-1) == "_") .AND. &
        &  is_number(field%info%name(endidx  :endidx  ))) THEN
        endidx = endidx - 1
        suffix_str = "*"
      END IF
    END IF

    get_var_basename = field%info%name(1:endidx)//suffix_str

  CONTAINS
    LOGICAL FUNCTION is_number(char)
      CHARACTER, INTENT(IN) :: char
      INTEGER :: ia
      ia = IACHAR(char)
      is_number = ia >= IACHAR('0') .AND. ia <= IACHAR('9')
    END FUNCTION is_number

  END FUNCTION get_var_basename



  !------------------------------------------------------------------------------------------------
  !> @return GRIB2 short name of a given variable.
  !
  !  In order to make sure that we really get the true GRIB2 shortName
  !  (which depends on the definition files used), we choose a rather
  !  awkward approach: A temporary file is created and read in
  !  immediately afterwards. This is slow!
  !
  SUBROUTINE identify_grb2_shortname(info, verticalAxisList, gribout_config, i_lctype, &
    &                                out_varnames_dict, name)

    TYPE (t_var_metadata),                      INTENT(IN)    :: info
    TYPE(t_verticalAxisList),                   INTENT(INOUT) :: verticalAxisList
    TYPE(t_gribout_config),                     INTENT(IN)    :: gribout_config
    INTEGER,                                    INTENT(IN)    :: i_lctype
    TYPE(t_dictionary),                         INTENT(IN)    :: out_varnames_dict
    CHARACTER(kind=c_char, LEN=cdi_max_name+1), INTENT(OUT)   :: name
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::identify_grb2_shortname"
    CHARACTER(LEN=*), PARAMETER :: tmp_filename_base = "tmpfile.grb"
    TYPE(t_verticalAxis),     POINTER        :: zaxis
    INTEGER                                  :: tmp_vlistID, tmp_gridID, tmp_zaxisID, tmp_varID,    &
      &                                         tmp_streamID, nmiss, tmp_taxisID, tmp_cdiInstID, &
      &                                         ierrstat
    REAL(wp)                                 :: tmp_var1(1)
    CHARACTER(LEN=LEN(tmp_filename_base)+32) :: tmp_filename

#ifndef HAVE_LIBGRIB_API
    name = ""
    RETURN
#endif

    ! retrieve vertical axis object
    zaxis => verticalAxisList%getEntry(icon_zaxis_type=info%vgrid)
    IF (.NOT. ASSOCIATED(zaxis)) THEN
      WRITE (message_text,*) "variable '", TRIM(info%name), "' :Zaxis no. ", info%vgrid, " undefined."
      CALL finish(routine, message_text)
    END IF
    tmp_zaxisID = zaxis%cdi_id
    ! create dummy variable list
    tmp_vlistID = vlistCreate()
    tmp_cdiInstID = institutInq(gribout_config%generatingCenter,          &
      &                         gribout_config%generatingSubcenter, '', '')
    ! define Institute
    CALL vlistDefInstitut(tmp_vlistID, tmp_cdiInstID)
    ! create a dummy time axis 
    tmp_taxisID = taxisCreate(TAXIS_ABSOLUTE) 
    ! assign the time axis to the variable list 
    CALL vlistDefTaxis(tmp_vlistID, tmp_taxisID) 
    ! create dummy grid (size=1):
    tmp_gridID  = gridCreate(GRID_UNSTRUCTURED, 1)
    ! add the variable in question to the vlist:
    tmp_varID = create_cdi_variable(tmp_vlistID, tmp_gridID, tmp_zaxisID,        &
      &                             info, 0._dp, FILETYPE_GRB2,                  &
      &                             gribout_config, i_lctype, out_varnames_dict)
    ! open temporary GRIB2 file, write variable:
    tmp_filename = tmp_filename_base//"."//int2string(get_my_global_mpi_id(),'(i0)')
    tmp_streamID = streamOpenWrite(TRIM(tmp_filename), FILETYPE_GRB2) 
    CALL streamDefVlist(tmp_streamID, tmp_vlistID) 
    tmp_var1(:) = 0._wp
    nmiss       = 0
    CALL streamWriteVarSlice(tmp_streamID, tmp_varID, 0, tmp_var1, nmiss) 
    CALL streamClose(tmp_streamID) 
    CALL vlistDestroy(tmp_vlistID)
    CALL taxisDestroy(tmp_taxisID)
    CALL gridDestroy(tmp_gridID) 
    ! Re-open file:
    tmp_streamID = streamOpenRead(TRIM(tmp_filename)) 
    ! Get the variable list of the dataset 
    tmp_vlistID = streamInqVlist(tmp_streamID) 
    CALL vlistInqVarName(tmp_vlistID, tmp_varID, name)   
    ! Close the input stream 
    CALL streamClose(tmp_streamID) 
    ierrstat = util_unlink(TRIM(tmp_filename))
  END SUBROUTINE identify_grb2_shortname


  !------------------------------------------------------------------------------------------------
  !> Print list of all output variables (LaTeX table formatting).
  !
  SUBROUTINE print_var_list(out_varnames_dict,   &
    &                       print_patch_id, iequations, gribout_config, &
    &                       i_lctype)

    TYPE(t_dictionary), TARGET, INTENT(IN) :: out_varnames_dict
    INTEGER,                INTENT(IN) :: iequations
    TYPE(t_gribout_config), TARGET, INTENT(IN) :: gribout_config
    INTEGER,                INTENT(IN) :: i_lctype
    INTEGER,                INTENT(IN) :: print_patch_id

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::print_var_list"

    CHARACTER(LEN=*), PARAMETER :: varprefix = "\varname{"
    INTEGER,          PARAMETER :: PREF = LEN_TRIM(varprefix)

    TYPE(output_var_search_state), TARGET :: search_state
    TYPE(t_level_selection),  POINTER              :: tmp_level_selection => NULL()
    TYPE(t_verticalAxisList), POINTER              :: it
    INTEGER                                        :: i, nout_vars, ierrstat
    ! ---------------------------------------------------------------------------

    search_state%out_varnames_dict => out_varnames_dict
    search_state%gribout_config => gribout_config
    search_state%patch_id = print_patch_id
    search_state%ilev_type = level_type_ml
    search_state%iout_var = 0
    search_state%i_lctype = i_lctype
    ! count the no. of output variables:
    nout_vars = total_number_of_variables(var_list_filter_output_patch_levtype,&
         search_state, var_filter_output)

    ! generate the CDI IDs for vertical axes:
    IF (iequations/=ihs_ocean) THEN ! atm
      CALL setup_ml_axes_atmo(search_state%tmp_verticalAxisList, &
        &                     tmp_level_selection, print_patch_id)
#ifndef __NO_JSBACH__
      IF (ANY(echam_phy_config(:)%ljsb)) &
           CALL setup_zaxes_jsbach(search_state%tmp_verticalAxisList)
#endif
    ELSE
      CALL setup_zaxes_oce(search_state%tmp_verticalAxisList, &
        &                  tmp_level_selection)
    END IF
    it => search_state%tmp_verticalAxisList
    IF (.NOT. ASSOCIATED(it%axis)) CALL finish(routine, "Internal error!")
    DO WHILE (ASSOCIATED(it))
      CALL it%axis%cdiZaxisCreate()
      it => it%next
    END DO

    ! allocate sufficient space
    ALLOCATE(search_state%out_vars(nout_vars), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    CALL var_lists_apply(output_var_print, search_state, &
      &                  var_list_filter_output_patch_levtype)
    CALL search_state%tmp_verticalAxisList%finalize()

    ! sort and remove duplicates
    CALL remove_duplicates(search_state%out_vars(1:nout_vars), nout_vars)
    CALL quicksort(search_state%out_vars(1:nout_vars))

    ! print table, but add a gap when new alphabetical letter starts:
    WRITE (0,*) "----------------------------------------------------------------------"
    WRITE (0,*) "List of output variables"
    WRITE (0,*) "----------------------------------------------------------------------"
    WRITE (0,*) " "
    DO i=1,nout_vars
      IF (i < nout_vars) THEN
        ! check for the initial character of the variable name:
        IF (   search_state%out_vars(i)(PREF+1:PREF+1) &
          & /= search_state%out_vars(i+1)(PREF+1:PREF+1)) THEN
          WRITE (0,*) TRIM(search_state%out_vars(i)), " \\[0.5em]"
        ELSE
          WRITE (0,*) TRIM(search_state%out_vars(i)), " \\"
        END IF
      ELSE
        WRITE (0,*) TRIM(search_state%out_vars(i)), " \\"
      END IF
    END DO
    WRITE (0,*) " "

    DEALLOCATE(search_state%out_vars, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE print_var_list

  !> add single variable and its grib2 shortname to lines of varname output
  SUBROUTINE output_var_print(field, state, var_list)
    TYPE(t_var_list_element), TARGET :: field
    CLASS(*), TARGET :: state
    TYPE(t_var_list_intrinsic), INTENT(in) :: var_list

    TYPE(t_var_list_element), POINTER :: container
    TYPE(t_cf_var), POINTER :: this_cf
    INTEGER :: iout_var
    CHARACTER(len=128) :: descr_string
    CHARACTER(kind=c_char, LEN = cdi_max_name + 1) :: name

    SELECT TYPE (state)
    TYPE is (output_var_search_state)
      ! Do not inspect element if output is disabled
      IF (.NOT. field%info%loutput) RETURN
      IF (field%info%post_op%lnew_cf) THEN
        this_cf => field%info%post_op%new_cf
      ELSE
        this_cf => field%info%cf
      END IF
      ! if no short is available and if the variable is a
      ! "reference" into another variable, then search for this
      ! source variable:
      IF ((LEN_TRIM(this_cf%long_name) == 0) .AND. field%info%lcontained) THEN
        container => get_var_container(state%patch_id, field)
        this_cf => container%info%cf
      END IF

      CALL identify_grb2_shortname(field%info,                        &
        &                          state%tmp_verticalAxisList, &
        &                          state%gribout_config, state%i_lctype, &
        &                          state%out_varnames_dict, name)
      name = tolower(name)
      IF ((name(1:3) == "var") .OR. (name(1:5) == "param")) THEN
        name = ""
      ELSE
        name = "\varname{"//TRIM(name)//"}"
      END IF

      descr_string = this_cf%long_name
      ! upcase first letter of description string:
      descr_string(1:1) = toupper(descr_string(1:1))

      iout_var = state%iout_var + 1
      state%iout_var = iout_var
      WRITE (state%out_vars(iout_var),'(6a)') &
        &         "\varname{", &
        &         tolower(get_var_basename(field)), &
        &         "} & ", &
        &         TRIM(name), &
        &         ' & ', &
        &         TRIM(descr_string)
    END SELECT
  END SUBROUTINE output_var_print

END MODULE mo_name_list_output_printvars
