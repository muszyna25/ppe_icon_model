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
MODULE mo_name_list_output_printvars

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char

  USE mo_cdi,                               ONLY: streamOpenWrite, FILETYPE_GRB2, gridCreate,                &
    &                                             GRID_UNSTRUCTURED, TAXIS_ABSOLUTE,                         &
    &                                             vlistDestroy, zaxisDestroy, streamClose, streamDefVlist,   &
    &                                             streamWriteVarSlice, gridDestroy, cdi_max_name,            &
    &                                             streamOpenRead, streamInqVlist, vlistinqvarname,           &
    &                                             taxisDestroy, vlistCreate, taxisCreate, institutInq,       &
    &                                             vlistDefTaxis, vlistDefInstitut
  USE mo_kind,                              ONLY: wp, dp
  USE mo_impl_constants,                    ONLY: ihs_ocean, SUCCESS
  USE mo_cf_convention,                     ONLY: t_cf_var
  USE mo_exception,                         ONLY: finish, message_text
  USE mo_linked_list,                       ONLY: t_var_list, t_list_element
  USE mo_var_metadata_types,                ONLY: t_var_metadata, VARNAME_LEN
  USE mo_gribout_config,                    ONLY: t_gribout_config
  USE mo_name_list_output_zaxes_types,      ONLY: t_verticalAxisList, t_verticalAxis
  USE mo_level_selection_types,             ONLY: t_level_selection
  USE mo_var_list_element,                  ONLY: t_var_list_element, level_type_ml
  USE mo_dictionary,                        ONLY: t_dictionary
  USE mo_util_sort,                         ONLY: quicksort
  USE mo_util_string,                       ONLY: remove_duplicates, toupper, tolower
  USE mo_util_file,                         ONLY: util_unlink
  USE mo_name_list_output_zaxes,            ONLY: setup_ml_axes_atmo, setup_pl_axis_atmo,         &
    &                                             setup_hl_axis_atmo, setup_il_axis_atmo,         &
    &                                             setup_zaxes_oce
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


CONTAINS

  !------------------------------------------------------------------------------------------------
  !> @return variable name without time level or tile suffix
  !
  FUNCTION get_var_basename(var)
    CHARACTER(LEN=VARNAME_LEN) :: get_var_basename
    TYPE(t_var_list_element)   :: var
    ! local variable
    INTEGER          :: endidx
    CHARACTER(LEN=1) :: suffix_str

    ! first cut off the time level suffix:
    endidx = INDEX(var%info%name,'.TL')
    IF (endidx==0) THEN
      endidx = LEN_TRIM(var%info%name)
    ELSE
      endidx = endidx - 1
    END IF

    suffix_str = " "

    ! condense two-digit suffices "01, 02, 03, ..." into "*"
    IF (endidx > 2) THEN
      IF (.NOT. is_number(var%info%name(endidx-2:endidx-2)) .AND. &
        &       is_number(var%info%name(endidx-1:endidx-1)) .AND. &
        &       is_number(var%info%name(endidx  :endidx  ))) THEN
        endidx = endidx - 2
        suffix_str = "*"
      END IF
    END IF

    ! condense one-digit suffices "_1, _2, _3, ..." into "_*"
    IF ((endidx > 1) .AND. (suffix_str == " ")) THEN
      IF ((var%info%name(endidx-1:endidx-1) == "_") .AND. &
        &  is_number(var%info%name(endidx  :endidx  ))) THEN
        endidx = endidx - 1
        suffix_str = "*"
      END IF
    END IF

    get_var_basename = TRIM(var%info%name(1:endidx))//suffix_str

  CONTAINS
    LOGICAL FUNCTION is_number(char)
      CHARACTER, INTENT(IN) :: char
      is_number = (IACHAR(char) - IACHAR('0')) <= 9
    END FUNCTION is_number

  END FUNCTION get_var_basename


  !------------------------------------------------------------------------------------------------
  !> @return container variable
  !
  FUNCTION get_var_container(var_lists, nvar_lists, print_patch_id, contained_elt)  RESULT(res)
    TYPE(t_list_element), POINTER :: res
    TYPE(t_var_list),              INTENT(IN) :: var_lists(:)
    INTEGER,                       INTENT(IN) :: nvar_lists
    INTEGER,                       INTENT(IN) :: print_patch_id
    TYPE(t_list_element), POINTER, INTENT(IN) :: contained_elt
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: i

    ! Unfortunately, there does not (yet) exist a link (pointer)
    ! between the contained element and the variable
    ! container. Therefore we need to loop over all variables and
    ! compare pointers.

    res => contained_elt
    VARLIST_LOOP : DO i = 1, nvar_lists
      IF (var_lists(i)%p%patch_id /= print_patch_id) CYCLE

      element => var_lists(i)%p%first_list_element
      DO
        IF (.NOT. ASSOCIATED(element)) EXIT
        IF (element%field%info%ncontained > 0) THEN
          IF (ASSOCIATED(element%field%r_ptr, contained_elt%field%r_ptr) .OR. &
            & ASSOCIATED(element%field%s_ptr, contained_elt%field%s_ptr) .OR. &
            & ASSOCIATED(element%field%i_ptr, contained_elt%field%i_ptr) .OR. &
            & ASSOCIATED(element%field%l_ptr, contained_elt%field%l_ptr)) THEN
            res => element
            EXIT VARLIST_LOOP
          END IF
        END IF

        element => element%next_list_element
      END DO
    END DO VARLIST_LOOP
  END FUNCTION get_var_container


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
    CHARACTER(LEN=*), PARAMETER :: tmp_filename = "tmpfile.grb"
    TYPE(t_verticalAxis),     POINTER :: zaxis
    INTEGER                           :: tmp_vlistID, tmp_gridID, tmp_zaxisID, tmp_varID,    &
      &                                  tmp_streamID, nmiss, i, tmp_taxisID, tmp_cdiInstID, &
      &                                  ierrstat
    REAL(wp)                          :: tmp_var1(1)

#ifndef HAVE_LIBGRIB_API
    name = ""
    RETURN
#endif

    ! retrieve vertical axis object
    zaxis => verticalAxisList%getEntry(icon_zaxis_type=info%vgrid)
    IF (.NOT. ASSOCIATED(zaxis)) THEN
      WRITE (message_text,'(a,i0,a)') 'Zaxis no. ', info%vgrid,' undefined.'
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
    tmp_streamID = streamOpenWrite(tmp_filename, FILETYPE_GRB2) 
    CALL streamDefVlist(tmp_streamID, tmp_vlistID) 
    tmp_var1(:) = 0._wp
    nmiss       = 0
    CALL streamWriteVarSlice(tmp_streamID, tmp_varID, 0, tmp_var1, nmiss) 
    CALL streamClose(tmp_streamID) 
    CALL vlistDestroy(tmp_vlistID)
    CALL taxisDestroy(tmp_taxisID)
    CALL gridDestroy(tmp_gridID) 
    ! Re-open file:
    tmp_streamID = streamOpenRead(tmp_filename) 
    ! Get the variable list of the dataset 
    tmp_vlistID = streamInqVlist(tmp_streamID) 
    CALL vlistInqVarName(tmp_vlistID, tmp_varID, name)   
    ! Close the input stream 
    CALL streamClose(tmp_streamID) 
    ierrstat = util_unlink(tmp_filename)
  END SUBROUTINE identify_grb2_shortname


  !------------------------------------------------------------------------------------------------
  !> Print list of all output variables (LaTeX table formatting).
  !
  SUBROUTINE print_var_list(var_lists, nvar_lists, out_varnames_dict,   &
    &                       print_patch_id, iequations, gribout_config, &
    &                       i_lctype)

    TYPE(t_var_list),       INTENT(IN) :: var_lists(:)
    INTEGER,                INTENT(IN) :: nvar_lists
    TYPE(t_dictionary),     INTENT(IN) :: out_varnames_dict
    INTEGER,                INTENT(IN) :: iequations
    TYPE(t_gribout_config), INTENT(IN) :: gribout_config
    INTEGER,                INTENT(IN) :: i_lctype
    INTEGER,                INTENT(IN) :: print_patch_id

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::print_var_list"
    INTEGER,          PARAMETER :: max_str_len = cdi_max_name + 1 + 128 + VARNAME_LEN + 99

    CHARACTER(LEN=*), PARAMETER :: varprefix = "\varname{"
    INTEGER,          PARAMETER :: PREF = LEN_TRIM(varprefix)

    TYPE(t_level_selection),  POINTER              :: tmp_level_selection => NULL()
    TYPE(t_verticalAxisList), TARGET               :: tmp_verticalAxisList
    TYPE(t_verticalAxisList), POINTER              :: it
    TYPE (t_var_metadata),    POINTER              :: info
    TYPE(t_list_element),     POINTER              :: element, src_element
    TYPE(t_cf_var),           POINTER              :: this_cf
    INTEGER                                        :: i, nout_vars, iout_var, ierrstat
    CHARACTER(kind=c_char, LEN = cdi_max_name + 1) :: name
    CHARACTER(LEN=max_str_len), ALLOCATABLE        :: out_vars(:)
    CHARACTER(len=128)                             :: descr_string
    ! ---------------------------------------------------------------------------

    ! generate the CDI IDs for vertical axes:
    IF (iequations/=ihs_ocean) THEN ! atm
      CALL setup_ml_axes_atmo(tmp_verticalAxisList, tmp_level_selection, print_patch_id)
#ifndef __NO_JSBACH__
      IF (ANY(echam_phy_config(:)%ljsb)) CALL setup_zaxes_jsbach(tmp_verticalAxisList)
#endif
    ELSE
      CALL setup_zaxes_oce(tmp_verticalAxisList, tmp_level_selection)
    END IF
    it => tmp_verticalAxisList
    DO
      IF (.NOT. ASSOCIATED(it%axis)) CALL finish(routine, "Internal error!")
      CALL it%axis%cdiZaxisCreate()
      IF (.NOT. ASSOCIATED(it%next))  EXIT
      it => it%next
    END DO

    ! count the no. of output variables:
    nout_vars = 0
    DO i = 1, nvar_lists
      IF (var_lists(i)%p%patch_id /= print_patch_id) CYCLE
      IF(.NOT. var_lists(i)%p%loutput) CYCLE
      IF (var_lists(i)%p%vlevel_type /= level_type_ml) CYCLE
      element => var_lists(i)%p%first_list_element
      DO
        IF (.NOT. ASSOCIATED(element)) EXIT
        info => element%field%info
        IF(.NOT. info%loutput) THEN
          element => element%next_list_element
          CYCLE
        END IF
        nout_vars = nout_vars + 1
        element => element%next_list_element
      END DO
    END DO

    ! allocate sufficient space
    ALLOCATE(out_vars(nout_vars), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    iout_var = 0
    DO i = 1, nvar_lists

      IF (var_lists(i)%p%patch_id /= print_patch_id) CYCLE
      ! Inspect only model level variables
      IF (var_lists(i)%p%vlevel_type /= level_type_ml) CYCLE
      ! Do not inspect element if output is disabled
      IF(.NOT. var_lists(i)%p%loutput) CYCLE

      element => var_lists(i)%p%first_list_element
      DO
        IF (.NOT. ASSOCIATED(element)) EXIT

        info => element%field%info

        ! Do not inspect element if output is disabled
        IF(.NOT. info%loutput) THEN
          element => element%next_list_element
          CYCLE
        END IF
        
        IF (info%post_op%lnew_cf) THEN
          this_cf => info%post_op%new_cf
        ELSE
          this_cf => info%cf
        END IF
        ! if no short is available and if the variable is a
        ! "reference" into another variable, then search for this
        ! source variable:
        IF ((LEN_TRIM(this_cf%long_name) == 0) .AND. info%lcontained) THEN
          src_element => get_var_container(var_lists, nvar_lists, print_patch_id, element)
          this_cf => src_element%field%info%cf
        END IF

        CALL identify_grb2_shortname(info, tmp_verticalAxisList,     &
          &                          gribout_config, i_lctype,       &
          &                          out_varnames_dict, name)
        name = TRIM(tolower(name))
        IF ((name(1:3) == "var") .OR. (name(1:5) == "param")) THEN
          name = ""
        ELSE
          name = "\varname{"//TRIM(name)//"}"
        END IF

        descr_string = this_cf%long_name
        ! upcase first letter of description string:
        descr_string(1:1) = toupper(descr_string(1:1))

        iout_var = iout_var + 1
        WRITE (out_vars(iout_var),'(5a)')                                             &
          &         "\varname{"//TRIM(tolower(get_var_basename(element%field)))//"}", &
          &         ' & ',                                                            &
          &         TRIM(name),                                                       &
          &         ' & ',                                                            &
          &         TRIM(descr_string)
        
        element => element%next_list_element
      ENDDO
    ENDDO
    CALL tmp_verticalAxisList%finalize()

    ! sort and remove duplicates
    CALL remove_duplicates(out_vars(1:nout_vars), nout_vars)
    CALL quicksort(out_vars(1:nout_vars))

    ! print table, but add a gap when new alphabetical letter starts:
    WRITE (0,*) "----------------------------------------------------------------------"
    WRITE (0,*) "List of output variables"
    WRITE (0,*) "----------------------------------------------------------------------"
    WRITE (0,*) " "
    DO i=1,nout_vars
      IF (i < nout_vars) THEN
        ! check for the initial character of the variable name:
        IF (out_vars(i)(PREF+1:PREF+1) /= out_vars(i+1)(PREF+1:PREF+1)) THEN
          WRITE (0,*) TRIM(out_vars(i)), " \\[0.5em]"
        ELSE
          WRITE (0,*) TRIM(out_vars(i)), " \\"
        END IF
      ELSE
        WRITE (0,*) TRIM(out_vars(i)), " \\"
      END IF
    END DO
    WRITE (0,*) " "

    DEALLOCATE(out_vars, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE print_var_list

END MODULE mo_name_list_output_printvars
