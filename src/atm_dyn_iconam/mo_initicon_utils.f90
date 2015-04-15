!>
!! This module contains the I/O routines for initicon
!!
!! @author Guenther Zaengl, DWD
!!
!!
!! @par Revision History
!! First version by Guenther Zaengl, DWD (2011-07-13)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_initicon_utils

  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: msg_level, iqv, iqc, iqi, iqr, iqs, check_uuid_gracefully
  USE mo_dynamics_config,     ONLY: nnow, nnow_rcf, nnew, nnew_rcf
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_diag, t_nh_prog
  USE mo_nonhydrostatic_config, ONLY: kstart_moist
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag, t_wtr_prog
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_initicon_types,      ONLY: t_initicon_state, alb_snow_var,                     &
                                    ana_varnames_dict, inventory_list_fg, inventory_list_ana
  USE mo_initicon_config,     ONLY: init_mode, nlev_in, nlevsoil_in, l_sst_in,          &
    &                               timeshift,                                          &
    &                               ana_varlist, ana_varnames_map_file, lread_ana,      &
    &                               lconsistency_checks, lp2cintp_incr, lp2cintp_sfcana
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, MODE_DWDANA, MODE_DWDANA_INC,      &
    &                               MODE_IAU, MODE_IAU_OLD, MODE_IFSANA,                &
    &                               MODE_COMBINED, MODE_COSMODE, MODE_ICONVREMAP, MODIS,&
    &                               min_rlcell_int, grf_bdywidth_c
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_radiation_config,    ONLY: albedo_type
  USE mo_physical_constants,  ONLY: tf_salt, tmelt
  USE mo_exception,           ONLY: message, finish, message_text, warning
  USE mo_grid_config,         ONLY: n_dom
  USE mo_mpi,                 ONLY: my_process_is_stdio, p_io, p_bcast, p_comm_work_test, p_comm_work
  USE mo_util_string,         ONLY: tolower, difference, add_to_list, one_of
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, ntiles_total, lseaice, llake, lmulti_snow,         &
    &                               isub_lake, frlnd_thrhld, frlake_thrhld, frsea_thrhld, nlev_snow
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_phyparam_soil,       ONLY: csalb_snow_min, csalb_snow_max, csalb_snow, crhosmin_ml, crhosmax_ml
  USE mo_nh_init_utils,       ONLY: hydro_adjust
  USE mo_seaice_nwp,          ONLY: frsi_min, seaice_coldinit_nwp
  USE mo_dictionary,          ONLY: dict_init, dict_finalize,                           &
    &                               dict_loadfile, dict_get, DICT_MAX_STRLEN, dict_resize
  USE mo_post_op,             ONLY: perform_post_op
  USE mo_var_metadata_types,  ONLY: t_var_metadata, POST_OP_NONE, VARNAME_LEN
  USE mo_linked_list,         ONLY: t_list_element
  USE mo_var_list,            ONLY: get_var_name, nvar_lists, var_lists, collect_group
  USE mo_var_list_element,    ONLY: level_type_ml
  USE mo_util_cdi_table,      ONLY: t_inventory_list, t_inventory_element, &
    &                               find_inventory_list_element
  USE mo_util_bool_table,     ONLY: init_bool_table, add_column, print_bool_table, &
    &                               t_bool_table
  USE mo_util_uuid,           ONLY: OPERATOR(==)
  USE mo_flake,               ONLY: flake_coldinit
  USE mo_time_config,         ONLY: time_config
  USE mtime,                  ONLY: newDatetime, datetime, OPERATOR(==), OPERATOR(+), &
    &                               deallocateDatetime
  USE mo_intp_data_strc,      ONLY: t_int_state, p_int_state
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell
  USE mo_statistics,          ONLY: time_avg

  IMPLICIT NONE


  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_initicon_utils'


  PUBLIC :: initicon_inverse_post_op
  PUBLIC :: check_input_validity
  PUBLIC :: create_input_groups
  PUBLIC :: copy_initicon2prog_atm
  PUBLIC :: copy_initicon2prog_sfc
  PUBLIC :: allocate_initicon
  PUBLIC :: allocate_extana_atm
  PUBLIC :: allocate_extana_sfc
  PUBLIC :: deallocate_initicon
  PUBLIC :: deallocate_extana_atm 
  PUBLIC :: deallocate_extana_sfc
  PUBLIC :: average_first_guess
  PUBLIC :: reinit_average_first_guess
  PUBLIC :: fill_tile_points

  CONTAINS




  !-------------
  !>
  !! SUBROUTINE initicon_inverse_post_op
  !! Perform inverse post_op on input field, if necessary 
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2013-07-05)
  !!
  !!
  SUBROUTINE initicon_inverse_post_op(varname, mapped_name, optvar_out2D, optvar_out3D)

    CHARACTER(len=*), INTENT(IN)      :: varname             !< var name of field to be read
    CHARACTER(LEN=DICT_MAX_STRLEN)    :: mapped_name         !< mapped input name
    REAL(wp), OPTIONAL, INTENT(INOUT) :: optvar_out2D(:,:)   !< 3D output field
    REAL(wp), OPTIONAL, INTENT(INOUT) :: optvar_out3D(:,:,:) !< 2D output field

    ! local variables
    INTEGER                         :: i              ! loop count
    TYPE(t_var_metadata), POINTER   :: info           ! variable metadata
    TYPE(t_list_element), POINTER   :: element
    CHARACTER(len=*), PARAMETER     :: routine = modname//':initicon_inverse_post_op'

    !-------------------------------------------------------------------------


    ! Check consistency of optional arguments
    !
    IF (PRESENT( optvar_out2D ) .AND. PRESENT( optvar_out3D )) THEN
      CALL finish(routine, 'Only one optional argument must be present')
    ELSE IF (.NOT.PRESENT( optvar_out2D ) .AND. .NOT.PRESENT( optvar_out3D )) THEN
      CALL finish(routine, 'One of 2 optional arguments must be present')
    ENDIF


    ! get metadata information for field to be read
    info => NULL()
    DO i = 1,nvar_lists
      ! loop only over model level variables
      IF (var_lists(i)%p%vlevel_type /= level_type_ml) CYCLE 

      element => NULL()
      DO
        IF(.NOT.ASSOCIATED(element)) THEN
          element => var_lists(i)%p%first_list_element
        ELSE
          element => element%next_list_element
        ENDIF
        IF(.NOT.ASSOCIATED(element)) EXIT

        ! Check for matching name (take care of suffix of
        ! time-dependent variables):
        IF (TRIM(varname) == TRIM(tolower(get_var_name(element%field)))) THEN
          info => element%field%info
          EXIT
        ENDIF
      END DO

      ! If info handle has been found, exit list loop
      IF (ASSOCIATED(info)) EXIT
    ENDDO

    IF (.NOT.ASSOCIATED(info)) THEN
      WRITE (message_text,'(a,a)') TRIM(varname), ' not found'
      CALL message('',message_text)
      CALL finish(routine, 'Varname does not match any of the ICON variable names')
    ENDIF

    ! perform post_op
    IF (info%post_op%ipost_op_type /= POST_OP_NONE) THEN
      IF(my_process_is_stdio() .AND. msg_level>10) THEN
        WRITE(message_text,'(a)') 'Inverse Post_op for: '//TRIM(mapped_name)
        CALL message(TRIM(routine), TRIM(message_text))
      ENDIF
      IF (PRESENT(optvar_out2D)) THEN
        CALL perform_post_op(info%post_op, optvar_out2D, opt_inverse=.TRUE.)
      ELSE IF (PRESENT(optvar_out3D)) THEN
        CALL perform_post_op(info%post_op, optvar_out3D, opt_inverse=.TRUE.)
      ENDIF
    ENDIF

  END SUBROUTINE initicon_inverse_post_op




  !>
  !! Check validity of input fields
  !!
  !! Check validity of input fields. So far, the following checks are performed:
  !!
  !! For First Guess:
  !! - Check validity of uuidOfHgrid: The uuidOfHGrid of the input fields must 
  !!   match the uuidOfHgrid of the horizontal grid file.
  !! - Check validity of first guess validity time: First guess validity time 
  !!   must comply with the model's initialization time (ini_datetime) minus 
  !!   dt_shift.
  !!
  !! For Analysis (increments)
  !! - Check validity of uuidOfHgrid: The uuidOfHGrid of the input fields must 
  !!   match the uuidOfHgrid of the horizontal grid file.
  !! - For MODE_DWDANA, MODE_COMBINED, and MODE_COSMODE check validity of 
  !!   analysis validity time:  The analysis field's validity time must match 
  !!   the model start time
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-07-28)
  !!
  SUBROUTINE check_input_validity(p_patch, inventory_list_fg, inventory_list_ana, &
    &                             grp_vars_fg, grp_vars_ana, ngrp_vars_fg,        &
    &                             ngrp_vars_ana)

    TYPE(t_patch)                , INTENT(IN) :: p_patch
    TYPE(t_inventory_list)       , INTENT(IN) :: inventory_list_fg
    TYPE(t_inventory_list)       , INTENT(IN) :: inventory_list_ana
    CHARACTER(LEN=VARNAME_LEN)   , INTENT(IN) :: grp_vars_fg(:)   ! vars (names) to be read from fg-file
    CHARACTER(LEN=VARNAME_LEN)   , INTENT(IN) :: grp_vars_ana(:)  ! vars (names) to be read from ana-file
    INTEGER                      , INTENT(IN) :: ngrp_vars_fg     ! number of fields in grp_vars_fg
    INTEGER                      , INTENT(IN) :: ngrp_vars_ana    ! number of fields in grp_vars_ana

    ! local
    TYPE(t_inventory_element), POINTER  :: this_list_element => NULL()
    INTEGER :: ivar                     ! loop counter
    LOGICAL :: lmatch_uuid
    LOGICAL :: lmatch_vtime

    TYPE(datetime),  POINTER :: mtime_inidatetime  ! INI-datetime in mtime format
    TYPE(datetime)           :: start_datetime     ! true start date (includes timeshift)
    !
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':check_input_validity'

  !-------------------------------------------------------------------


    ! initialization
    !
    lmatch_uuid  = .FALSE. 
    lmatch_vtime = .FALSE.

    ! get ini-datetime in mtime-format
    mtime_inidatetime => newDatetime(time_config%ini_datetime%year,       &
      &                              time_config%ini_datetime%month,      &
      &                              time_config%ini_datetime%day,        &
      &                              time_config%ini_datetime%hour,       &
      &                              time_config%ini_datetime%minute,     &
      &                              INT(time_config%ini_datetime%second),&
      &                              ms=0)

    ! add timeshift to INI-datetime to get true starting time
    start_datetime = mtime_inidatetime + timeshift%mtime_shift



    ! Loop over all required input fields and perform some sanity checks.
    ! This is done separately for the analysis- and first guess fields.

    !
    ! First guess
    !
    DO ivar = 1,ngrp_vars_fg
      ! find matching list element
      this_list_element => find_inventory_list_element (inventory_list_fg, &
        &                                              TRIM(grp_vars_fg(ivar)) )


      !************************!
      !  Check uuidOfHGrid     !
      !************************!
      lmatch_uuid = (this_list_element%field%uuidOfHGrid == p_patch%grid_uuid)

      IF (.NOT. lmatch_uuid) THEN
        WRITE(message_text,'(a)') 'Non-matching uuidOfHGrid for first guess field '&
          &                       //TRIM(grp_vars_fg(ivar))//'.'
        IF (check_uuid_gracefully) THEN
          CALL warning(routine, TRIM(message_text))
        ELSE
          CALL finish(routine, TRIM(message_text))
        END IF
      ENDIF


      !**************************************!
      !  Check Validity time of first guess  !
      !**************************************!

      ! check correctness of validity-time
      lmatch_vtime = (this_list_element%field%vdatetime == start_datetime)

!       write(0,*) "vdatetime, start_datetime: ", this_list_element%field%vdatetime, start_datetime
!       write(0,*) "mtime_inidatetime, mtime_shift: ", mtime_inidatetime, timeshift%mtime_shift

      IF (.NOT. lmatch_vtime) THEN
        WRITE(message_text,'(a)') 'Non-matching validity datetime for first guess field '&
          &                       //TRIM(grp_vars_fg(ivar))//'.'
        CALL finish(routine, TRIM(message_text))
      ENDIF

    ENDDO


    !
    ! Analysis (increments)
    !
    DO ivar = 1,ngrp_vars_ana
      ! find matching list element
      this_list_element => find_inventory_list_element (inventory_list_ana, &
        &                                              TRIM(grp_vars_ana(ivar)) )

      !************************!
      !  Check uuidOfHGrid     !
      !************************!
      lmatch_uuid = (this_list_element%field%uuidOfHGrid == p_patch%grid_uuid)

      IF (.NOT. lmatch_uuid) THEN
        WRITE(message_text,'(a)') 'Non-matching uuidOfHGrid for analysis field '&
          &                       //TRIM(grp_vars_ana(ivar))//'.'
        IF (check_uuid_gracefully) THEN
          CALL warning(routine, TRIM(message_text))
        ELSE
          CALL finish(routine, TRIM(message_text))
        END IF
      ENDIF


      !**************************************!
      !  Check Validity time of analysis     !
      !**************************************!
      SELECT CASE (init_mode)
      CASE(MODE_DWDANA,MODE_COMBINED,MODE_COSMODE)
        !
        ! analysis field's validity time must match the model start time
        !
        ! check correctness of validity-time
        lmatch_vtime = (this_list_element%field%vdatetime == start_datetime)

        ! write(0,*) "vname: ", TRIM(grp_vars_ana(ivar))
        ! write(0,*) "vdatetime, start_datetime: ", this_list_element%field%vdatetime, start_datetime
        ! write(0,*) "mtime_inidatetime, mtime_shift: ", mtime_inidatetime, timeshift%mtime_shift

        IF (.NOT. lmatch_vtime) THEN
          WRITE(message_text,'(a)') 'Non-matching validity datetime for analysis field '&
            &                       //TRIM(grp_vars_ana(ivar))//'.'
          CALL finish(routine, TRIM(message_text))
        ENDIF
      CASE default
        !
      END SELECT
    ENDDO


    ! cleanup
    CALL deallocateDatetime(mtime_inidatetime)

  END SUBROUTINE check_input_validity



  !-------------
  !>
  !! SUBROUTINE create_input_groups
  !! Generates groups 'grp_vars_fg' and 'grp_vars_ana', which contain those fields which 
  !! must be read from the FG- and ANA-File, respectively.
  !! Both groups are based on two of a bunch of available ICON-internal output groups, depending on 
  !! which input mode is used
  !! groups for MODE_DWD     : mode_dwd_fg_in, mode_dwd_ana_in
  !! groups for MODE_IAU     : mode_iau_fg_in, mode_iau_ana_in
  !! groups for MODE_IAU_OLD : mode_iau_old_fg_in, mode_iau_old_ana_in
  !! groups for MODE_COMBINED: mode_combined_in
  !! groups for MODE_COSMODE : mode_cosmode_in
  !!
  !! In a first step it is checked, whether the ANA-File contains all members of the group 'grp_vars_ana'.
  !! If a member is missing, it is checked (based on the group grp_vars_ana_mandatory) whether the 
  !! ANA-Field is mandatory or not. If the field is mandatory, i.e. if it is part of the group 
  !! grp_vars_ana_mandatory provided via Namelist, the model aborts. If it is not mandatory, the model 
  !! tries to fall back to the corresponding FG-Field. This is done as follows:
  !! The missing ANA-Field is removed from the group 'grp_vars_ana' and added to the group 
  !! 'in_grp_vars_fg'. A warning is issued, however the model does not abort. In a second step it is 
  !! checked, whether the FG-File contains all members of the group 'grp_vars_fg'. If this is not the 
  !! case, the model aborts.
  !! At the end, a table is printed that shows which variables are part of which input group, meaning 
  !! which field will be read from which file.
  !!
  !! Special case: lread_ana=.FALSE.  : In this case, ICON will be started from first guess fields only
  !!                                    The analysis group varlist is re-set to 0 accordingly.
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2013-07-08)
  !!
  !!
  SUBROUTINE create_input_groups(p_patch, grp_vars_fg, ngrp_vars_fg, grp_vars_ana, ngrp_vars_ana, &
    &                            grp_vars_fg_default, ngrp_vars_fg_default,              &
    &                            grp_vars_ana_default, ngrp_vars_ana_default,            &
    &                            init_mode)

    TYPE(t_patch)             , INTENT(IN)    :: p_patch          ! current patch
    CHARACTER(LEN=VARNAME_LEN), INTENT(INOUT) :: grp_vars_fg(:)   ! vars (names) to be read from fg-file
    CHARACTER(LEN=VARNAME_LEN), INTENT(INOUT) :: grp_vars_ana(:)  ! vars (names) to be read from ana-file
    CHARACTER(LEN=VARNAME_LEN), INTENT(INOUT) :: grp_vars_fg_default(:)   ! default vars fg-file
    CHARACTER(LEN=VARNAME_LEN), INTENT(INOUT) :: grp_vars_ana_default(:)  ! default vars ana-file
    INTEGER                   , INTENT(OUT)   :: ngrp_vars_fg     ! number of fields in grp_vars_fg
    INTEGER                   , INTENT(OUT)   :: ngrp_vars_ana    ! number of fields in grp_vars_ana
    INTEGER                   , INTENT(OUT)   :: ngrp_vars_fg_default  ! default number grp_vars_fg
    INTEGER                   , INTENT(OUT)   :: ngrp_vars_ana_default ! default number grp_vars_ana
    INTEGER                   , INTENT(IN)    :: init_mode        ! initialization mode

    ! local variables
    INTEGER                    :: jg                              ! patch id
    CHARACTER(LEN=VARNAME_LEN) :: grp_vars_anafile(200)           ! ana-file inventory group
    CHARACTER(LEN=VARNAME_LEN) :: grp_vars_fgfile(200)            ! fg-file inventory group
    CHARACTER(LEN=VARNAME_LEN) :: grp_name                        ! group name
    INTEGER :: ivar, mpi_comm
    INTEGER :: index, is_one_of

    CHARACTER(LEN=*), PARAMETER :: routine = modname//':create_input_groups'
    TYPE(t_bool_table) :: bool_table

    ! list of mandatory analysis fields (provided via Namelist)
    CHARACTER(LEN=VARNAME_LEN) :: grp_vars_ana_mandatory(SIZE(grp_vars_ana_default))
    INTEGER :: nvars_ana_mandatory

    ! additional list for LOG printout
    CHARACTER(LEN=VARNAME_LEN) :: grp_vars_fg_default_grib2(SIZE(grp_vars_fg_default))
    CHARACTER(LEN=VARNAME_LEN) :: grp_vars_ana_default_grib2(SIZE(grp_vars_ana_default))
    CHARACTER(LEN=VARNAME_LEN) :: grp_vars_fg_grib2(SIZE(grp_vars_fg))
    CHARACTER(LEN=VARNAME_LEN) :: grp_vars_ana_grib2(SIZE(grp_vars_ana))
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: ana_default_txt, ana_this_txt

    CHARACTER(LEN=200) :: buffer_miss_ana   ! buffer for names of missing mandatory analysis fields
    CHARACTER(LEN=200) :: buffer_miss_fg    ! buffer for names of missing mandatory first guess fields
    LOGICAL :: lmiss_ana                    ! True, if there are missing mandatory analysis fields
    LOGICAL :: lmiss_fg                     ! True, if there are missing mandatory first guess fields

    INTEGER :: nelement                     ! element counter
    TYPE(t_inventory_element), POINTER :: current_element  ! pointer to linked list element
    !-------------------------------------------------------------------------

    IF(my_process_is_stdio()) THEN


      ! Initialization
      jg = p_patch%id
      lmiss_ana = .FALSE.
      lmiss_fg  = .FALSE.
      nvars_ana_mandatory = 0

      !===================
      ! 1: Collect groups
      !====================

      SELECT CASE(init_mode)
        CASE(MODE_DWDANA, MODE_DWDANA_INC,MODE_ICONVREMAP)
          ! Collect group 'grp_vars_fg_default' from mode_dwd_fg_in
          !
          grp_name ='mode_dwd_fg_in' 
          CALL collect_group(TRIM(grp_name), grp_vars_fg_default, ngrp_vars_fg_default,    &
            &                loutputvars_only=.FALSE.,lremap_lonlat=.FALSE.)


          ! When starting from analysis increments, we need both 
          ! the full FG field and the analysis increment
          ! maybe we should create separate groups for MODE_DWDANA_INC
          IF (init_mode == MODE_DWDANA_INC) THEN
             CALL add_to_list(grp_vars_fg_default, ngrp_vars_fg_default, (/'qv'/) , 1)
          ENDIF


          ! Collect group 'grp_vars_ana_default' from mode_dwd_ana_in
          !
          grp_name ='mode_dwd_ana_in' 
          CALL collect_group(TRIM(grp_name), grp_vars_ana_default, ngrp_vars_ana_default,    &
            &                loutputvars_only=.FALSE.,lremap_lonlat=.FALSE.)

          ! initialize grp_vars_fg and grp_vars_ana which will be the groups that control 
          ! the reading stuff
          !
          IF (lread_ana) THEN
            ! initialize grp_vars_fg and grp_vars_ana with grp_vars_fg_default and grp_vars_ana_default

            grp_vars_fg (1:ngrp_vars_fg_default) = grp_vars_fg_default (1:ngrp_vars_fg_default)
            grp_vars_ana(1:ngrp_vars_ana_default)= grp_vars_ana_default(1:ngrp_vars_ana_default)
            ngrp_vars_fg  = ngrp_vars_fg_default
            ngrp_vars_ana = ngrp_vars_ana_default
          ELSE
            ! lump together grp_vars_fg_default and grp_vars_ana_default
            !
            ! grp_vars_fg = grp_vars_fg_default + grp_vars_ana_default
            ngrp_vars_fg = 0
            CALL add_to_list(grp_vars_fg, ngrp_vars_fg, grp_vars_fg_default(1:ngrp_vars_fg_default)  , &
              &              ngrp_vars_fg_default)
            CALL add_to_list(grp_vars_fg, ngrp_vars_fg, grp_vars_ana_default(1:ngrp_vars_ana_default), &
              &              ngrp_vars_ana_default)

            ! Remove fields 'u', 'v', 'temp', 'pres'
            CALL difference(grp_vars_fg, ngrp_vars_fg, (/'u   ','v   ','temp','pres'/), 4)

            ! grp_vars_ana = --
            ngrp_vars_ana = 0
          ENDIF

        CASE(MODE_IAU)
          ! Collect group 'grp_vars_fg_default' from mode_dwd_fg_in
          !
          grp_name ='mode_iau_fg_in' 
          CALL collect_group(TRIM(grp_name), grp_vars_fg_default, ngrp_vars_fg_default,    &
            &                loutputvars_only=.FALSE.,lremap_lonlat=.FALSE.)

          ! Collect group 'grp_vars_ana_default' from mode_dwd_ana_in
          !
          grp_name ='mode_iau_ana_in' 
          CALL collect_group(TRIM(grp_name), grp_vars_ana_default, ngrp_vars_ana_default,    &
            &                loutputvars_only=.FALSE.,lremap_lonlat=.FALSE.)

          ! initialize grp_vars_fg and grp_vars_ana which will be the groups that control 
          ! the reading stuff
          !
          IF (.NOT. (lp2cintp_incr(jg) .AND. lp2cintp_sfcana(jg)) ) THEN
            ! initialize grp_vars_fg and grp_vars_ana with grp_vars_fg_default and grp_vars_ana_default

            grp_vars_fg (1:ngrp_vars_fg_default) = grp_vars_fg_default (1:ngrp_vars_fg_default)
            grp_vars_ana(1:ngrp_vars_ana_default)= grp_vars_ana_default(1:ngrp_vars_ana_default)
            ngrp_vars_fg  = ngrp_vars_fg_default
            ngrp_vars_ana = ngrp_vars_ana_default
          ELSE
            ! lump together grp_vars_fg_default and grp_vars_ana_default
            !
            ! grp_vars_fg = grp_vars_fg_default + grp_vars_ana_default
            ngrp_vars_fg = 0
            CALL add_to_list(grp_vars_fg, ngrp_vars_fg, grp_vars_fg_default(1:ngrp_vars_fg_default)  , &
              &              ngrp_vars_fg_default)
            CALL add_to_list(grp_vars_fg, ngrp_vars_fg, grp_vars_ana_default(1:ngrp_vars_ana_default), &
              &              ngrp_vars_ana_default)

            ! Remove fields 'u', 'v', 'temp', 'pres'
            CALL difference(grp_vars_fg, ngrp_vars_fg, (/'u   ','v   ','temp','pres'/), 4)

            ! grp_vars_ana = --
            ngrp_vars_ana = 0
          ENDIF

        CASE(MODE_IAU_OLD)
          ! Collect group 'grp_vars_fg_default' from mode_iau_old_fg_in
          !
          grp_name ='mode_iau_old_fg_in' 
          CALL collect_group(TRIM(grp_name), grp_vars_fg_default, ngrp_vars_fg_default,    &
            &                loutputvars_only=.FALSE.,lremap_lonlat=.FALSE.)

          ! Collect group 'grp_vars_ana_default' from mode_iau_old_ana_in
          !
          grp_name ='mode_iau_old_ana_in' 
          CALL collect_group(TRIM(grp_name), grp_vars_ana_default, ngrp_vars_ana_default,    &
            &                loutputvars_only=.FALSE.,lremap_lonlat=.FALSE.)

          ! initialize grp_vars_fg and grp_vars_ana which will be the groups that control 
          ! the reading stuff
          !
          IF (.NOT. (lp2cintp_incr(jg) .AND. lp2cintp_sfcana(jg)) ) THEN
            ! initialize grp_vars_fg and grp_vars_ana with grp_vars_fg_default and grp_vars_ana_default

            grp_vars_fg (1:ngrp_vars_fg_default) = grp_vars_fg_default (1:ngrp_vars_fg_default)
            grp_vars_ana(1:ngrp_vars_ana_default)= grp_vars_ana_default(1:ngrp_vars_ana_default)
            ngrp_vars_fg  = ngrp_vars_fg_default
            ngrp_vars_ana = ngrp_vars_ana_default
          ELSE
            ! lump together grp_vars_fg_default and grp_vars_ana_default
            !
            ! grp_vars_fg = grp_vars_fg_default + grp_vars_ana_default
            ngrp_vars_fg = 0
            CALL add_to_list(grp_vars_fg, ngrp_vars_fg, grp_vars_fg_default(1:ngrp_vars_fg_default)  , &
              &              ngrp_vars_fg_default)
            CALL add_to_list(grp_vars_fg, ngrp_vars_fg, grp_vars_ana_default(1:ngrp_vars_ana_default), &
              &              ngrp_vars_ana_default)

            ! Remove fields 'u', 'v', 'temp', 'pres'
            CALL difference(grp_vars_fg, ngrp_vars_fg, (/'u   ','v   ','temp','pres'/), 4)

            ! grp_vars_ana = --
            ngrp_vars_ana = 0
          ENDIF

        CASE(MODE_COMBINED,MODE_COSMODE)

          IF (init_mode == MODE_COMBINED) THEN
            grp_name ='mode_combined_in' 
          ELSE
            grp_name ='mode_cosmode_in' 
          ENDIF

          ! Collect group 'grp_vars_fg_default' from grp_name
          !
          CALL collect_group(TRIM(grp_name), grp_vars_fg_default, ngrp_vars_fg_default,    &
            &                loutputvars_only=.FALSE.,lremap_lonlat=.FALSE.)

          ! remove W_SO from default list and replace it by SMI
          CALL difference (grp_vars_fg_default, ngrp_vars_fg_default, (/'w_so'/), 1)
          CALL add_to_list(grp_vars_fg_default, ngrp_vars_fg_default, (/'smi'/) , 1)

          ! initialize grp_vars_fg which will be the group that controls the reading stuff
          !
          ! initialize grp_vars_fg with grp_vars_fg_default
          grp_vars_fg (1:ngrp_vars_fg_default) = grp_vars_fg_default (1:ngrp_vars_fg_default)

          ngrp_vars_fg = ngrp_vars_fg_default

          ! no analysis group
          ! ngrp_vars_ana_[default] = --
          ngrp_vars_ana_default = 0
          ngrp_vars_ana         = 0

        CASE DEFAULT

      END SELECT



      !===============================================================================
      ! 2: generate list of mandatory analysis fields (i.e. for which no fall back to
      !    FG fields is allowed )
      !===============================================================================

      IF( lread_ana .AND. ana_varlist(1) /= ' ' ) THEN
        ! translate GRIB2 varname to internal netcdf varname
        ! If requested GRIB2 varname is not found in the dictionary 
        ! (i.e. due to typos) -> Model abort
        DO ivar=1,SIZE(ana_varlist)
          IF (ana_varlist(ivar) /= ' ') THEN
            nvars_ana_mandatory = nvars_ana_mandatory + 1
            ! translate GRIB2 -> NetCDF
            grp_vars_ana_mandatory(ivar) = TRIM(dict_get(ana_varnames_dict, ana_varlist(ivar), &
              &                            linverse=.TRUE.))
          ELSE
            EXIT
          ENDIF
        ENDDO
      END IF


      !========================================================
      ! 3: Generate file inventory lists for FG- and ANA-files
      !========================================================

      ! get ANA-file varnames from inventory list (surface fields, only)
      ! Translation to internal names has already been performed in
      ! complete_inventory_list
      !
      IF (lread_ana .AND. .NOT. (lp2cintp_incr(jg) .AND. lp2cintp_sfcana(jg)) ) THEN  ! skip, when starting from first guess, only
        nelement = 0
        current_element => inventory_list_ana(jg)%p%first_list_element
        DO WHILE (ASSOCIATED(current_element))
          nelement = nelement + 1
          grp_vars_anafile(nelement) = TRIM(current_element%field%name)
          current_element => current_element%next_list_element
        ENDDO
      ENDIF


      ! get FG-file varnames from inventory list
      ! Translation to internal names has already been performed in
      ! complete_inventory_list
      !
      nelement = 0
      current_element => inventory_list_fg(jg)%p%first_list_element
      DO WHILE (ASSOCIATED(current_element))
        nelement = nelement + 1
        grp_vars_fgfile(nelement) = TRIM(current_element%field%name)
        current_element => current_element%next_list_element
      ENDDO



      !======================================
      ! 4: Check for missing input fields
      !======================================

      ! Check, whether the ANA-file inventory list contains all required analysis fields.
      ! If not, check whether the missing field is mandatory. If so, issue an error and abort. If 
      ! the field is not mandatory, remove the corresponding variable name from the group 
      ! 'grp_vars_ana' and issue a warning. The missing field is added to the group 'grp_vars_fg' 
      ! and thus the model tries to read it from the FG-File as fall back.
      !
      IF (lread_ana .AND. .NOT. (lp2cintp_incr(jg) .AND. lp2cintp_sfcana(jg)) ) THEN
        DO ivar=1,ngrp_vars_ana_default
          index = one_of(TRIM(grp_vars_ana_default(ivar)),grp_vars_anafile(:))

          IF ( index == -1) THEN  ! variable not found
            ! Check whether this field is mandatory, or whether we may fall back to 
            ! the first guess
            is_one_of = one_of(TRIM(grp_vars_ana_default(ivar)),grp_vars_ana_mandatory(1:nvars_ana_mandatory))

            IF ( is_one_of == -1) THEN  ! analysis field is not mandatory
              ! fall back to first guess
              !
              WRITE(message_text,'(a)') 'Field '//TRIM(grp_vars_ana_default(ivar))//' not found in ANA-input file.'
              CALL message(routine, TRIM(message_text))
              WRITE(message_text,'(a)') 'Field '//TRIM(grp_vars_ana_default(ivar))//' will be read from FG-input, instead.'
              CALL message(routine, TRIM(message_text))

              ! remove missing field from analysis input-group grp_vars_ana
              CALL difference(grp_vars_ana, ngrp_vars_ana, grp_vars_ana_default(ivar:ivar), 1)

              ! add missing field to the FG-group grp_vars_fg
              CALL add_to_list(grp_vars_fg, ngrp_vars_fg, grp_vars_ana_default(ivar:ivar), 1)

            ELSE  ! analysis field is mandatory

              ! add missing field to buffer
              IF (.NOT. lmiss_ana) THEN
                buffer_miss_ana = TRIM(grp_vars_ana_default(ivar))//', '
                lmiss_ana = .TRUE.
              ELSE
                IF ((LEN_TRIM(buffer_miss_ana)+LEN_TRIM(grp_vars_ana_default(ivar))+2)<= LEN(buffer_miss_ana)) THEN
                  buffer_miss_ana = TRIM(buffer_miss_ana)//TRIM(grp_vars_ana_default(ivar))//', '
                ELSE
                  CYCLE
                ENDIF
              ENDIF

            ENDIF
          ENDIF
        ENDDO
      ENDIF  ! lread_ana  



      ! Check, whether the FG-file inventory group contains all fields which are needed for a 
      ! successfull model start. If not, then stop the model and issue an error. 
      DO ivar=1,ngrp_vars_fg
        index = one_of(TRIM(grp_vars_fg(ivar)),grp_vars_fgfile(:))

        IF ( index == -1) THEN   ! variable not found

          ! add missing field to buffer
          IF (.NOT. lmiss_fg) THEN
            buffer_miss_fg = TRIM(grp_vars_fg(ivar))//', '
            lmiss_fg = .TRUE.
          ELSE
            IF ((LEN_TRIM(buffer_miss_fg)+LEN_TRIM(grp_vars_fg(ivar))+2)<= LEN(buffer_miss_fg)) THEN
              buffer_miss_fg = TRIM(buffer_miss_fg)//TRIM(grp_vars_fg(ivar))//', '
            ELSE
              CYCLE
            ENDIF
          ENDIF

          ! remove missing field from first guess input-group grp_vars_fg
          CALL difference(grp_vars_fg, ngrp_vars_fg, grp_vars_fg(ivar:ivar), 1)
        ENDIF
      ENDDO




      !====================
      ! 5: Printout table
      !====================

      !
      ! For printout, translate variable names from Netcdf (internal) to GRIB2
      ! print both GRIB2 (internal)
      DO ivar = 1,ngrp_vars_fg_default
        grp_vars_fg_default_grib2(ivar) = TRIM(dict_get(ana_varnames_dict,        &
          &                               TRIM(grp_vars_fg_default(ivar)),        &
          &                               linverse=.FALSE.))//" ("//              &
          &                               TRIM(grp_vars_fg_default(ivar))//")"
      ENDDO
      DO ivar = 1,ngrp_vars_fg
        grp_vars_fg_grib2(ivar) = TRIM(dict_get(ana_varnames_dict,     &
          &                       TRIM(grp_vars_fg(ivar)),             &
          &                       linverse=.FALSE.))//" ("//           &
          &                       TRIM(grp_vars_fg(ivar))//")"
      ENDDO
      DO ivar = 1,ngrp_vars_ana_default
        grp_vars_ana_default_grib2(ivar) = TRIM(dict_get(ana_varnames_dict,       &
          &                                TRIM(grp_vars_ana_default(ivar)),      &
          &                                linverse=.FALSE.))//" ("//             &
          &                                TRIM(grp_vars_ana_default(ivar))//")"
      ENDDO
      DO ivar = 1,ngrp_vars_ana
        grp_vars_ana_grib2(ivar) = TRIM(dict_get(ana_varnames_dict,    &
          &                        TRIM(grp_vars_ana(ivar)),           &
          &                        linverse=.FALSE.))//" ("//          &
          &                        TRIM(grp_vars_ana(ivar))//")"
      ENDDO
      WRITE(message_text,'(a,i2)') 'INIT_MODE ', init_mode
      CALL message(message_text, 'Required input fields: Source of FG and ANA fields')
      CALL init_bool_table(bool_table)
      IF ((init_mode == MODE_DWDANA_INC) .OR. (init_mode == MODE_IAU) .OR. (init_mode == MODE_IAU_OLD) ) THEN
        ana_default_txt = "ANA_inc (expected)"
        ana_this_txt    = "ANA_inc (this run)"
      ELSE
        ana_default_txt = "ANA (expected)"
        ana_this_txt    = "ANA (this run)"
      ENDIF
      CALL add_column(bool_table, "FG (expected)", grp_vars_fg_default_grib2,  ngrp_vars_fg_default)
      CALL add_column(bool_table, "FG (this run)",           grp_vars_fg_grib2,          ngrp_vars_fg)
      CALL add_column(bool_table, TRIM(ana_default_txt),grp_vars_ana_default_grib2, ngrp_vars_ana_default)
      CALL add_column(bool_table, TRIM(ana_this_txt)   ,grp_vars_ana_grib2,         ngrp_vars_ana)
      CALL print_bool_table(bool_table)


      !
      ! abort, if any mandatory first guess or analysis field is missing
      !
      IF (lmiss_ana) THEN
        WRITE(message_text,'(a)') 'Field(s) '//TRIM(buffer_miss_ana)// &
          &                       ' mandatory, but not found in ANA-file.'
        CALL finish(routine, TRIM(message_text))
      ENDIF
      !
      IF (lmiss_fg) THEN

        WRITE(message_text,'(a)') 'Field(s) '//TRIM(buffer_miss_fg)// &
          &                       ' missing in FG-input file.'
        CALL finish(routine, TRIM(message_text))
      ENDIF  



      ! additional sanity checks for input fields
      !
      IF ( lconsistency_checks ) THEN
        CALL check_input_validity(p_patch, inventory_list_fg(jg), inventory_list_ana(jg), &
          &                       grp_vars_fg, grp_vars_ana, ngrp_vars_fg, ngrp_vars_ana)
      ENDIF

    ENDIF  ! my_process_is_stdio()


    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test 
    ELSE
      mpi_comm = p_comm_work
    ENDIF
    CALL p_bcast(grp_vars_fg,  p_io, mpi_comm)
    CALL p_bcast(ngrp_vars_fg, p_io, mpi_comm)
    CALL p_bcast(grp_vars_ana, p_io, mpi_comm)
    CALL p_bcast(ngrp_vars_ana,p_io, mpi_comm)

  END SUBROUTINE create_input_groups

  !>
  !! SUBROUTINE fill_tile_points
  !! Used in the case of a tile coldstart, i.e. initializing a run with tiles with first guess
  !! data not containing tiles (or only tile-averaged variables)
  !! Specifically, this routine fills sub-grid scale (previously nonexistent) land and water points
  !! with appropriate data from neighboring grid points where possible
  !!
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD (2015-01-16)
  !!
  !!
  SUBROUTINE fill_tile_points(p_patch, p_lnd_state, ext_data)

    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data),     INTENT(INOUT) :: ext_data(:)

    TYPE(t_lnd_prog),  POINTER :: lnd_prog
    TYPE(t_lnd_diag),  POINTER :: lnd_diag
    TYPE(t_wtr_prog),  POINTER :: wtr_prog

    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

    INTEGER :: jg, jb, jk, jc, jt, ji
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

    REAL(wp), DIMENSION(nproma,10) :: lpmask, fpmask, spmask
    REAL(wp), DIMENSION(nproma)    :: lpcount, fpcount, spcount
    REAL(wp) :: wgt(10), fr_lnd, fr_lk, fr_sea, th_notile

    !-------------------------------------------------------------------------

    th_notile = 0.5_wp

    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      lnd_prog => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
      lnd_diag => p_lnd_state(jg)%diag_lnd
      wtr_prog => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))

      iidx     => p_int_state(jg)%rbf_c2grad_idx
      iblk     => p_int_state(jg)%rbf_c2grad_blk

      i_startblk = p_patch(jg)%cells%start_block(grf_bdywidth_c+1)
      i_endblk   = p_patch(jg)%cells%end_block(min_rlcell_int)

      DO ji = 2, 10
        SELECT CASE (ji)
        CASE(2,5,8) ! direct neighbors
          wgt(ji) = 1._wp
        CASE DEFAULT ! indirect neighbors
          wgt(ji) = 0.5_wp
        END SELECT
      ENDDO

!$OMP PARALLEL
!$OMP DO PRIVATE(ji,jb,jk,jc,jt,i_startidx,i_endidx,lpmask,fpmask,spmask,lpcount,fpcount,spcount,fr_lnd,fr_lk,fr_sea)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell_int)

        lpmask(:,:) = 0._wp
        fpmask(:,:) = 0._wp
        spmask(:,:) = 0._wp

        lpcount(:) = 0._wp
        fpcount(:) = 0._wp
        spcount(:) = 0._wp

        ! Prepare control variables to determine for which grid points neighbor filling is possible
        DO jc = i_startidx, i_endidx

          fr_lnd = ext_data(jg)%atm%fr_land(jc,jb)
          IF (fr_lnd > frlnd_thrhld .AND. fr_lnd < th_notile) THEN
            lpmask(jc,1) = 1._wp
            DO ji = 2,10
              fr_lnd = ext_data(jg)%atm%fr_land(iidx(ji,jc,jb),iblk(ji,jc,jb))
              IF (fr_lnd > th_notile) THEN
                lpmask(jc,ji) = wgt(ji)
                lpcount(jc) = lpcount(jc)+wgt(ji)
              ENDIF
            ENDDO
            IF (lpcount(jc) == 0._wp) lpmask(jc,1) = 0._wp
          ENDIF

          fr_lk = ext_data(jg)%atm%fr_lake(jc,jb)
          IF (fr_lk > frlake_thrhld .AND. fr_lk < th_notile) THEN
            fpmask(jc,1) = 1._wp
            DO ji = 2,10
              fr_lk = ext_data(jg)%atm%fr_lake(iidx(ji,jc,jb),iblk(ji,jc,jb))
              IF (fr_lk > th_notile) THEN
                fpmask(jc,ji) = wgt(ji)
                fpcount(jc) = fpcount(jc)+wgt(ji)
              ENDIF
            ENDDO
            IF (fpcount(jc) == 0._wp) fpmask(jc,1) = 0._wp
          ENDIF

          fr_sea = 1._wp - (ext_data(jg)%atm%fr_lake(jc,jb)+ext_data(jg)%atm%fr_land(jc,jb))
          IF (fr_sea > frsea_thrhld .AND. fr_sea < th_notile) THEN
            spmask(jc,1) = 1._wp
            DO ji = 2,10
              fr_sea = 1._wp - (ext_data(jg)%atm%fr_lake(iidx(ji,jc,jb),iblk(ji,jc,jb)) + &
                                ext_data(jg)%atm%fr_land(iidx(ji,jc,jb),iblk(ji,jc,jb))   )
              IF (fr_sea > th_notile) THEN
                spmask(jc,ji) = wgt(ji)
                spcount(jc) = spcount(jc)+wgt(ji)
              ENDIF
            ENDDO
            IF (spcount(jc) == 0._wp) spmask(jc,1) = 0._wp
          ENDIF

        ENDDO

        ! Apply neighbor filling
        DO jc = i_startidx, i_endidx

          ! a) ocean points
          IF (spmask(jc,1) == 1._wp) THEN
            CALL ngb_search(lnd_diag%fr_seaice, iidx, iblk, spmask, spcount, jc, jb)
            CALL ngb_search(wtr_prog%t_ice, iidx, iblk, spmask, spcount, jc, jb)
            CALL ngb_search(wtr_prog%h_ice, iidx, iblk, spmask, spcount, jc, jb)
          ENDIF

          ! b) lake points
          IF (fpmask(jc,1) == 1._wp) THEN
            CALL ngb_search(wtr_prog%t_mnw_lk, iidx, iblk, fpmask, fpcount, jc, jb)
            CALL ngb_search(wtr_prog%t_wml_lk, iidx, iblk, fpmask, fpcount, jc, jb)
            CALL ngb_search(wtr_prog%h_ml_lk,  iidx, iblk, fpmask, fpcount, jc, jb)
            CALL ngb_search(wtr_prog%t_bot_lk, iidx, iblk, fpmask, fpcount, jc, jb)
            CALL ngb_search(wtr_prog%c_t_lk,   iidx, iblk, fpmask, fpcount, jc, jb)
            CALL ngb_search(wtr_prog%t_b1_lk,  iidx, iblk, fpmask, fpcount, jc, jb)
            CALL ngb_search(wtr_prog%h_b1_lk,  iidx, iblk, fpmask, fpcount, jc, jb)
          ENDIF
        ENDDO

        ! c) land points
        DO jt = 1, ntiles_total

          ! single-layer fields
          DO jc = i_startidx, i_endidx
            IF (lpmask(jc,1) == 1._wp) THEN
              CALL ngb_search(lnd_diag%freshsnow_t(:,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
              CALL ngb_search(lnd_prog%w_snow_t   (:,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
              CALL ngb_search(lnd_prog%w_i_t      (:,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
              CALL ngb_search(lnd_diag%h_snow_t   (:,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
              CALL ngb_search(lnd_prog%t_snow_t   (:,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
              CALL ngb_search(lnd_prog%rho_snow_t (:,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)

              CALL ngb_search(lnd_prog%t_so_t(:,nlev_soil+1,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
            ENDIF
          ENDDO

          ! soil fields
          DO jk = 1, nlev_soil
            DO jc = i_startidx, i_endidx
              IF (lpmask(jc,1) == 1._wp) THEN
                CALL ngb_search(lnd_prog%t_so_t(:,jk,:,jt),     iidx, iblk, lpmask, lpcount, jc, jb)
                CALL ngb_search(lnd_prog%w_so_t(:,jk,:,jt),     iidx, iblk, lpmask, lpcount, jc, jb)
                CALL ngb_search(lnd_prog%w_so_ice_t(:,jk,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
              ENDIF
            ENDDO
          ENDDO

          IF (lmulti_snow) THEN ! multi-layer snow fields
            DO jk = 1, nlev_snow
              DO jc = i_startidx, i_endidx
                IF (lpmask(jc,1) == 1._wp) THEN
                  CALL ngb_search(lnd_prog%t_snow_mult_t(:,jk,:,jt),   iidx, iblk, lpmask, lpcount, jc, jb)
                  CALL ngb_search(lnd_prog%rho_snow_mult_t(:,jk,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
                  CALL ngb_search(lnd_prog%wtot_snow_t(:,jk,:,jt),     iidx, iblk, lpmask, lpcount, jc, jb)
                  CALL ngb_search(lnd_prog%wliq_snow_t(:,jk,:,jt),     iidx, iblk, lpmask, lpcount, jc, jb)
                  CALL ngb_search(lnd_prog%dzh_snow_t(:,jk,:,jt),      iidx, iblk, lpmask, lpcount, jc, jb)

                  IF (jk == 1) CALL ngb_search(lnd_prog%t_snow_mult_t(:,nlev_snow+1,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
                ENDIF
              ENDDO
            ENDDO
          ENDIF

        ENDDO

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ENDDO  ! jg

    CONTAINS

    SUBROUTINE ngb_search (fld, iidx, iblk, mask, cnt, jc, jb)

      REAL(wp), INTENT(INOUT) :: fld(:,:)
      REAL(wp), INTENT(IN)    :: mask(:,:), cnt(:)
      INTEGER,  INTENT(IN)    :: iidx(:,:,:), iblk(:,:,:), jc, jb

      INTEGER :: ji

      fld(jc,jb) = 0._wp
      DO ji = 2, 10
        fld(jc,jb) = fld(jc,jb) + fld(iidx(ji,jc,jb),iblk(ji,jc,jb))*mask(jc,ji)
      ENDDO
      fld(jc,jb) = fld(jc,jb)/cnt(jc)

    END SUBROUTINE ngb_search

  END SUBROUTINE fill_tile_points


  !>
  !! SUBROUTINE copy_initicon2prog_atm
  !! Copies atmospheric fields interpolated by init_icon to the
  !! prognostic model state variables 
  !!
  !! Required input: initicon state
  !! Output is written on fields of NH state
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-28)
  !! Modification by Daniel Reinert, DWD (2012-12-19)
  !! - encapsulated surface specific part
  !!
  !!
  SUBROUTINE copy_initicon2prog_atm(p_patch, initicon, p_nh_state)

    TYPE(t_patch),          INTENT(IN) :: p_patch(:)
    TYPE(t_initicon_state), INTENT(IN) :: initicon(:)

    TYPE(t_nh_state),      INTENT(INOUT) :: p_nh_state(:)

    INTEGER :: jg, jb, jk, jc, je
    INTEGER :: nblks_c, npromz_c, nblks_e, npromz_e, nlen, nlev, nlevp1, ntl, ntlr

!$OMP PARALLEL PRIVATE(jg,nblks_c,npromz_c,nblks_e,npromz_e,nlev,nlevp1,ntl,ntlr)
    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      nblks_c   = p_patch(jg)%nblks_c
      npromz_c  = p_patch(jg)%npromz_c
      nblks_e   = p_patch(jg)%nblks_e
      npromz_e  = p_patch(jg)%npromz_e
      nlev      = p_patch(jg)%nlev
      nlevp1    = p_patch(jg)%nlevp1
      ntl       = nnow(jg)
      ntlr      = nnow_rcf(jg)

!$OMP DO PRIVATE(jb,jk,je,nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_e

        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF

        ! Wind speed
        DO jk = 1, nlev
          DO je = 1, nlen
            p_nh_state(jg)%prog(ntl)%vn(je,jk,jb) = initicon(jg)%atm%vn(je,jk,jb)
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jb,jk,jc,nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c

        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF

        ! 3D fields
        DO jk = 1, nlev
          DO jc = 1, nlen
            ! Dynamic prognostic variables on cell points
            p_nh_state(jg)%prog(ntl)%w(jc,jk,jb)       = initicon(jg)%atm%w(jc,jk,jb)
            p_nh_state(jg)%prog(ntl)%theta_v(jc,jk,jb) = initicon(jg)%atm%theta_v(jc,jk,jb)
            p_nh_state(jg)%prog(ntl)%exner(jc,jk,jb)   = initicon(jg)%atm%exner(jc,jk,jb)
            p_nh_state(jg)%prog(ntl)%rho(jc,jk,jb)     = initicon(jg)%atm%rho(jc,jk,jb)

            ! Water vapor
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqv) = initicon(jg)%atm%qv(jc,jk,jb)
          ENDDO
        ENDDO

        ! Cloud and precipitation hydrometeors - these are supposed to be zero in the region where
        ! moisture physics is turned off
        DO jk = kstart_moist(jg), nlev
          DO jc = 1, nlen
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqc) = initicon(jg)%atm%qc(jc,jk,jb)
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqi) = initicon(jg)%atm%qi(jc,jk,jb)
            IF ( iqr /= 0 ) THEN
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqr) = initicon(jg)%atm%qr(jc,jk,jb)
            END IF
            IF ( iqs /= 0 ) THEN
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqs) = initicon(jg)%atm%qs(jc,jk,jb)
            END IF
          ENDDO
        ENDDO

        ! w at surface level
        DO jc = 1, nlen
          p_nh_state(jg)%prog(ntl)%w(jc,nlevp1,jb)      = initicon(jg)%atm%w(jc,nlevp1,jb)
          p_nh_state(jg)%prog(nnew(jg))%w(jc,nlevp1,jb) = initicon(jg)%atm%w(jc,nlevp1,jb)
        ENDDO

        IF (init_mode == MODE_ICONVREMAP) THEN ! copy also TKE field
          DO jk = 1, nlevp1
            DO jc = 1, nlen
              p_nh_state(jg)%prog(ntlr)%tke(jc,jk,jb) = initicon(jg)%atm%tke(jc,jk,jb)
            ENDDO
          ENDDO
        ENDIF

      ENDDO  ! jb
!$OMP END DO NOWAIT

    ENDDO  ! jg
!$OMP END PARALLEL

    ! Finally, compute exact hydrostatic adjustment for thermodynamic fields
    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE
      ntl = nnow(jg)

      CALL hydro_adjust(p_patch(jg), p_nh_state(jg)%metrics,                                  &
                        p_nh_state(jg)%prog(ntl)%rho,     p_nh_state(jg)%prog(ntl)%exner,     &
                        p_nh_state(jg)%prog(ntl)%theta_v )

    ENDDO

  END SUBROUTINE copy_initicon2prog_atm




  !>
  !! SUBROUTINE average_first_guess
  !! Averages atmospheric variables needed as first guess for data assimilation 
  !!
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD (2014-11-24)
  !! Modification by Daniel Reinert, DWD (2014-12-17)
  !! - make use of time_avg function, which makes finalize-option obsolete
  !!
  !!
  SUBROUTINE average_first_guess(p_patch, p_int, p_diag, p_prog_dyn, p_prog)

    TYPE(t_patch),          INTENT(IN) :: p_patch
    TYPE(t_int_state),      INTENT(IN) :: p_int

    TYPE(t_nh_diag),     INTENT(INOUT) :: p_diag
    TYPE(t_nh_prog),     INTENT(INOUT) :: p_prog_dyn, p_prog

    INTEGER :: jb, jk, jc
    INTEGER :: nlev, rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx
    REAL(wp):: wgt                     ! time average weight

    CHARACTER(len=*), PARAMETER     :: routine = modname//':average_first_guess'
    !------------------------------------------------------------------------------

    CALL rbf_vec_interpol_cell(p_prog_dyn%vn, p_patch, p_int, p_diag%u, p_diag%v, &
                               opt_rlend=min_rlcell_int)

    nlev = p_patch%nlev

    rl_start = 1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)


    p_diag%nsteps_avg = p_diag%nsteps_avg + 1

    ! compute weight
    wgt = 1._wp/REAL(p_diag%nsteps_avg(1),wp)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          p_diag%u_avg(jc,jk,jb)    = time_avg(p_diag%u_avg(jc,jk,jb)   , &
            &                                  p_diag%u(jc,jk,jb)       , &
            &                                  wgt)
          p_diag%v_avg(jc,jk,jb)    = time_avg(p_diag%v_avg(jc,jk,jb)   , &
            &                                  p_diag%v(jc,jk,jb)       , &
            &                                  wgt)
          p_diag%temp_avg(jc,jk,jb) = time_avg(p_diag%temp_avg(jc,jk,jb), &
            &                                  p_diag%temp(jc,jk,jb)    , &
            &                                  wgt)
          p_diag%pres_avg(jc,jk,jb) = time_avg(p_diag%pres_avg(jc,jk,jb), &
            &                                  p_diag%pres(jc,jk,jb)    , &
            &                                  wgt)
          p_diag%qv_avg(jc,jk,jb)   = time_avg(p_diag%qv_avg(jc,jk,jb)  , &
            &                                  p_prog%tracer(jc,jk,jb,iqv), &
            &                                  wgt)
        ENDDO
      ENDDO
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL


    ! debug output
    IF(my_process_is_stdio() .AND. msg_level>=13) THEN
      WRITE(message_text,'(a,I3)') 'step ', p_diag%nsteps_avg(1)
      CALL message(TRIM(routine), TRIM(message_text))
    ENDIF

  END SUBROUTINE average_first_guess


  !>
  !! SUBROUTINE reinit_average_first_guess
  !! Re-Initialization routine for SUBROUTINE average_first_guess.
  !! Ensures that the average is centered in time. 
  !!
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD (2015-02-10)
  !!
  !!
  SUBROUTINE reinit_average_first_guess(p_patch, p_diag, p_prog)

    TYPE(t_patch),       INTENT(IN)    :: p_patch

    TYPE(t_nh_diag),     INTENT(INOUT) :: p_diag
    TYPE(t_nh_prog),     INTENT(IN)    :: p_prog

    INTEGER :: jb, jk, jc
    INTEGER :: nlev, rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx

    CHARACTER(len=*), PARAMETER     :: routine = modname//':reinit_average_first_guess'

    !------------------------------------------------------------------------------

    nlev = p_patch%nlev

    rl_start = 1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_diag%u_avg(jc,jk,jb)    = p_diag%u(jc,jk,jb)
          p_diag%v_avg(jc,jk,jb)    = p_diag%v(jc,jk,jb)
          p_diag%temp_avg(jc,jk,jb) = p_diag%temp(jc,jk,jb)
          p_diag%pres_avg(jc,jk,jb) = p_diag%pres(jc,jk,jb)
          p_diag%qv_avg(jc,jk,jb)   = p_prog%tracer(jc,jk,jb,iqv)
        ENDDO
      ENDDO
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

    p_diag%nsteps_avg = 1

    ! debug output
    IF(my_process_is_stdio() .AND. msg_level>=13) THEN
      WRITE(message_text,'(a,I3)') 'step ', p_diag%nsteps_avg(1)
      CALL message(TRIM(routine), TRIM(message_text))
    ENDIF

  END SUBROUTINE reinit_average_first_guess



  !-------------
  !>
  !! SUBROUTINE copy_initicon2prog_sfc
  !! Copies surface fields interpolated by init_icon to the prognostic model 
  !! state variables. 
  !!
  !! Required input: initicon state
  !! Output is written on fields of land state
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-28)
  !! Modification by Daniel Reinert, DWD (2012-12-19)
  !! - encapsulated surface specific part
  !! Modification by Daniel Reinert, DWD (2013-07-09)
  !! - moved sea-ice coldstart into separate subroutine
  !!
  !!
  SUBROUTINE copy_initicon2prog_sfc(p_patch, initicon, p_lnd_state, ext_data)

    TYPE(t_patch),          INTENT(IN) :: p_patch(:)
    TYPE(t_initicon_state), INTENT(IN) :: initicon(:)

    TYPE(t_lnd_state),     INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data), INTENT(   IN) :: ext_data(:)

    INTEGER  :: jg, jb, jc, jt, js, jp, ic
    INTEGER  :: nblks_c, npromz_c, nlen
    REAL(wp) :: zfrice_thrhld, zminsnow_alb, t_fac

!$OMP PARALLEL PRIVATE(jg,nblks_c,npromz_c)
    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      nblks_c   = p_patch(jg)%nblks_c
      npromz_c  = p_patch(jg)%npromz_c


!$OMP DO PRIVATE(jb,jc,nlen,jt,js,jp,ic,zfrice_thrhld,zminsnow_alb,t_fac) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c

        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF


        ! ground temperature
        DO jc = 1, nlen
          p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_g(jc,jb) = initicon(jg)%sfc%tskin(jc,jb)
          p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_g(jc,jb) = initicon(jg)%sfc%tskin(jc,jb)
        ENDDO

        ! Fill also SST and sea ice fraction fields over ocean points; SST is limited to 30 deg C
        ! Note: missing values of the sea ice fraction, which may occur due to differing land-sea masks, 
        ! are indicated with -999.9; non-ocean points are filled with zero for both fields
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, ext_data(jg)%atm%sp_count(jb)
          jc = ext_data(jg)%atm%idx_lst_sp(ic,jb)
          IF ( l_sst_in .AND. initicon(jg)%sfc%sst(jc,jb) > 10._wp  ) THEN
            p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = initicon(jg)%sfc%sst(jc,jb)              
          ELSE
           p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = MIN(303.15_wp,initicon(jg)%sfc%tskin(jc,jb))
          ENDIF
          !
          ! In case of missing sea ice fraction values, we make use of the sea 
          ! surface temperature (tskin over ocean points). For tskin<=tf_salt, 
          ! we set the sea ice fraction to one. For tskin>tf_salt, we set it to 0.
          ! Note: tf_salt=271.45K is the salt-water freezing point
          !
          IF ( initicon(jg)%sfc%seaice(jc,jb) > -999.0_wp ) THEN
            p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = initicon(jg)%sfc%seaice(jc,jb) 
          ELSE    ! missing value
            IF ( initicon(jg)%sfc%tskin(jc,jb) <= tf_salt ) THEN
              p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 1._wp     ! sea ice point
            ELSE
              p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 0._wp     ! water point
            ENDIF
          ENDIF

          ! For fr_seaice in ]0,frsi_min[, set fr_seaice to 0
          ! For fr_seaice in ]1-frsi_min,1[, set fr_seaice to 1. This will ensure in 
          ! init_sea_lists, that sea-ice and water fractions sum up exactly to the total 
          ! sea fraction.
          IF (p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) < frsi_min ) THEN
             p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 0._wp
          ENDIF
          IF (p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) > (1._wp-frsi_min) ) THEN
             p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 1._wp
          ENDIF

        ENDDO
        ! In addition, write skin temperature to lake points, limited to 33 deg C. These will
        ! be used to initialize lake points until something more reasonable becomes available
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, ext_data(jg)%atm%fp_count(jb)
          jc = ext_data(jg)%atm%idx_lst_fp(ic,jb)
          p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = MIN(306.15_wp,initicon(jg)%sfc%tskin(jc,jb))
        ENDDO

        IF ( atm_phy_nwp_config(jg)%inwp_surface > 0 ) THEN
          DO jt = 1, ntiles_total
            DO jc = 1, nlen
              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_snow_t(jc,jb,jt)           = &
                &                                                initicon(jg)%sfc%tsnow   (jc,jb)

               ! Initialize freshsnow
               ! for seapoints, freshsnow is set to 0
              IF(alb_snow_var == 'ALB_SNOW') THEN

                IF (albedo_type == MODIS) THEN
                  IF (ext_data(jg)%atm%alb_dif(jc,jb) > csalb_snow_min) THEN
                    t_fac = MIN(1._wp,0.1_wp*(tmelt-p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_snow_t(jc,jb,jt)))
                    zminsnow_alb = (1._wp-t_fac)*csalb_snow_min + t_fac*ext_data(jg)%atm%alb_dif(jc,jb)
                  ELSE
                    zminsnow_alb = csalb_snow_min
                  ENDIF
                ELSE
                  IF (ext_data(jg)%atm%lc_class_t(jc,jb,jt) == ext_data(jg)%atm%i_lc_snow_ice) THEN
                    t_fac = MIN(1._wp,0.1_wp*(tmelt-p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_snow_t(jc,jb,jt)))
                    zminsnow_alb = (1._wp-t_fac)*csalb_snow_min + t_fac*csalb_snow
                  ELSE
                    zminsnow_alb = csalb_snow_min
                  ENDIF
                ENDIF

                p_lnd_state(jg)%diag_lnd%freshsnow_t(jc,jb,jt)    =  MAX(0._wp,MIN(1._wp, &
            &                           (initicon(jg)%sfc%snowalb (jc,jb)-zminsnow_alb)   &
            &                          /(csalb_snow_max-zminsnow_alb)))                   &
            &                          * REAL(NINT(ext_data(jg)%atm%fr_land(jc,jb)),wp) 
              ELSE
                p_lnd_state(jg)%diag_lnd%freshsnow_t(jc,jb,jt)    =  MAX(0._wp,MIN(1._wp, &
            &                     1._wp - ((initicon(jg)%sfc%snowalb (jc,jb)-crhosmin_ml) &
            &                    /(crhosmax_ml-crhosmin_ml))))                            &
            &                    * REAL(NINT(ext_data(jg)%atm%fr_land(jc,jb)),wp)
              END IF


              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_snow_t(jc,jb,jt)           = &
                &                                                initicon(jg)%sfc%snowweq (jc,jb)
              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%rho_snow_t(jc,jb,jt)         = &
                &                                                initicon(jg)%sfc%snowdens(jc,jb) 
              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_i_t(jc,jb,jt)              = &
                &                                                initicon(jg)%sfc%skinres (jc,jb)

              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_snow_t(jc,jb,jt)   = &
                &                                                initicon(jg)%sfc%tsnow   (jc,jb)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%w_snow_t(jc,jb,jt)   = &
                &                                                initicon(jg)%sfc%snowweq (jc,jb)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%rho_snow_t(jc,jb,jt) = &
                &                                                initicon(jg)%sfc%snowdens(jc,jb)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%w_i_t(jc,jb,jt)      = &
                &                                                initicon(jg)%sfc%skinres (jc,jb)
            ENDDO
          ENDDO

          ! Multi-layer surface fields
          DO jt = 1, ntiles_total

            DO js = 0, nlev_soil
              jp = js+1 ! indexing for the ICON state field starts at 1
              DO jc = 1, nlen
                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(jc,jp,jb,jt)= &
                  &                                              initicon(jg)%sfc%tsoil(jc,js,jb)
                p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_so_t(jc,jp,jb,jt)= &
                  &                                              initicon(jg)%sfc%tsoil(jc,js,jb)
              ENDDO
            ENDDO

            ! For soil water, no comparable layer shift exists
            DO js = 1, nlev_soil
              DO jc = 1, nlen
                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_so_t(jc,js,jb,jt)= &
                  &                                              initicon(jg)%sfc%wsoil(jc,js,jb)
                p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%w_so_t(jc,js,jb,jt)= &
                  &                                              initicon(jg)%sfc%wsoil(jc,js,jb)
              ENDDO
            ENDDO

          ENDDO



          ! Coldstart for sea-ice parameterization scheme
          ! Sea-ice surface temperature is initialized with tskin from IFS.
          ! Since the seaice index list is not yet available at this stage, we loop over 
          ! all sea points and initialize points with fr_seaice>threshold. 
          ! The threshold is 0.5 without tiles and frsi_min with tiles.
          ! Note that exactly the same threshold values must be used as in init_sea_lists. 
          ! If not, you will see what you get.
          !@Pilar: This should still work out for you, since the non-sea-ice points are 
          !        now initialized during warmstart initialization 
          !        in mo_nwp_sfc_utils:nwp_surface_init
          !
          IF (lseaice) THEN
            IF ( ntiles_total == 1 ) THEN  ! no tile approach
              zfrice_thrhld = 0.5_wp
            ELSE
              zfrice_thrhld = frsi_min
            ENDIF

            CALL seaice_coldinit_nwp(nproma, zfrice_thrhld,                               &
              &         frsi    = p_lnd_state(jg)%diag_lnd%fr_seaice(:,jb),               &
              &         temp_in = initicon(jg)%sfc%tskin(:,jb),                           &
              &         tice_p  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_ice    (:,jb), &
              &         hice_p  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ice    (:,jb), &
              &         tsnow_p = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_snow_si(:,jb), &
              &         hsnow_p = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_snow_si(:,jb), &
              &         tice_n  = p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%t_ice    (:,jb), &
              &         hice_n  = p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%h_ice    (:,jb), &
              &         tsnow_n = p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%t_snow_si(:,jb), &
              &         hsnow_n = p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%h_snow_si(:,jb)  )

          ENDIF  ! leseaice

          ! Cold-start initialization of the fresh-water lake model FLake.
          ! The procedure is the same as in "int2lm".
          ! Note that no lake ice is assumed at the cold start.

          ! Make use of sfc%ls_mask in order to identify potentially problematic points, 
          ! where depth_lk>0 (lake point in ICON) but ls_mask >0.5 (land point in IFS).
          ! At these points, tskin should not be used to initialize the water temperature.

          IF (llake) THEN
            CALL flake_coldinit(                                        &
              &     nflkgb      = ext_data(jg)%atm%fp_count    (jb), &  ! in
              &     idx_lst_fp  = ext_data(jg)%atm%idx_lst_fp(:,jb), &  ! in
              &     depth_lk    = ext_data(jg)%atm%depth_lk  (:,jb), &  ! in
              &     tskin       = initicon(jg)%sfc%tskin     (:,jb), &  ! in
              &     t_snow_lk_p = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_snow_lk(:,jb), &
              &     h_snow_lk_p = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_snow_lk(:,jb), &
              &     t_ice_p     = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_ice    (:,jb), &
              &     h_ice_p     = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ice    (:,jb), &
              &     t_mnw_lk_p  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_mnw_lk (:,jb), &
              &     t_wml_lk_p  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_wml_lk (:,jb), & 
              &     t_bot_lk_p  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_bot_lk (:,jb), &
              &     c_t_lk_p    = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%c_t_lk   (:,jb), &
              &     h_ml_lk_p   = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ml_lk  (:,jb), &
              &     t_b1_lk_p   = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_b1_lk  (:,jb), &
              &     h_b1_lk_p   = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_b1_lk  (:,jb), &
              &     t_g_lk_p    = p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_g_t    (:,jb,isub_lake) )
          ENDIF  ! llake

        ENDIF   ! inwp_surface > 0
      ENDDO
!$OMP END DO NOWAIT

    ENDDO
!$OMP END PARALLEL

  END SUBROUTINE copy_initicon2prog_sfc




  !-------------
  !>
  !! SUBROUTINE allocate_initicon
  !! Allocates the components of the initicon data type
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE allocate_initicon (p_patch, initicon)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)

    ! Local variables: loop control and dimensions
    INTEGER :: jg, nlev, nlevp1, nblks_c, nblks_e, mpi_comm

!-------------------------------------------------------------------------

    ! Loop over model domains
    DO jg = 1, n_dom

      nlev = p_patch(jg)%nlev
      nlevp1 = nlev + 1
      nblks_c = p_patch(jg)%nblks_c
      nblks_e = p_patch(jg)%nblks_e


      ! basic init_icon data
      ALLOCATE(initicon(jg)%topography_c    (nproma,nblks_c),        &
               initicon(jg)%z_ifc           (nproma,nlevp1,nblks_c), &
               initicon(jg)%z_mc            (nproma,nlev  ,nblks_c) )
      ! allocate groups for list of fields that must be read during initialization
      ALLOCATE(initicon(jg)%grp_vars_fg (200)        , &
               initicon(jg)%grp_vars_ana(200)        , &
               initicon(jg)%grp_vars_fg_default (200), &
               initicon(jg)%grp_vars_ana_default(200)  )

      ! Allocate atmospheric output data
      IF ( ANY((/MODE_IFSANA, MODE_DWDANA, MODE_COSMODE, MODE_COMBINED, MODE_ICONVREMAP/)==init_mode) ) THEN
        ALLOCATE(initicon(jg)%atm%vn        (nproma,nlev  ,nblks_e), &
                 initicon(jg)%atm%u         (nproma,nlev  ,nblks_c), &
                 initicon(jg)%atm%v         (nproma,nlev  ,nblks_c), &
                 initicon(jg)%atm%w         (nproma,nlevp1,nblks_c), &
                 initicon(jg)%atm%temp      (nproma,nlev  ,nblks_c), &
                 initicon(jg)%atm%exner     (nproma,nlev  ,nblks_c), &
                 initicon(jg)%atm%pres      (nproma,nlev  ,nblks_c), &
                 initicon(jg)%atm%rho       (nproma,nlev  ,nblks_c), &
                 initicon(jg)%atm%theta_v   (nproma,nlev  ,nblks_c), &
                 initicon(jg)%atm%qv        (nproma,nlev  ,nblks_c), &
                 initicon(jg)%atm%qc        (nproma,nlev  ,nblks_c), &
                 initicon(jg)%atm%qi        (nproma,nlev  ,nblks_c), &
                 initicon(jg)%atm%qr        (nproma,nlev  ,nblks_c), &
                 initicon(jg)%atm%qs        (nproma,nlev  ,nblks_c)  )

        initicon(jg)%atm%linitialized = .TRUE.
      ENDIF

      IF ( init_mode == MODE_ICONVREMAP ) THEN
        ALLOCATE(initicon(jg)%atm%tke(nproma,nlevp1,nblks_c))
      ENDIF

      ! Allocate surface output data
      ! always allocate sst (to be on the safe side)
      ALLOCATE(initicon(jg)%sfc%sst(nproma,nblks_c))

      ! initialize with 0 in order to avoid checks in the parent-to-child interpolation
      ! routine
      initicon(jg)%sfc%sst(:,:) = 0._wp

      IF ( init_mode == MODE_IFSANA ) THEN
        ALLOCATE(initicon(jg)%sfc%tskin    (nproma,nblks_c             ), &
                 initicon(jg)%sfc%tsnow    (nproma,nblks_c             ), &
                 initicon(jg)%sfc%snowalb  (nproma,nblks_c             ), &
                 initicon(jg)%sfc%snowweq  (nproma,nblks_c             ), &
                 initicon(jg)%sfc%snowdens (nproma,nblks_c             ), &
                 initicon(jg)%sfc%skinres  (nproma,nblks_c             ), &
                 initicon(jg)%sfc%ls_mask  (nproma,nblks_c             ), &
                 initicon(jg)%sfc%seaice   (nproma,nblks_c             ), &
                 initicon(jg)%sfc%tsoil    (nproma,0:nlev_soil,nblks_c ), &
                 initicon(jg)%sfc%wsoil    (nproma,  nlev_soil,nblks_c )  )
                 ! note the flipped dimensions with respect to sfc_in!
        initicon(jg)%sfc%linitialized = .TRUE.
      ENDIF


      ! atmospheric assimilation increments
      IF ( ANY((/MODE_DWDANA_INC, MODE_IAU, MODE_IAU_OLD/) == init_mode) ) THEN
        ALLOCATE(initicon(jg)%atm_inc%temp (nproma,nlev,nblks_c      ), &
                 initicon(jg)%atm_inc%pres (nproma,nlev,nblks_c      ), &
                 initicon(jg)%atm_inc%u    (nproma,nlev,nblks_c      ), &
                 initicon(jg)%atm_inc%v    (nproma,nlev,nblks_c      ), &
                 initicon(jg)%atm_inc%vn   (nproma,nlev,nblks_e      ), &
                 initicon(jg)%atm_inc%qv   (nproma,nlev,nblks_c      )  )
      
        initicon(jg)%atm_inc%linitialized = .TRUE.
      ENDIF

      ! surface assimilation increments
      IF ( (init_mode == MODE_IAU) .OR. (init_mode == MODE_IAU_OLD) ) THEN
        ALLOCATE(initicon(jg)%sfc_inc%w_so (nproma,nlev_soil,nblks_c ) )

        ! initialize with 0, since some increments are only read 
        ! for specific times
!$OMP PARALLEL WORKSHARE
        initicon(jg)%sfc_inc%w_so(:,:,:) = 0._wp
!$OMP END PARALLEL WORKSHARE

        initicon(jg)%sfc_inc%linitialized = .TRUE.

        ! allocate additional fields for MODE_IAU
        IF (init_mode == MODE_IAU) THEN
          ALLOCATE(initicon(jg)%sfc_inc%h_snow    (nproma,nblks_c ), &
            &      initicon(jg)%sfc_inc%freshsnow (nproma,nblks_c )  )

        ! initialize with 0, since some increments are only read 
        ! for specific times
!$OMP PARALLEL WORKSHARE
        initicon(jg)%sfc_inc%h_snow   (:,:) = 0._wp
        initicon(jg)%sfc_inc%freshsnow(:,:) = 0._wp
!$OMP END PARALLEL WORKSHARE
        ENDIF  ! MODE_IAU

      ENDIF



    ENDDO ! loop over model domains

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test 
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! read the map file into dictionary data structure:
    CALL dict_init(ana_varnames_dict, lcase_sensitive=.FALSE.)
    IF(ana_varnames_map_file /= ' ') THEN
      IF (my_process_is_stdio()) THEN 
        CALL dict_loadfile(ana_varnames_dict, TRIM(ana_varnames_map_file))
      END IF
      CALL p_bcast(ana_varnames_dict%nmax_entries,     p_io, mpi_comm)
      CALL p_bcast(ana_varnames_dict%nentries,         p_io, mpi_comm)
      CALL p_bcast(ana_varnames_dict%lcase_sensitive,  p_io, mpi_comm)
      IF (.NOT. my_process_is_stdio()) THEN 
        CALL dict_resize(ana_varnames_dict, ana_varnames_dict%nmax_entries)
      END IF
      CALL p_bcast(ana_varnames_dict%array(1,:), p_io, mpi_comm)
      CALL p_bcast(ana_varnames_dict%array(2,:), p_io, mpi_comm)
    END IF

  END SUBROUTINE allocate_initicon


  !-------------
  !>
  !! SUBROUTINE allocate_extana_atm
  !! Allocates fields for reading in external analysis data
  !!
  SUBROUTINE allocate_extana_atm (jg, nblks_c, nblks_e, initicon)
    INTEGER,                INTENT(IN)    :: jg, nblks_c, nblks_e
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)
    ! Local variables: loop control and dimensions
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = modname//':allocate_extana_atm'

    IF (nlev_in == 0) THEN
      CALL finish(routine, "Number of input levels <nlev_in> not yet initialized.")
    END IF

    ! Allocate atmospheric input data
    ALLOCATE( &
      initicon(jg)%atm_in%psfc    (nproma,        nblks_c),   &
      initicon(jg)%atm_in%phi_sfc (nproma,        nblks_c),   &
      initicon(jg)%atm_in%pres    (nproma,nlev_in,nblks_c),   &
      initicon(jg)%atm_in%z3d     (nproma,nlev_in,nblks_c),   &
      initicon(jg)%atm_in%temp    (nproma,nlev_in,nblks_c),   &
      initicon(jg)%atm_in%u       (nproma,nlev_in,nblks_c),   &
      initicon(jg)%atm_in%v       (nproma,nlev_in,nblks_c),   &
      initicon(jg)%atm_in%vn      (nproma,nlev_in,nblks_e),   &
      initicon(jg)%atm_in%w       (nproma,nlev_in,nblks_c),   &
      initicon(jg)%atm_in%omega   (nproma,nlev_in,nblks_c),   &
      initicon(jg)%atm_in%qv      (nproma,nlev_in,nblks_c),   &
      initicon(jg)%atm_in%qc      (nproma,nlev_in,nblks_c),   &
      initicon(jg)%atm_in%qi      (nproma,nlev_in,nblks_c),   &
      initicon(jg)%atm_in%qr      (nproma,nlev_in,nblks_c),   &
      initicon(jg)%atm_in%qs      (nproma,nlev_in,nblks_c)    )

    IF (init_mode == MODE_COSMODE .OR. init_mode == MODE_ICONVREMAP) THEN
      ALLOCATE( &
        initicon(jg)%atm_in%z3d_ifc (nproma,nlev_in+1,nblks_c), &
        initicon(jg)%atm_in%w_ifc   (nproma,nlev_in+1,nblks_c)  )
    ENDIF

    IF (init_mode == MODE_ICONVREMAP) THEN
      ALLOCATE( &
        initicon(jg)%atm_in%rho     (nproma,nlev_in  ,nblks_c), &
        initicon(jg)%atm_in%theta_v (nproma,nlev_in  ,nblks_c), &
        initicon(jg)%atm_in%tke     (nproma,nlev_in  ,nblks_c), &
        initicon(jg)%atm_in%tke_ifc (nproma,nlev_in+1,nblks_c)  )
    ENDIF

    initicon(jg)%atm_in%linitialized = .TRUE.
  END SUBROUTINE allocate_extana_atm


  !-------------
  !>
  !! SUBROUTINE allocate_extana_sfc
  !! Allocates fields for reading in external analysis data
  !!
  SUBROUTINE allocate_extana_sfc (jg, nblks_c, initicon)
    INTEGER,                INTENT(IN)    :: jg, nblks_c
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)

    ! Allocate surface input data
    ! The extra soil temperature levels are not read in; they are only used to simplify vertical interpolation
    ALLOCATE(initicon(jg)%sfc_in%phi      (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%tskin    (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%sst      (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%tsnow    (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%snowalb  (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%snowweq  (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%snowdens (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%skinres  (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%ls_mask  (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%seaice   (nproma,nblks_c                ), &
      initicon(jg)%sfc_in%tsoil    (nproma,nblks_c,0:nlevsoil_in+1), &
      initicon(jg)%sfc_in%wsoil    (nproma,nblks_c,0:nlevsoil_in+1)  )
    initicon(jg)%sfc_in%linitialized = .TRUE.
  END SUBROUTINE allocate_extana_sfc


  !-------------
  !>
  !! SUBROUTINE deallocate_initicon
  !! Deallocates the components of the initicon data type
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE deallocate_initicon (initicon)

    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)

    ! Local variables: loop control
    INTEGER :: jg

!-------------------------------------------------------------------------

    ! Loop over model domains
    DO jg = 1, n_dom

      ! basic init_icon data
      DEALLOCATE(initicon(jg)%topography_c,     &
                 initicon(jg)%z_ifc,            &
                 initicon(jg)%z_mc              )
      ! deallocate groups for list of fields that must be read during initialization
      DEALLOCATE(initicon(jg)%grp_vars_fg,         &
                 initicon(jg)%grp_vars_fg_default, &
                 initicon(jg)%grp_vars_ana,        &
                 initicon(jg)%grp_vars_ana_default )

      ! atmospheric output data
      IF (initicon(jg)%atm%linitialized) THEN
        DEALLOCATE(initicon(jg)%atm%vn,      &
                   initicon(jg)%atm%u,       &
                   initicon(jg)%atm%v,       &
                   initicon(jg)%atm%w,       &
                   initicon(jg)%atm%temp,    &
                   initicon(jg)%atm%exner,   &
                   initicon(jg)%atm%pres,    &  
                   initicon(jg)%atm%rho,     &
                   initicon(jg)%atm%theta_v, &
                   initicon(jg)%atm%qv,      &
                   initicon(jg)%atm%qc,      &
                   initicon(jg)%atm%qi,      &
                   initicon(jg)%atm%qr,      &
                   initicon(jg)%atm%qs       )
      ENDIF

      IF ( init_mode == MODE_ICONVREMAP ) THEN
        DEALLOCATE(initicon(jg)%atm%tke)
      ENDIF

      ! always allocated (hack!)
      DEALLOCATE(initicon(jg)%sfc%sst)

      ! surface output data
      IF (initicon(jg)%sfc%linitialized) THEN
        DEALLOCATE(initicon(jg)%sfc%tskin,    &
                   initicon(jg)%sfc%tsnow,    &
                   initicon(jg)%sfc%snowalb,  &
                   initicon(jg)%sfc%snowweq,  &
                   initicon(jg)%sfc%snowdens, &
                   initicon(jg)%sfc%skinres,  &
                   initicon(jg)%sfc%ls_mask,  &
                   initicon(jg)%sfc%seaice,   &
                   initicon(jg)%sfc%tsoil,    &
                   initicon(jg)%sfc%wsoil     )
      ENDIF


      ! atmospheric assimilation increments
      IF (initicon(jg)%atm_inc%linitialized) THEN
        DEALLOCATE(initicon(jg)%atm_inc%temp,    &
                   initicon(jg)%atm_inc%pres,    &
                   initicon(jg)%atm_inc%u   ,    &
                   initicon(jg)%atm_inc%v   ,    &
                   initicon(jg)%atm_inc%vn  ,    &
                   initicon(jg)%atm_inc%qv       )
      ENDIF


      ! surface assimilation increments
      IF ( initicon(jg)%sfc_inc%linitialized ) THEN
        DEALLOCATE(initicon(jg)%sfc_inc%w_so )
        IF (ALLOCATED(initicon(jg)%sfc_inc%h_snow))    DEALLOCATE(initicon(jg)%sfc_inc%h_snow )
        IF (ALLOCATED(initicon(jg)%sfc_inc%freshsnow)) DEALLOCATE(initicon(jg)%sfc_inc%freshsnow )
      ENDIF


    ENDDO ! loop over model domains

    ! destroy variable name dictionaries:
    CALL dict_finalize(ana_varnames_dict)

  END SUBROUTINE deallocate_initicon


  !-------------
  !>
  !! SUBROUTINE deallocate_extana_atm
  !! Deallocates the components of the initicon data type
  !!
  SUBROUTINE deallocate_extana_atm (initicon)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)
    ! Local variables: loop control
    INTEGER :: jg

    ! Loop over model domains
    DO jg = 1, n_dom
      IF (.NOT. initicon(jg)%atm_in%linitialized) CYCLE

      ! atmospheric input data
      DEALLOCATE(initicon(jg)%atm_in%psfc,    &
                 initicon(jg)%atm_in%phi_sfc, &
                 initicon(jg)%atm_in%pres,    &
                 initicon(jg)%atm_in%temp,    &
                 initicon(jg)%atm_in%u,       &
                 initicon(jg)%atm_in%v,       &
                 initicon(jg)%atm_in%vn,      &
                 initicon(jg)%atm_in%w,       &
                 initicon(jg)%atm_in%z3d,     &
                 initicon(jg)%atm_in%omega,   &
                 initicon(jg)%atm_in%qv,      &
                 initicon(jg)%atm_in%qc,      &
                 initicon(jg)%atm_in%qi,      &
                 initicon(jg)%atm_in%qr,      &
                 initicon(jg)%atm_in%qs )

      IF (init_mode == MODE_COSMODE .OR. init_mode == MODE_ICONVREMAP) THEN
        DEALLOCATE( &
                 initicon(jg)%atm_in%z3d_ifc, &
                 initicon(jg)%atm_in%w_ifc    )
      ENDIF

      IF (init_mode == MODE_ICONVREMAP) THEN
        DEALLOCATE( &
          initicon(jg)%atm_in%rho,     &
          initicon(jg)%atm_in%theta_v, &
          initicon(jg)%atm_in%tke,     &
          initicon(jg)%atm_in%tke_ifc  )
      ENDIF

      initicon(jg)%atm_in%linitialized = .FALSE.
    ENDDO ! loop over model domains

  END SUBROUTINE deallocate_extana_atm


  !-------------
  !>
  !! SUBROUTINE deallocate_extana_sfc
  !! Deallocates the components of the initicon data type
  !!
  SUBROUTINE deallocate_extana_sfc (initicon)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)
    ! Local variables: loop control
    INTEGER :: jg

    ! Loop over model domains
    DO jg = 1, n_dom
      IF (.NOT. initicon(jg)%sfc_in%linitialized) CYCLE
      ! surface input data
      DEALLOCATE(initicon(jg)%sfc_in%phi,      &
                 initicon(jg)%sfc_in%tskin,    &
                 initicon(jg)%sfc_in%sst,      &
                 initicon(jg)%sfc_in%tsnow,    &
                 initicon(jg)%sfc_in%snowalb,  &
                 initicon(jg)%sfc_in%snowweq,  &
                 initicon(jg)%sfc_in%snowdens, &
                 initicon(jg)%sfc_in%skinres,  &
                 initicon(jg)%sfc_in%ls_mask,  &
                 initicon(jg)%sfc_in%seaice,   &
                 initicon(jg)%sfc_in%tsoil,    &
                 initicon(jg)%sfc_in%wsoil     )
      initicon(jg)%sfc_in%linitialized = .FALSE.
    ENDDO ! loop over model domains

  END SUBROUTINE deallocate_extana_sfc

END MODULE mo_initicon_utils

