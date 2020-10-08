!>
!! Utility funtions for handling field specific meta information
!!
!! Contains utility funtions which are used for defining variable specific
!! meta information. These have nothing to do with var lists itself. That's
!! why they have been moved here from mo_var_list.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2014-01-22)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_var_metadata

  USE mo_kind,               ONLY: wp, sp, dp
  USE mo_exception,          ONLY: finish
  USE mo_cf_convention,      ONLY: t_cf_var
  USE mo_grib2,              ONLY: t_grib2_var
  USE mo_var_metadata_types, ONLY: t_hor_interp_meta, t_vert_interp_meta, &
    &  t_post_op_meta, VINTP_TYPE_LIST, t_var_metadata, t_union_vals, &
    &  t_var_metadata_dynamic
  USE mo_util_string,        ONLY: toupper
  USE mo_action_types,       ONLY: t_var_action
  USE mo_util_texthash,      ONLY: text_hash_c
  USE mo_tracer_metadata_types, ONLY: t_tracer_meta
  USE mo_tracer_metadata,    ONLY: create_tracer_metadata
  USE mo_util_string,        ONLY: tolower
  USE mo_impl_constants,     ONLY: TIMELEVEL_SUFFIX, MAX_TIME_LEVELS, vname_len

  IMPLICIT NONE
  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_var_metadata'

  PUBLIC :: create_hor_interp_metadata, create_vert_interp_metadata
  PUBLIC :: post_op, vintp_types, vintp_type_id, set_var_metadata
  PUBLIC :: set_var_metadata_dyn
  PUBLIC :: get_var_name ! return plain variable name (without timelevel)
  PUBLIC :: get_var_timelevel ! return variable timelevel (or "-1")
  PUBLIC :: get_timelevel_string ! return the default string with timelevel encoded

CONTAINS
  !------------------------------------------------------------------------------------------------
  !
  ! Quasi-constructor for horizontal interpolation meta data
  !
  ! Fills data structure with default values (unless set otherwise).
  FUNCTION create_hor_interp_metadata(hor_intp_type, fallback_type, lonlat_id)    &
    & RESULT(him)
    TYPE(t_hor_interp_meta) :: him
    INTEGER, INTENT(IN), OPTIONAL :: hor_intp_type, fallback_type, lonlat_id

    ! supersede with user definitions
    IF (PRESENT(hor_intp_type)) him%hor_intp_type = hor_intp_type
    IF (PRESENT(fallback_type)) him%fallback_type = fallback_type
    IF (PRESENT(lonlat_id))     him%lonlat_id     = lonlat_id
  END FUNCTION create_hor_interp_metadata

  !------------------------------------------------------------------------------------------------
  ! HANDLING OF VERTICAL INTERPOLATION MODES
  !------------------------------------------------------------------------------------------------

  !> Implements a (somewhat randomly chosen) one-to-one mapping
  !  between a string and an integer ID number between 1 and
  !  MAX_VINTP_TYPES.
  !
  FUNCTION vintp_type_id(in_str)
    INTEGER                      :: vintp_type_id, ivintp_type
    CHARACTER(LEN=*), INTENT(IN) :: in_str
    CHARACTER(len=*), PARAMETER :: routine = modname//"::vintp_type_id"
    CHARACTER(len=len_trim(in_str)) :: in_str_upper
    INTEGER :: n

    vintp_type_id = 0
    in_str_upper = toupper(in_str)
    n = SIZE(VINTP_TYPE_LIST)
    LOOP_VINTP_TYPES : DO ivintp_type=1,n
      IF (in_str_upper == VINTP_TYPE_LIST(ivintp_type)) THEN
        vintp_type_id = ivintp_type
        EXIT LOOP_VINTP_TYPES
      END IF
    END DO LOOP_VINTP_TYPES
    ! paranoia:
    IF ((vintp_type_id < 1) .OR. (vintp_type_id > n)) &
      &  CALL finish(routine, "Invalid vertical interpolation type!")
  END FUNCTION vintp_type_id

  !> Utility function with *a lot* of optional string parameters v1,
  !  v2, v3, v4, ...; mapping those onto a
  !  LOGICAL(:) according to the "vintp_type_id"
  !  function.
  !
  FUNCTION vintp_types(v01, v02, v03, v04, v05, v06, v07, v08, v09, v10)
    LOGICAL :: vintp_types(SIZE(VINTP_TYPE_LIST))
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: &
      &   v01, v02, v03, v04, v05, v06, v07, v08, v09, v10

    vintp_types(:) = .FALSE.
    IF (PRESENT(v01)) vintp_types(vintp_type_id(v01)) = .TRUE.
    IF (PRESENT(v02)) vintp_types(vintp_type_id(v02)) = .TRUE.
    IF (PRESENT(v03)) vintp_types(vintp_type_id(v03)) = .TRUE.
    IF (PRESENT(v04)) vintp_types(vintp_type_id(v04)) = .TRUE.
    IF (PRESENT(v05)) vintp_types(vintp_type_id(v05)) = .TRUE.
    IF (PRESENT(v06)) vintp_types(vintp_type_id(v06)) = .TRUE.
    IF (PRESENT(v07)) vintp_types(vintp_type_id(v07)) = .TRUE.
    IF (PRESENT(v08)) vintp_types(vintp_type_id(v08)) = .TRUE.
    IF (PRESENT(v09)) vintp_types(vintp_type_id(v09)) = .TRUE.
    IF (PRESENT(v10)) vintp_types(vintp_type_id(v10)) = .TRUE.
  END FUNCTION vintp_types

  !------------------------------------------------------------------------------------------------
  !
  ! Quasi-constructor for vertical interpolation meta data
  !
  ! Fills data structure with default values (unless set otherwise).
  FUNCTION create_vert_interp_metadata(vert_intp_type, vert_intp_method, l_hires_intp,       &
    & l_restore_fricred, l_loglin, l_extrapol, l_satlimit, l_restore_pbldev, l_pd_limit,     &
    & lower_limit) RESULT(vim)
    TYPE(t_vert_interp_meta) :: vim
    LOGICAL, INTENT(IN), OPTIONAL :: vert_intp_type(SIZE(VINTP_TYPE_LIST))
    INTEGER, INTENT(IN), OPTIONAL :: vert_intp_method 
    LOGICAL, INTENT(IN), OPTIONAL :: l_hires_intp, l_restore_fricred, &
      & l_loglin, l_extrapol, l_satlimit, l_restore_pbldev, l_pd_limit
    REAL(wp), INTENT(IN), OPTIONAL :: lower_limit

    ! supersede with user definitions
    IF (PRESENT(vert_intp_type))    vim%vert_intp_type    = vert_intp_type
    IF (PRESENT(vert_intp_method))  vim%vert_intp_method  = vert_intp_method
    IF (PRESENT(l_hires_intp))      vim%l_hires_intp      = l_hires_intp
    IF (PRESENT(l_restore_fricred)) vim%l_restore_fricred = l_restore_fricred
    IF (PRESENT(l_loglin))          vim%l_loglin          = l_loglin
    IF (PRESENT(l_extrapol))        vim%l_extrapol        = l_extrapol
    IF (PRESENT(l_satlimit))        vim%l_satlimit        = l_satlimit
    IF (PRESENT(l_restore_pbldev))  vim%l_restore_pbldev  = l_restore_pbldev
    IF (PRESENT(l_pd_limit))        vim%l_pd_limit        = l_pd_limit
    IF (PRESENT(lower_limit))       vim%lower_limit       = lower_limit
  END FUNCTION create_vert_interp_metadata

  !----------------------------------------------------------------------------------------
  !
  FUNCTION post_op(ipost_op_type, new_cf, new_grib2, arg1)
    TYPE(t_post_op_meta) :: post_op
    INTEGER,           INTENT(IN), OPTIONAL :: ipost_op_type    !< type of post-processing operation
    TYPE(t_cf_var),    INTENT(IN), OPTIONAL :: new_cf           !< CF information of modified field
    TYPE(t_grib2_var), INTENT(IN), OPTIONAL :: new_grib2        !< GRIB2 information of modified field
    CLASS(*),          INTENT(IN), OPTIONAL :: arg1             !< post-op argument (e.g. scaling factor)

    IF (PRESENT(ipost_op_type)) post_op%ipost_op_type = ipost_op_type
    IF (PRESENT(arg1)) THEN
      SELECT TYPE(arg1)
      TYPE is (INTEGER)
        post_op%arg1%ival = arg1
      TYPE is (REAL(dp))
        post_op%arg1%rval = arg1
      TYPE is (REAL(sp))
        post_op%arg1%sval = arg1
      END SELECT
    ENDIF
    IF (PRESENT(new_cf)) THEN
      post_op%lnew_cf = .TRUE.
      post_op%new_cf  = new_cf
    END IF
    IF (PRESENT(new_grib2)) THEN
      post_op%lnew_grib2 = .TRUE.
      post_op%new_grib2  = new_grib2
    END IF
  END FUNCTION post_op

  !------------------------------------------------------------------------------------------------
  ! Set parameters of list element already created
  ! (private routine within this module)
  !
  ! Set each parameter in data type var_metadata if the respective
  ! optional parameter is present.
  SUBROUTINE set_var_metadata (info,                                           &
         &                     name, hgrid, vgrid, cf, grib2, ldims,           &
         &                     loutput, lcontainer, lrestart, lrestart_cont,   &
         &                     initval, isteptype, resetval, lmiss, missval,   &
         &                     tlev_source, vert_interp,                       &
         &                     hor_interp, in_group, verbose,                  &
         &                     l_pp_scheduler_task, post_op, action_list,      &
         &                     var_class, data_type, idx_tracer, idx_diag,     &
         &                     lopenacc)
    TYPE(t_var_metadata),    INTENT(inout)        :: info          ! memory info struct.
    CHARACTER(len=*),        INTENT(in), OPTIONAL :: name          ! variable name
    INTEGER,                 INTENT(in), OPTIONAL :: hgrid         ! horizontal grid type used
    INTEGER,                 INTENT(in), OPTIONAL :: vgrid         ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in), OPTIONAL :: cf            ! CF convention
    TYPE(t_grib2_var),       INTENT(in), OPTIONAL :: grib2         ! GRIB2
    INTEGER,                 INTENT(in)           :: ldims(:)      ! used dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput       ! into output var_list
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer    ! true if container
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart      ! restart file flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont ! continue on restart
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: initval       ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype     ! type of statistical processing
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: resetval      ! reset value
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss         ! missing value flag
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: missval       ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source   ! actual TL for TL dependent vars
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp   ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp    ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(:)          ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(in), OPTIONAL :: post_op       !< "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(in), OPTIONAL :: action_list   !< regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class     ! variable class/species
    INTEGER,                 INTENT(IN), OPTIONAL :: data_type     ! variable data type
    INTEGER,                 INTENT(IN), OPTIONAL :: idx_tracer    ! index of tracer in tracer container 
    INTEGER,                 INTENT(IN), OPTIONAL :: idx_diag      ! index of tracer in diagnostics container 
    LOGICAL,                 INTENT(IN), OPTIONAL :: lopenacc      ! variable data type
    !LOGICAL :: lverbose

    ! set flags from optional parameters
    !lverbose = .FALSE.
    IF (PRESENT(name)) THEN
      info%name      = name
      info%key = text_hash_c(TRIM(name))
      info%key_notl = text_hash_c(tolower(strip_timelev(name)))
    END IF
    !IF (PRESENT(verbose))       lverbose             = verbose
    IF (PRESENT(data_type))     info%data_type       = data_type
    ! set components describing the 'Content of the field'
    IF (PRESENT(var_class))     info%var_class       = var_class
    IF (PRESENT(cf))            info%cf              = cf
    IF (PRESENT(grib2))         info%grib2           = grib2
    IF (PRESENT(hgrid))         info%hgrid           = hgrid
    IF (PRESENT(vgrid))         info%vgrid           = vgrid
    info%used_dimensions(:SIZE(ldims)) = ldims(:)
    IF (SIZE(ldims) .LT. 5) info%used_dimensions(SIZE(ldims)+1:) = 1
    IF (PRESENT(loutput))       info%loutput         = loutput
    IF (PRESENT(lcontainer))    info%lcontainer      = lcontainer
    IF (info%lcontainer) THEN
      info%ncontained   =  0
      info%var_ref_pos  = -1 ! UNDEFINED
    END IF
    IF (PRESENT(resetval))      info%resetval      = resetval
    IF (PRESENT(isteptype))     info%isteptype     = isteptype
    IF (PRESENT(lmiss))         info%lmiss         = lmiss
    IF (PRESENT(missval))       info%missval       = missval
    IF (PRESENT(lrestart))      info%lrestart      = lrestart
    IF (PRESENT(lrestart_cont)) info%lrestart_cont = lrestart_cont
    IF (PRESENT(initval))       info%initval       = initval
    IF (PRESENT(tlev_source))   info%tlev_source   = tlev_source
    ! set flags concerning vertical interpolation
    IF (PRESENT(vert_interp))   info%vert_interp   = vert_interp
    IF (PRESENT(hor_interp))    info%hor_interp    = hor_interp
    IF (PRESENT(in_group)) &
      & info%in_group(:SIZE(in_group)) = in_group(:)
    IF (PRESENT(l_pp_scheduler_task)) &
      & info%l_pp_scheduler_task = l_pp_scheduler_task
    IF (PRESENT(post_op))       info%post_op       = post_op
    IF (PRESENT(action_list))   info%action_list   = action_list
    ! indices of tracer in tracer container and in diagnostic container
    IF (PRESENT(idx_tracer))    info%idx_tracer    = idx_tracer
    IF (PRESENT(idx_diag))      info%idx_diag      = idx_diag
    IF (PRESENT(lopenacc))      info%lopenacc      = lopenacc
    ! perform consistency checks on variable's meta-data:
    CALL check_metadata_consistency()
    !LK    IF (lverbose) CALL print_var_metadata (info)
  CONTAINS

    SUBROUTINE check_metadata_consistency()
      CHARACTER(*), PARAMETER :: routine = modname//':check_metadata_consistency'

      IF (info%lrestart .AND. info%lcontainer) &
        & CALL finish(routine//' - '//TRIM(info%name), &
          & 'Container variables are not restartable! Use var references instead.')
    ! ... put other consistency checks here ...
    END SUBROUTINE check_metadata_consistency

  END SUBROUTINE set_var_metadata
  !------------------------------------------------------------------------------------------------
  !
  ! Set dynamic metadata, i.e. polymorphic tracer metadata
  !
  SUBROUTINE set_var_metadata_dyn(this_info_dyn,tracer_info)
    TYPE(t_var_metadata_dynamic),INTENT(INOUT) :: this_info_dyn
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info

    IF (PRESENT(tracer_info)) THEN
      ALLOCATE(this_info_dyn%tracer, source=tracer_info)
    ELSE
      ALLOCATE(t_tracer_meta :: this_info_dyn%tracer)
      SELECT TYPE(tm => this_info_dyn%tracer)
      TYPE IS(t_tracer_meta)
        tm = create_tracer_metadata(lis_tracer=.FALSE.)
      END SELECT
    ENDIF
  END SUBROUTINE set_var_metadata_dyn

  !------------------------------------------------------------------------------------------------
  !> @return Plain variable name (i.e. without TIMELEVEL_SUFFIX)
  CHARACTER(LEN=vname_len) FUNCTION get_var_name(info)
    TYPE(t_var_metadata), INTENT(IN) :: info

    get_var_name = strip_timelev(info%name)
  END FUNCTION get_var_name

  CHARACTER(LEN=vname_len) FUNCTION strip_timelev(vname)
    CHARACTER(*), INTENT(IN)   :: vname
    INTEGER :: idx

    idx = INDEX(vname,TIMELEVEL_SUFFIX)
    IF (idx .EQ. 0) THEN
      strip_timelev = vname
    ELSE
      strip_timelev = vname(1:idx-1)
    END IF
  END FUNCTION strip_timelev

  !------------------------------------------------------------------------------------------------
  ! construct string for timelevel encoding into variable names
  CHARACTER(len=4) FUNCTION get_timelevel_string(timelevel)
    INTEGER, INTENT(IN) :: timelevel

    WRITE(get_timelevel_string,'("'//TIMELEVEL_SUFFIX//'",i1)') timelevel
  END FUNCTION get_timelevel_string

  !------------------------------------------------------------------------------------------------
  !> @return time level (extracted from time level suffix) or "-1"
  INTEGER FUNCTION get_var_timelevel(vname) RESULT(tl)
    CHARACTER(*), INTENT(IN) :: vname
    CHARACTER(*), PARAMETER :: routine = modname//':get_var_timelevel'

    tl = INDEX(vname,TIMELEVEL_SUFFIX)
    IF (tl .EQ. 0) THEN
      tl = -1
    ELSE
      tl = ICHAR(vname(tl+3:tl+3)) - ICHAR('0')
      IF (tl .LE. 0 .OR. tl .GT. MAX_TIME_LEVELS) &
        & CALL finish(routine, 'Illegal time level in '//TRIM(vname))
    END IF
  END FUNCTION get_var_timelevel
END MODULE mo_var_metadata

