!! Group information for ICON variables.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_var_groups

  USE mo_impl_constants,        ONLY: VARNAME_LEN, TIMELEVEL_SUFFIX
  USE mo_exception,             ONLY: finish
  USE mo_util_string,           ONLY: toupper, int2string
  USE mo_fortran_tools,         ONLY: resize_arr_c1d
  USE mo_util_sort,             ONLY: quicksort


  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_var_groups'


  ! maximum number of variable groups supported by a single info state
  INTEGER, PARAMETER :: MAX_GROUPS = 120





  ! ---------------------------------------------------------------
  ! STATICALLY DEFINED VARIABLE GROUPS
  ! ---------------------------------------------------------------

  ! List of *static* variable groups.
  ! 
  ! A variable can have any combination of this which means that it is
  ! part of each of these different variable sets.
  ! A variable is added to an existing group by setting the meta-data
  ! information "in_group" as follows
  !
  !   CALL add_var( p_prog_list, ..., in_group=groups("nh_prog_vars") )
  !
  ! It is also possible to add a variable to more than one group:
  !
  !   CALL add_var( diag_list, ...,   &
  !                 in_group=groups("multisnow_vars", "snow_vars"))
  !
  ! New groups can be added by extending the VAR_GROUPS list.
  !
  ! Note that the statically defined group list "var_groups" is
  ! non-public. Its contents are copied to a dynamically growing list
  ! "var_groups_dyn".

  CHARACTER(len=VARNAME_LEN), PARAMETER :: VAR_GROUPS_STATIC(61) = &
    (/ "ALL                   ",  &
    &  "ATMO_ML_VARS          ",  &
    &  "ATMO_PL_VARS          ",  &
    &  "ATMO_ZL_VARS          ",  &
    &  "NH_PROG_VARS          ",  &
    &  "ATMO_DERIVED_VARS     ",  &
    &  "RAD_VARS              ",  &
    &  "PRECIP_VARS           ",  &
    &  "CLOUD_DIAG            ",  &
    &  "PBL_VARS              ",  &
    &  "PHYS_TENDENCIES       ",  &
    &  "PROG_TIMEMEAN         ",  &
    &  "ECHAM_TIMEMEAN        ",  &
    &  "TRACER_TIMEMEAN       ",  &
    &  "ATMO_TIMEMEAN         ",  &
    &  "LAND_VARS             ",  &
    &  "LAND_TILE_VARS        ",  &
    &  "MULTISNOW_VARS        ",  &
    &  "ADDITIONAL_PRECIP_VARS",  &
    &  "SNOW_VARS             ",  &
    &  "DWD_FG_ATM_VARS       ",  &  ! DWD First Guess (atmosphere) 
    &  "DWD_FG_SFC_VARS       ",  &  ! DWD First Guess (surface/soil)
    &  "DWD_FG_SFC_VARS_T     ",  &  ! DWD First Guess (surface/soil) tiles
    &  "MODE_DWD_FG_IN        ",  &  ! Input first guess fields for MODE_DWD
    &  "MODE_DWD_ANA_IN       ",  &  ! Input analysis fields for MODE_DWD
    &  "MODE_IAU_FG_IN        ",  &  ! First guess input for IAU
    &  "MODE_IAU_ANA_IN       ",  &  ! Analysis input for IAU
    &  "MODE_IAU_ANAATM_IN    ",  &  ! Atmospheric analysis input for (old/new) IAU
    &  "MODE_IAU_OLD_FG_IN    ",  &  ! First guess input for old IAU mode
    &  "MODE_IAU_OLD_ANA_IN   ",  &  ! Analysis input for old IAU mode
    &  "MODE_COMBINED_IN      ",  &  ! Input fields for MODE_COMBINED
    &  "MODE_COSMO_IN         ",  &  ! Input fields for MODE_COSMO
    &  "OCE_PROG              ",  &
    &  "OCE_DIAG              ",  &
    &  "OCE_EDDY              ",  &
    &  "OCE_DEFAULT           ",  &
    &  "OCEAN_MOC             ",  &  ! meant o hold all kinds of overturning fields (atl, pac, global,...)
    &  "OCEAN_FLOWS           ",  &  ! meant o hold all through flows
    &  "HAMOCC_BASE           ",  &
    &  "HAMOCC_TEND           ",  &
    &  "HAMOCC_MONI           ",  &
    &  "HAMOCC_SED            ",  &
    &  "oce_essentials        ",  &
    &  "oce_force_essentials  ",  &
    &  "OCE_AUX               ",  &
    &  "OCEAN_MONITOR         ",  &
    &  "OCE_GEOMETRY          ",  &
    &  "OCE_PHYSICS           ",  &
    &  "OCE_COEFFS            ",  &
    &  "ICE_DEFAULT           ",  &
    &  "ICE_BUDGETS           ",  &
    &  "ICE_DIAG              ",  &
    &  "LATBC_PREFETCH_VARS   ",  &
    &  "mode_iniana           ",  &  ! Variable set needed to initialize ICON (MODE_ICONVREMAP)
    &  "icon_lbc_vars         ",  &  ! Variable set needed for ICON-LAM lateral boundary conditions
    &  "ART_AEROSOL           ",  &  ! ICON-ART fields for aerosol particles
    &  "ART_CHEMISTRY         ",  &  ! ICON-ART fields for chemical tracers
    &  "ART_PASSIVE           ",  &  ! ICON-ART fields for passive tracers
    &  "ART_DIAGNOSTICS       ",  &  ! ICON-ART fields for diagnostic fields
    &  "ART_ROUTINE_DIAG      ",  &  ! ICON-ART fields for routine diagnostic fields
    &  "RTTOV                 " /)


  ! ---------------------------------------------------------------
  ! DYNAMICALLY DEFINED VARIABLE GROUPS
  ! ---------------------------------------------------------------

  TYPE t_var_groups
    ! variable group names. The list ordering is important, since the
    ! entry #i (= ID #i) corresponds to the i'th LOGICAL entry in the
    ! "in_group" of the t_var_metadata type below.
    !
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: name(:)

  CONTAINS

    PROCEDURE :: init              => t_var_groups_init
    PROCEDURE :: finalize          => t_var_groups_finalize
    PROCEDURE :: group_id          => t_var_groups_group_id
    PROCEDURE :: add               => t_var_groups_add
    PROCEDURE :: alphabetical_list => t_var_groups_alphabetical_list
  END TYPE t_var_groups


  ! List of variable groups - including the statically defined group
  ! names from VAR_GROUPS as well as dynamic variable groups, used,
  ! e.g., for tiles.
  !
  TYPE(t_var_groups) :: var_groups_dyn


  INTERFACE groups
    MODULE PROCEDURE groups_arg
    MODULE PROCEDURE groups_vec
  END INTERFACE


  PUBLIC :: MAX_GROUPS
  PUBLIC :: t_var_groups
  PUBLIC :: var_groups_dyn
  PUBLIC :: groups


CONTAINS

  !----------------------------------------------------------------------------------------
  !> Copy all statically defined variable group names into the
  !  "var_groups_dyn" data structure. The statically defined group
  !  list is non-public.
  SUBROUTINE t_var_groups_init(var_groups)
    CLASS(t_var_groups), INTENT(INOUT) :: var_groups

    CHARACTER(*), PARAMETER :: routine = modname//"::t_var_groups_init"
    INTEGER :: istat ! status

    IF (ALLOCATED(var_groups%name)) THEN
      CALL finish(routine, "Internal error!")
    END IF
    ALLOCATE(var_groups%name(SIZE(VAR_GROUPS_STATIC)), STAT=istat)
    IF (istat /= 0)  CALL finish(routine, 'ALLOCATE failed!')
    ! add static groups
    var_groups%name(1:SIZE(VAR_GROUPS_STATIC)) = VAR_GROUPS_STATIC(1:SIZE(VAR_GROUPS_STATIC))
  END SUBROUTINE t_var_groups_init


  !----------------------------------------------------------------------------------------
  !> Destructor.
  !
  SUBROUTINE t_var_groups_finalize(var_groups)
    CLASS(t_var_groups), INTENT(INOUT) :: var_groups

    CHARACTER(*), PARAMETER :: routine = modname//"::t_var_groups_finalize"
    INTEGER :: istat ! status

    IF (ALLOCATED(var_groups%name)) THEN
      DEALLOCATE(var_groups%name, STAT=istat)
      IF (istat /= 0)  CALL finish(routine, 'DEALLOCATE failed!')
    END IF
  END SUBROUTINE t_var_groups_finalize


  !----------------------------------------------------------------------------------------
  !> Implements a (somewhat randomly chosen) one-to-one mapping
  !  between a string and an integer ID number between 1 and
  !  size(var_groups%name).
  !
  FUNCTION t_var_groups_group_id(var_groups, in_str, opt_lcheck)  RESULT(group_id)
    INTEGER                            :: group_id
    CLASS(t_var_groups), INTENT(INOUT) :: var_groups
    CHARACTER(LEN=*) ,   INTENT(IN)    :: in_str
    LOGICAL, OPTIONAL,   INTENT(IN)    :: opt_lcheck           
    !
    ! Local
    CHARACTER(*), PARAMETER :: routine = modname//"::t_var_groups_group_id"
    LOGICAL :: lcheck
    INTEGER :: max_size, igrp

    IF (.NOT. ALLOCATED(var_groups%name))  CALL var_groups%init()

    ! search the variable groups (which includes the statically
    ! defined groups and the dynamically defined ones):
    group_id = 0
    LOOP_GROUPS : DO igrp=1,SIZE(var_groups%name)
      IF (toupper(TRIM(in_str)) == toupper(TRIM(var_groups%name(igrp)))) THEN
        group_id = igrp
        EXIT LOOP_GROUPS
      END IF
    END DO LOOP_GROUPS

    ! If the group does not exist, create it.
    IF (group_id == 0) THEN
      !
      ! increase dynamic groups array by one element
      CALL resize_arr_c1d(var_groups%name,1)
      !
      ! add new group
      var_groups%name(SIZE(var_groups%name)) = toupper(TRIM(in_str))
      !
      ! return its group ID (including offset from static groups array)
      group_id = SIZE(var_groups%name)
    ENDIF

    ! paranoia:
    lcheck = .TRUE.
    IF (PRESENT(opt_lcheck))  lcheck = opt_lcheck
    IF (lcheck) THEN
      max_size = SIZE(var_groups%name)
      IF ((group_id < 1) .OR. (group_id > max_size)) &
        &  CALL finish(routine, "Invalid group ID: "//TRIM(in_str)//" = "//TRIM(int2string(group_id,"(i0)")))
    ENDIF
  END FUNCTION t_var_groups_group_id


  !----------------------------------------------------------------------------------------
  !> Add new (tile) member to variable group
  !
  !  Adds new tile member to variable-specific tile-group. 
  !  If the group does not exist, a group (named after the 
  !  corresponding container) is added to the dynamic variable 
  !  groups list first.
  ! 
  !  @par Revision History
  !  Initial revision by Daniel Reinert, DWD (2015-01-29)
  ! 
  SUBROUTINE t_var_groups_add(var_groups, group_name, in_group_new, opt_in_group)
    CLASS(t_var_groups), INTENT(INOUT) :: var_groups
    CHARACTER(len=*) ,   INTENT(in)    :: group_name
    LOGICAL          ,   INTENT(out)   :: in_group_new(:)
    LOGICAL, OPTIONAL,   INTENT(in)    :: opt_in_group(:)
    !
    ! Local
    CHARACTER(*), PARAMETER :: routine = modname//"::t_var_groups_group_add"
    INTEGER  :: idx, grp_id
    CHARACTER(len=LEN(group_name)) ::  group_name_plain

    ! check whether a group with name 'group_name_plain' exists and return its ID.
    !
    ! remove time level string from group name
    idx = INDEX(group_name, TIMELEVEL_SUFFIX)
    IF (idx > 0) THEN
      group_name_plain = TRIM(group_name(1:idx-1))
    ELSE
      group_name_plain = TRIM(group_name)
    ENDIF
    grp_id = var_groups%group_id(TRIM(group_name_plain), opt_lcheck=.FALSE.)

    ! If the group does not exist, create it.
    IF (grp_id == 0) THEN
      !
      ! increase dynamic groups array by one element
      CALL resize_arr_c1d(var_groups%name,1)
      !
      ! add new group
      var_groups%name(SIZE(var_groups%name)) = toupper(TRIM(group_name_plain))
      !
      ! return its group ID (including offset from static groups array)
      grp_id = var_groups%group_id(TRIM(group_name_plain))
    ENDIF

    !
    ! update in_group metainfo
    in_group_new(:) = groups()   ! initialization
    IF (PRESENT(opt_in_group)) THEN
      in_group_new(1:SIZE(opt_in_group)) = opt_in_group(:)
    ENDIF
    !
    IF (grp_id > MAX_GROUPS)  CALL finish(routine, TRIM(group_name))
    in_group_new(grp_id) = .TRUE.

  END SUBROUTINE t_var_groups_add


  !----------------------------------------------------------------------------------------
  !> Returns an alphabetically sorted list of all groups.
  ! 
  FUNCTION t_var_groups_alphabetical_list(var_groups)  RESULT(sorted_list)
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: sorted_list(:)
    CLASS(t_var_groups), INTENT(IN) :: var_groups
    INTEGER :: i

    ALLOCATE(CHARACTER(VARNAME_LEN) :: sorted_list(SIZE(var_groups%name)))
    DO i=1,SIZE(var_groups%name(:))
      sorted_list(i) = toupper(var_groups%name(i))
    END DO
    CALL quicksort(sorted_list)
  END FUNCTION t_var_groups_alphabetical_list


  !----------------------------------------------------------------------------------------
  !> Utility function with *a lot* of optional string parameters g1,
  !  g2, g3, g4, ...; mapping those onto a
  !  LOGICAL(DIMENSION=MAX_GROUPS) according to the "group_id"
  !  function.
  !
  FUNCTION groups_arg(g01, g02, g03, g04, g05, g06, g07, g08, g09, g10, g11, g12, g13)
    LOGICAL :: groups_arg(MAX_GROUPS)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: &
      &   g01, g02, g03, g04, g05, g06, g07, g08, g09, g10, g11, g12, g13

    groups_arg(:) = .FALSE.
    groups_arg(var_groups_dyn%group_id("ALL")) = .TRUE.
    IF (PRESENT(g01)) groups_arg(var_groups_dyn%group_id(g01)) = .TRUE.
    IF (PRESENT(g02)) groups_arg(var_groups_dyn%group_id(g02)) = .TRUE.
    IF (PRESENT(g03)) groups_arg(var_groups_dyn%group_id(g03)) = .TRUE.
    IF (PRESENT(g04)) groups_arg(var_groups_dyn%group_id(g04)) = .TRUE.
    IF (PRESENT(g05)) groups_arg(var_groups_dyn%group_id(g05)) = .TRUE.
    IF (PRESENT(g06)) groups_arg(var_groups_dyn%group_id(g06)) = .TRUE.
    IF (PRESENT(g07)) groups_arg(var_groups_dyn%group_id(g07)) = .TRUE.
    IF (PRESENT(g08)) groups_arg(var_groups_dyn%group_id(g08)) = .TRUE.
    IF (PRESENT(g09)) groups_arg(var_groups_dyn%group_id(g09)) = .TRUE.
    IF (PRESENT(g10)) groups_arg(var_groups_dyn%group_id(g10)) = .TRUE.
    IF (PRESENT(g11)) groups_arg(var_groups_dyn%group_id(g11)) = .TRUE.
    IF (PRESENT(g12)) groups_arg(var_groups_dyn%group_id(g12)) = .TRUE.
    IF (PRESENT(g13)) groups_arg(var_groups_dyn%group_id(g13)) = .TRUE.
  END FUNCTION groups_arg


  !----------------------------------------------------------------------------------------
  !> The same, but provide list of groups as one character vector of group names.
  !  Attention: the strings passed in group_list must be of length VARNAME_LEN !
  !
  FUNCTION groups_vec(group_list)
    LOGICAL :: groups_vec(MAX_GROUPS)
    CHARACTER(LEN=VARNAME_LEN), INTENT(IN) :: group_list(:)

    CHARACTER(*), PARAMETER :: routine = modname//"::groups_vec"
    INTEGER :: i, grp_id

    groups_vec(:) = .FALSE.
    groups_vec(var_groups_dyn%group_id("ALL")) = .TRUE.
    DO i=1,SIZE(group_list)
      IF (TRIM(group_list(i)) == "ALL") CYCLE
      grp_id = var_groups_dyn%group_id(TRIM(group_list(i)))
      IF (grp_id > MAX_GROUPS)  CALL finish(routine, TRIM(group_list(i)))
      groups_vec(grp_id) = .TRUE.
    END DO
  END FUNCTION groups_vec

END MODULE mo_var_groups
