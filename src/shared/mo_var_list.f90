MODULE mo_var_list

  USE mo_kind,             ONLY: dp
  USE mo_cf_convention,    ONLY: t_cf_var
  USE mo_grib2,            ONLY: t_grib2_var
  USE mo_var_metadata,     ONLY: t_var_metadata, t_union_vals
  USE mo_var_list_element, ONLY: t_var_list_element
  USE mo_linked_list,      ONLY: t_var_list, t_list_element, &
       &                         new_list, delete_list,      &
       &                         append_list_element,        &
       &                         find_list_element,          &
       &                         remove_list_element 
  USE mo_exception,        ONLY: message, message_text, finish

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: new_var_list              ! get a pointer to a new output var_list
  PUBLIC :: delete_var_list           ! delete an output var_list
  PUBLIC :: delete_var_lists          ! delete all output var_lists
  PUBLIC :: get_var_list              ! get a pointer to an existing output var_list
  PUBLIC :: set_var_list              ! set default parameters of an output var_list

  PUBLIC :: default_var_list_settings ! set default settings for a whole list

  PUBLIC :: var_lists                 ! vector of output var_lists
  PUBLIC :: nvar_lists                ! number of output var_lists defined so far

  PUBLIC :: add_var                   ! create/allocate a new var_list list entry
  PUBLIC :: get_var                   ! obtain reference to existing list entry

  INTERFACE add_var  ! create a new list entry
     MODULE PROCEDURE add_var_list_element_r4d
     MODULE PROCEDURE add_var_list_element_r3d 
     MODULE PROCEDURE add_var_list_element_r2d 
     MODULE PROCEDURE add_var_list_element_r1d 
     MODULE PROCEDURE add_var_list_element_i4d
     MODULE PROCEDURE add_var_list_element_i3d 
     MODULE PROCEDURE add_var_list_element_i2d 
     MODULE PROCEDURE add_var_list_element_i1d 
     MODULE PROCEDURE add_var_list_element_l4d
     MODULE PROCEDURE add_var_list_element_l3d 
     MODULE PROCEDURE add_var_list_element_l2d 
     MODULE PROCEDURE add_var_list_element_l1d 
  END INTERFACE
  
  INTERFACE get_var  ! obtain reference to a list entry
     MODULE PROCEDURE get_var_list_element_r4d
     MODULE PROCEDURE get_var_list_element_r3d
     MODULE PROCEDURE get_var_list_element_r2d
     MODULE PROCEDURE get_var_list_element_r1d
!     MODULE PROCEDURE get_var_list_element_i4d
!     MODULE PROCEDURE get_var_list_element_i3d
!     MODULE PROCEDURE get_var_list_element_i2d
!     MODULE PROCEDURE get_var_list_element_i1d
!     MODULE PROCEDURE get_var_list_element_l4d
!     MODULE PROCEDURE get_var_list_element_l3d
!     MODULE PROCEDURE get_var_list_element_l2d
!     MODULE PROCEDURE get_var_list_element_l1d
  END INTERFACE

  INTERFACE assign_if_present  ! purely internal
     MODULE PROCEDURE assign_if_present_character
     MODULE PROCEDURE assign_if_present_logical
     MODULE PROCEDURE assign_if_present_integer
     MODULE PROCEDURE assign_if_present_integers
     MODULE PROCEDURE assign_if_present_real
     MODULE PROCEDURE assign_if_present_cf
     MODULE PROCEDURE assign_if_present_grib2
     MODULE PROCEDURE assign_if_present_union_vals
  END INTERFACE

  INTEGER, PARAMETER           :: max_var_lists  = 64    ! max number of output var_lists
  INTEGER,                SAVE :: nvar_lists =  0    ! var_lists allocated so far
  TYPE(t_var_list), TARGET, SAVE :: var_lists(max_var_lists) ! memory buffer array

CONTAINS
  !------------------------------------------------------------------------------------------------
  !
  ! Create a new memory buffer / output var_list
  ! Get a pointer to the new var_list
  !
  SUBROUTINE new_var_list (this_list, name, output_type, restart_type, &
       post_suf, rest_suf, init_suf, lpost, lrestart, linitial)
    !
    TYPE(t_var_list),   POINTER, INTENT(out)          :: this_list ! anchor
    CHARACTER(len=*),          INTENT(in)           :: name         ! name of output var_list
    INTEGER,                   INTENT(in), OPTIONAL :: output_type  ! 'GRIB2' or 'NetCDF'
    INTEGER,                   INTENT(in), OPTIONAL :: restart_type ! 'GRIB2' or 'NetCDF'
    CHARACTER(len=*),          INTENT(in), OPTIONAL :: post_suf     ! suffix of output file
    CHARACTER(len=*),          INTENT(in), OPTIONAL :: rest_suf     ! suffix of restart file
    CHARACTER(len=*),          INTENT(in), OPTIONAL :: init_suf     ! suffix of initial file
    LOGICAL,                   INTENT(in), OPTIONAL :: lpost        ! write to  output file
    LOGICAL,                   INTENT(in), OPTIONAL :: lrestart     ! write to restart file
    LOGICAL,                   INTENT(in), OPTIONAL :: linitial     ! read from initial file
    !
    INTEGER :: i
    !
    CALL message('new_var_list','Adding new var_list '//name)
    !
    ! name must be unique
    !
    NULLIFY (this_list)
    !
    ! look if name exists already in list
    !
    IF (ANY(var_lists(1:nvar_lists)%name == name)) THEN
      CALL finish('new_list', 'output var_list '//name//' already used.')
    ENDIF
    !
    ! find next free entry
    !
    DO i = 1, nvar_lists    
      IF (var_lists(i)%name == '') THEN
        this_list => var_lists(i)
        EXIT
      ENDIF
    END DO
    !
    IF(.NOT. ASSOCIATED(this_list)) THEN
      nvar_lists = nvar_lists + 1
      IF (nvar_lists > max_var_lists) THEN
        CALL finish('new_list', 'var_lists container overflow, increase "max_var_lists"')
      ENDIF
      this_list => var_lists(nvar_lists)
    ENDIF
    !
    CALL new_list (this_list)
    !
    ! set default list characteristics
    !
    this_list%name     = name
    this_list%post_suf = '_'//name
    this_list%rest_suf = this_list%post_suf
    this_list%init_suf = this_list%post_suf
    !
    ! set non-default list characteristics
    !
    CALL assign_if_present(this_list%output_type,  output_type)
    CALL assign_if_present(this_list%restart_type, restart_type)
    CALL assign_if_present(this_list%post_suf,     post_suf) 
    CALL assign_if_present(this_list%rest_suf,     rest_suf) 
    CALL assign_if_present(this_list%init_suf,     init_suf) 
    CALL assign_if_present(this_list%lpost,        lpost)
    CALL assign_if_present(this_list%lrestart,     lrestart)
    CALL assign_if_present(this_list%linitial,     linitial)
    !
  END SUBROUTINE new_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! Get a reference to a memory buffer/output var_list
  !
  SUBROUTINE get_var_list (this_list, name)
    !
    TYPE(t_var_list),  POINTER, INTENT(out) :: this_list ! pointer
    CHARACTER(len=*),         INTENT(in)  :: name      ! name of output var_list
    !
    INTEGER :: i
    !
    NULLIFY (this_list)
    !
    DO i = 1, nvar_lists
      IF (var_lists(i)%name == name) THEN
        this_list => var_lists(i)
        EXIT
      ENDIF
    END DO
    !
  END SUBROUTINE get_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! Change parameters of an already existent output var_list
  !
  SUBROUTINE set_var_list (this_list, output_type, restart_type, &
       post_suf, rest_suf, init_suf, lpost, lrestart, linitial)
    !
    TYPE(t_var_list),   INTENT(inout)        :: this_list   ! output var_list to change
    INTEGER,          INTENT(in), OPTIONAL :: output_type    ! 'GRIB' or 'NetCDF'
    INTEGER,          INTENT(in), OPTIONAL :: restart_type   ! 'GRIB' or 'NetCDF'
    CHARACTER(len=*), INTENT(in), OPTIONAL :: post_suf       ! suffix of output  file
    CHARACTER(len=*), INTENT(in), OPTIONAL :: rest_suf       ! suffix of restart file
    CHARACTER(len=*), INTENT(in), OPTIONAL :: init_suf       ! suffix of initial file
    LOGICAL,          INTENT(in), OPTIONAL :: lpost          ! in standard output file
    LOGICAL,          INTENT(in), OPTIONAL :: lrestart       ! in standard restartfile
    LOGICAL,          INTENT(in), OPTIONAL :: linitial       ! in standard initialfile
    !
    CALL assign_if_present(this_list%output_type,  output_type)
    CALL assign_if_present(this_list%restart_type, restart_type)
    CALL assign_if_present(this_list%post_suf,     post_suf) 
    CALL assign_if_present(this_list%rest_suf,     rest_suf) 
    CALL assign_if_present(this_list%init_suf,     init_suf) 
    CALL assign_if_present(this_list%lpost,        lpost)
    CALL assign_if_present(this_list%lrestart,     lrestart)
    CALL assign_if_present(this_list%linitial,     linitial)
    !
  END SUBROUTINE set_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! Delete an output var_list, nullify the associated pointer
  !
  SUBROUTINE delete_var_list(this_list)
    !
    TYPE(t_var_list), POINTER :: this_list
    !
    IF (ASSOCIATED(this_list)) THEN
      CALL delete_list(this_list)
      NULLIFY(this_list)
    ENDIF
    !
  END SUBROUTINE delete_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! Delete all output var_lists 
  !
  SUBROUTINE delete_var_lists
    !
    TYPE(t_var_list), POINTER :: this_list
    !
    INTEGER :: i
    !
    DO i = 1, nvar_lists
      this_list => var_lists(i)
      CALL delete_var_list (this_list)
    END DO
    !
  END SUBROUTINE delete_var_lists
  !------------------------------------------------------------------------------------------------
  !
  ! Get a copy of the metadata concerning a var_list element
  !
  SUBROUTINE get_var_list_element_info (this_list, name, info)
    !
    TYPE(t_var_list),     INTENT(in)  :: this_list ! list
    CHARACTER(len=*),   INTENT(in)  :: name         ! name of variable
    TYPE(t_var_metadata), INTENT(out) :: info         ! variable meta data
    !
    TYPE(t_list_element), POINTER :: element
    ! 
    element => find_list_element (this_list, name)
    IF (ASSOCIATED (element)) THEN
      info = element%field%info
    ENDIF
    !
  END SUBROUTINE get_var_list_element_info
  !------------------------------------------------------------------------------------------------
  !
  ! Set default meta data of output var_list
  !
  SUBROUTINE default_var_list_settings (this_list, &
       name, cf, grib2, lrestart, lpost, laccu, resetval, lmiss, missval, &
       lrestart_cont, grid_representation)
    !
    TYPE(t_var_list),   INTENT(inout)         :: this_list        ! output var_list
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: name          ! variable name
    TYPE(t_cf_var),     INTENT(in), OPTIONAL :: cf         
    TYPE(t_grib2_var),  INTENT(in), OPTIONAL :: grib2
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart      ! restart file flag
    LOGICAL,            INTENT(in), OPTIONAL :: lpost         ! into output var_list
    LOGICAL,            INTENT(in), OPTIONAL :: laccu         ! accumulation flag
    TYPE(t_union_vals), INTENT(in), OPTIONAL :: resetval      ! reset value
    LOGICAL,            INTENT(in), OPTIONAL :: lmiss         ! missing value flag
    TYPE(t_union_vals), INTENT(in), OPTIONAL :: missval       ! missing value
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart_cont ! continue on restart
    INTEGER,            INTENT(in), OPTIONAL :: grid_representation
    !
    ! more complicated ... haven't got default_info anymore ...
    !    CALL set_var_list_metadata (this_list%default_info                                &
    !       name, cf, grib2, lrestart, lpost, laccu, reset, lmiss, missval, &
    !       lrestart_cont, grid_representation, verbose=.FALSE.)
    !
  END SUBROUTINE default_var_list_settings
  !------------------------------------------------------------------------------------------------
  !
  ! Change metadata of a var_list element
  !
  SUBROUTINE set_var_list_element_metadata (this_list,                 &
       name, cf, grib2, lrestart, lpost, laccu, resetval, lmiss, missval, &
       lrestart_cont, initval, grid_representation, verbose)
    !
    TYPE(t_var_list),   INTENT(inout)        :: this_list  ! memory info struct.
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: name          ! variable name
    TYPE(t_cf_var),     INTENT(in), OPTIONAL :: cf         
    TYPE(t_grib2_var),  INTENT(in), OPTIONAL :: grib2
    LOGICAL,            INTENT(in), OPTIONAL :: lpost         ! into output var_list
    LOGICAL,            INTENT(in), OPTIONAL :: laccu         ! accumulation flag
    TYPE(t_union_vals), INTENT(in), OPTIONAL :: resetval      ! reset value
    LOGICAL,            INTENT(in), OPTIONAL :: lmiss         ! missing value flag
    TYPE(t_union_vals), INTENT(in), OPTIONAL :: missval       ! missing value
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart_cont ! continue on restart
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart      ! restart file flag
    TYPE(t_union_vals), INTENT(in), OPTIONAL :: initval       ! value if var not available
    INTEGER,            INTENT(in), OPTIONAL :: grid_representation
    LOGICAL,            INTENT(in), OPTIONAL :: verbose
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    !
    IF (ASSOCIATED (element)) THEN
      CALL set_var_metadata (element%field%info,                                     &
       name=name, grid_representation=grid_representation, cf=cf, grib2=grib2,       &
       lpost=lpost, lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
       laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,                 &
       verbose=verbose)
    ENDIF
    !
  END SUBROUTINE set_var_list_element_metadata
  !------------------------------------------------------------------------------------------------
  !
  ! Set parameters of list element already created
  ! (private routine within this module)
  !
  ! Set each parameter in data type var_metadata if the respective 
  ! optional parameter is present.
  !
  SUBROUTINE set_var_metadata (info,                  &
         name, grid_representation, cf, grib2, ldims, &
         lpost, lrestart, lrestart_cont, initval,     &
         laccu, resetval, lmiss, missval, verbose)
    !
    TYPE(t_var_metadata), INTENT(inout)        :: info          ! memory info struct.
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: name          ! variable name
    INTEGER,            INTENT(in), OPTIONAL :: grid_representation
    TYPE(t_cf_var),       INTENT(in), OPTIONAL :: cf         
    TYPE(t_grib2_var),    INTENT(in), OPTIONAL :: grib2
    INTEGER,            INTENT(in), OPTIONAL :: ldims(:)      ! dimensions to allocate  
    LOGICAL,            INTENT(in), OPTIONAL :: lpost         ! into output var_list
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart      ! restart file flag
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart_cont ! continue on restart
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: initval       ! value if var not available
    LOGICAL,            INTENT(in), OPTIONAL :: laccu         ! accumulation flag
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: resetval      ! reset value
    LOGICAL,            INTENT(in), OPTIONAL :: lmiss         ! missing value flag
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: missval       ! missing value
    LOGICAL,            INTENT(in), OPTIONAL :: verbose
    !
    LOGICAL :: lverbose
    !
    ! set flags from optional parameters
    !
    lverbose = .FALSE.; IF(PRESENT(verbose)) lverbose = verbose
    !
    ! set components describing the 'Content of the field'
    !
    CALL assign_if_present (info%name,  name)
    CALL assign_if_present (info%cf,    cf)
    CALL assign_if_present (info%grib2, grib2)
    !
    ! set grid type
    !
    CALL assign_if_present (info%grid_representation, grid_representation)
    !
    ! set flags concerning I/O
    !
    CALL assign_if_present (info%lpost,         lpost)
    CALL assign_if_present (info%resetval,      resetval)
    CALL assign_if_present (info%laccu,         laccu)
    CALL assign_if_present (info%lmiss,         lmiss)
    CALL assign_if_present (info%missval,       missval)
    CALL assign_if_present (info%lrestart,      lrestart)
    CALL assign_if_present (info%lrestart_cont, lrestart_cont)
    CALL assign_if_present (info%initval,       initval)
    !
    ! printout (optional)
    !
!LK    IF (lverbose) CALL print_var_metadata (info)
    !
  END SUBROUTINE set_var_metadata
  !------------------------------------------------------------------------------------------------
  !
  ! Create a list new entry
  !
  ! Specific routines for pointers of different rank
  !
  !================================================================================================ 
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 4d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_r4d(this_list, name, ptr, &
       grid_representation, cf, grib2, ldims, lpost,       &
       lrestart, lrestart_cont, initval, laccu, resetval,  &
       lmiss, missval, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(dp),             POINTER              :: ptr(:,:,:,:)        ! reference to field
    INTEGER,              INTENT(in)           :: grid_representation ! grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(4)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    TYPE(t_union_vals),        INTENT(in), OPTIONAL :: initval             ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    TYPE(t_union_vals),        INTENT(in), OPTIONAL :: resetval            ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    TYPE(t_union_vals),        INTENT(in), OPTIONAL :: missval             ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(dp),             POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    !
    ! and set meta data
    !
    IF (PRESENT(p4)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, grid_representation=grid_representation,              & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = new_list_element%field%info%used_dimensions(4)
      ALLOCATE(new_list_element%field%r_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r4d','allocation of array '//TRIM(name)//'failed')
      ENDIF
      NULLIFY(new_list_element%field%i_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%memory_used = this_list%memory_used+8*SIZE(new_list_element%field%r_ptr)
    ELSE
      new_list_element%field%r_ptr => p4
    ENDIF
    ptr => new_list_element%field%r_ptr(:,:,:,:)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%r_ptr = new_list_element%field%info%missval%rval
    ELSE
      new_list_element%field%r_ptr = 0.0_dp
    END IF
    !
  END SUBROUTINE add_var_list_element_r4d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 3d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_r3d(this_list, name, ptr, &
       grid_representation, cf, grib2, ldims, lpost,       &
       lrestart, lrestart_cont, initval, laccu, resetval,  &
       lmiss, missval, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),   INTENT(in)           :: name                ! name of variable
    REAL(dp),           POINTER              :: ptr(:,:,:)          ! reference to field
    INTEGER,            INTENT(in)           :: grid_representation ! grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,            INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions
    LOGICAL,            INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: initval             ! value if var not available
    LOGICAL,            INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: resetval            ! reset value (after accumulation)
    LOGICAL,            INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: missval             ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(dp),           POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,            INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    !
    ! and set meta data
    !
    IF (PRESENT(p4)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, grid_representation=grid_representation,              & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = 1
      ALLOCATE (new_list_element%field%r_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r3d','allocation of array '//TRIM(name)//'failed')
      ENDIF
      NULLIFY(new_list_element%field%i_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%memory_used = this_list%memory_used+8*SIZE(new_list_element%field%r_ptr)
    ELSE
      new_list_element%field%r_ptr => p4
    ENDIF
    ptr => new_list_element%field%r_ptr(:,:,:,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%r_ptr = new_list_element%field%info%missval%rval
    ELSE
      new_list_element%field%r_ptr = 0.0_dp
    END IF
    !
  END SUBROUTINE add_var_list_element_r3d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 2d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_r2d(this_list, name, ptr, &
       grid_representation, cf, grib2, ldims, lpost,       &
       lrestart, lrestart_cont, initval, laccu, resetval,  &
       lmiss, missval, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),   INTENT(in)           :: name                ! name of variable
    REAL(dp),           POINTER              :: ptr(:,:)            ! reference to field
    INTEGER,            INTENT(in)           :: grid_representation ! grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,            INTENT(in), OPTIONAL :: ldims(2)            ! local dimensions
    LOGICAL,            INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: initval             ! value if var not available
    LOGICAL,            INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: resetval            ! reset value (after accumulation)
    LOGICAL,            INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: missval             ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(dp),           POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,            INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    !
    ! and set meta data
    !
    IF (PRESENT(p4)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, grid_representation=grid_representation,              & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = 1
      idims(4) = 1
      ALLOCATE (new_list_element%field%r_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r2d','allocation of array '//TRIM(name)//'failed')
      ENDIF
      NULLIFY(new_list_element%field%i_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%memory_used = this_list%memory_used+8*SIZE(new_list_element%field%r_ptr)
    ELSE
      new_list_element%field%r_ptr => p4
    ENDIF
    ptr => new_list_element%field%r_ptr(:,:,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%r_ptr = new_list_element%field%info%missval%rval
    ELSE
      new_list_element%field%r_ptr = 0.0_dp
    END IF
    !
  END SUBROUTINE add_var_list_element_r2d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 1d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_r1d(this_list, name, ptr, &
       grid_representation, cf, grib2, ldims, lpost,       &
       lrestart, lrestart_cont, initval, laccu, resetval,  &
       lmiss, missval, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),   INTENT(in)           :: name                ! name of variable
    REAL(dp),           POINTER              :: ptr(:)              ! reference to field
    INTEGER,            INTENT(in)           :: grid_representation ! grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,            INTENT(in), OPTIONAL :: ldims(1)            ! local dimensions
    LOGICAL,            INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: initval             ! value if var not available
    LOGICAL,            INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: resetval            ! reset value (after accumulation)
    LOGICAL,            INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: missval             ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(dp),           POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,            INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    !
    ! and set meta data
    !
    IF (PRESENT(p4)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, grid_representation=grid_representation,              & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = 1
      idims(3) = 1
      idims(4) = 1
      ALLOCATE (new_list_element%field%r_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r1d','allocation of array '//TRIM(name)//'failed')
      ENDIF
      NULLIFY(new_list_element%field%i_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%memory_used = this_list%memory_used+8*SIZE(new_list_element%field%r_ptr)
    ELSE
      new_list_element%field%r_ptr => p4
    ENDIF
    ptr => new_list_element%field%r_ptr(:,1,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%r_ptr = new_list_element%field%info%missval%rval
    ELSE
      new_list_element%field%r_ptr = 0.0_dp
    END IF
    !
  END SUBROUTINE add_var_list_element_r1d
  !
  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 4d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_i4d(this_list, name, ptr, &
       grid_representation, cf, grib2, ldims, lpost,       &
       lrestart, lrestart_cont, initval, laccu, resetval,  &
       lmiss, missval, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),   INTENT(in)           :: name                ! name of variable
    INTEGER,            POINTER              :: ptr(:,:,:,:)        ! reference to field
    INTEGER,            INTENT(in)           :: grid_representation ! grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,            INTENT(in), OPTIONAL :: ldims(4)            ! local dimensions
    LOGICAL,            INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: initval             ! value if var not available
    LOGICAL,            INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: resetval            ! reset value (after accumulation)
    LOGICAL,            INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: missval             ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    INTEGER,            POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,            INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    !
    ! and set meta data
    !
    IF (PRESENT(p4)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, grid_representation=grid_representation,              & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = new_list_element%field%info%used_dimensions(4)
      ALLOCATE (new_list_element%field%i_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_i4d','allocation of array '//TRIM(name)//'failed')
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%memory_used = this_list%memory_used+4*SIZE(new_list_element%field%i_ptr)
    ELSE
      new_list_element%field%i_ptr => p4
    ENDIF
    ptr => new_list_element%field%i_ptr(:,:,:,:)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%i_ptr = new_list_element%field%info%missval%ival
    ELSE
      new_list_element%field%i_ptr = 0
    END IF
    !
  END SUBROUTINE add_var_list_element_i4d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 3d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_i3d(this_list, name, ptr, &
       grid_representation, cf, grib2, ldims, lpost,       &
       lrestart, lrestart_cont, initval, laccu, resetval,  &
       lmiss, missval, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),   INTENT(in)           :: name                ! name of variable
    INTEGER,            POINTER              :: ptr(:,:,:)          ! reference to field
    INTEGER,            INTENT(in)           :: grid_representation ! grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,            INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions
    LOGICAL,            INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: initval             ! value if var not available
    LOGICAL,            INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: resetval            ! reset value (after accumulation)
    LOGICAL,            INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: missval             ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    INTEGER,            POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,            INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    !
    ! and set meta data
    !
    IF (PRESENT(p4)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, grid_representation=grid_representation,              & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = 1
      ALLOCATE (new_list_element%field%i_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r3d','allocation of array '//TRIM(name)//'failed')
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%memory_used = this_list%memory_used+4*SIZE(new_list_element%field%i_ptr)
    ELSE
      new_list_element%field%i_ptr => p4
    ENDIF
    ptr => new_list_element%field%i_ptr(:,:,:,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%i_ptr = new_list_element%field%info%missval%ival
    ELSE
      new_list_element%field%i_ptr = 0
    END IF
    !
  END SUBROUTINE add_var_list_element_i3d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 2d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_i2d(this_list, name, ptr, &
       grid_representation, cf, grib2, ldims, lpost,       &
       lrestart, lrestart_cont, initval, laccu, resetval,  &
       lmiss, missval, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    INTEGER,              POINTER              :: ptr(:,:)            ! reference to field
    INTEGER,              INTENT(in)           :: grid_representation ! grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(2)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: initval             ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: resetval            ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: missval             ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    INTEGER,              POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    !
    ! and set meta data
    !
    IF (PRESENT(p4)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, grid_representation=grid_representation,              & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = 1
      idims(4) = 1
      ALLOCATE (new_list_element%field%i_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r2d','allocation of array '//TRIM(name)//'failed')
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%memory_used = this_list%memory_used+4*SIZE(new_list_element%field%i_ptr)
    ELSE
      new_list_element%field%i_ptr => p4
    ENDIF
    ptr => new_list_element%field%i_ptr(:,:,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%i_ptr = new_list_element%field%info%missval%ival
    ELSE
      new_list_element%field%i_ptr = 0
    END IF
    !
  END SUBROUTINE add_var_list_element_i2d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 1d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_i1d(this_list, name, ptr, &
       grid_representation, cf, grib2, ldims, lpost,       &
       lrestart, lrestart_cont, initval, laccu, resetval,  &
       lmiss, missval, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    INTEGER,              POINTER              :: ptr(:)              ! reference to field
    INTEGER,              INTENT(in)           :: grid_representation ! grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(1)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: initval             ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: resetval            ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: missval             ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    INTEGER,              POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    !
    ! and set meta data
    !
    IF (PRESENT(p4)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, grid_representation=grid_representation,              & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = 1
      idims(3) = 1
      idims(4) = 1
      ALLOCATE (new_list_element%field%i_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r1d','allocation of array '//TRIM(name)//'failed')
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%memory_used = this_list%memory_used+4*SIZE(new_list_element%field%i_ptr)
    ELSE
      new_list_element%field%i_ptr => p4
    ENDIF
    ptr => new_list_element%field%i_ptr(:,1,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%i_ptr = new_list_element%field%info%missval%ival
    ELSE
      new_list_element%field%i_ptr = 0
    END IF
    !
  END SUBROUTINE add_var_list_element_i1d
  !
  !================================================================================================
  ! LOGICAL SECTION -------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 4d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_l4d(this_list, name, ptr, &
       grid_representation, cf, grib2, ldims, lpost,       &
       lrestart, lrestart_cont, initval, laccu, resetval,  &
       lmiss, missval, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    LOGICAL,              POINTER              :: ptr(:,:,:,:)        ! reference to field
    INTEGER,              INTENT(in)           :: grid_representation ! grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(4)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: initval             ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: resetval            ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: missval             ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    LOGICAL,              POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    !
    ! and set meta data
    !
    IF (PRESENT(p4)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, grid_representation=grid_representation,              & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = new_list_element%field%info%used_dimensions(4)
      ALLOCATE (new_list_element%field%l_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r4d','allocation of array '//TRIM(name)//'failed')
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%i_ptr)      
      this_list%memory_used = this_list%memory_used+4*SIZE(new_list_element%field%l_ptr)
    ELSE
      new_list_element%field%l_ptr => p4
    ENDIF
    ptr => new_list_element%field%l_ptr(:,:,:,:)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%l_ptr = new_list_element%field%info%missval%lval
    ELSE
      new_list_element%field%l_ptr = .FALSE.
    END IF
    !
  END SUBROUTINE add_var_list_element_l4d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 3d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_l3d(this_list, name, ptr, &
       grid_representation, cf, grib2, ldims, lpost,       &
       lrestart, lrestart_cont, initval, laccu, resetval,  &
       lmiss, missval, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    LOGICAL,              POINTER              :: ptr(:,:,:)          ! reference to field
    INTEGER,              INTENT(in)           :: grid_representation ! grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: initval             ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: resetval            ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: missval             ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    LOGICAL,              POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    !
    ! and set meta data
    !
    IF (PRESENT(p4)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, grid_representation=grid_representation,              & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = 1
      ALLOCATE (new_list_element%field%l_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r3d','allocation of array '//TRIM(name)//'failed')
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%i_ptr)      
      this_list%memory_used = this_list%memory_used+4*SIZE(new_list_element%field%l_ptr)
    ELSE
      new_list_element%field%l_ptr => p4
    ENDIF
    ptr => new_list_element%field%l_ptr(:,:,:,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%l_ptr = new_list_element%field%info%missval%lval
    ELSE
      new_list_element%field%l_ptr = .FALSE.
    END IF
    !
  END SUBROUTINE add_var_list_element_l3d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 2d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_l2d(this_list, name, ptr, &
       grid_representation, cf, grib2, ldims, lpost,       &
       lrestart, lrestart_cont, initval, laccu, resetval,  &
       lmiss, missval, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    LOGICAL,              POINTER              :: ptr(:,:)            ! reference to field
    INTEGER,              INTENT(in)           :: grid_representation ! grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(2)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: initval             ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: resetval            ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: missval             ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    LOGICAL,              POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    !
    ! and set meta data
    !
    IF (PRESENT(p4)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, grid_representation=grid_representation,              & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = 1
      idims(4) = 1
      ALLOCATE (new_list_element%field%l_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r2d','allocation of array '//TRIM(name)//'failed')
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%i_ptr)      
      this_list%memory_used = this_list%memory_used+4*SIZE(new_list_element%field%l_ptr)
    ELSE
      new_list_element%field%l_ptr => p4
    ENDIF
    ptr => new_list_element%field%l_ptr(:,:,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%l_ptr = new_list_element%field%info%missval%lval
    ELSE
      new_list_element%field%l_ptr = .FALSE.
    END IF
    !
  END SUBROUTINE add_var_list_element_l2d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 1d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_l1d(this_list, name, ptr, &
       grid_representation, cf, grib2, ldims, lpost,       &
       lrestart, lrestart_cont, initval, laccu, resetval,  &
       lmiss, missval, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    LOGICAL,              POINTER              :: ptr(:)              ! reference to field
    INTEGER,              INTENT(in)           :: grid_representation ! grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(1)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: initval             ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: resetval            ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: missval             ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    LOGICAL,              POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    !
    ! and set meta data
    !
    IF (PRESENT(p4)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, grid_representation=grid_representation,              & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = 1
      idims(3) = 1
      idims(4) = 1
      ALLOCATE (new_list_element%field%l_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r1d','allocation of array '//TRIM(name)//'failed')
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%i_ptr)      
      this_list%memory_used = this_list%memory_used+4*SIZE(new_list_element%field%l_ptr)
    ELSE
      new_list_element%field%l_ptr => p4
    ENDIF
    ptr => new_list_element%field%l_ptr(:,1,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%l_ptr = new_list_element%field%info%missval%lval
    ELSE
      new_list_element%field%l_ptr = .FALSE.
    END IF
    !
  END SUBROUTINE add_var_list_element_l1d
  !================================================================================================
  !------------------------------------------------------------------------------------------------
  !
  ! remove one element from the list
  ! the element is identified by its name
  !
  SUBROUTINE remove_var_list_element (this_list, name)
    TYPE(t_var_list),   INTENT(inout) :: this_list
    CHARACTER(len=*), INTENT(in)    :: name   
    TYPE(t_list_element), POINTER :: ptr

    IF (this_list%first_list_element%field%info%name == name) THEN
      CALL remove_list_element (this_list, this_list%first_list_element)
      RETURN
    ELSE
      ptr => this_list%first_list_element
      DO
        IF (.NOT.ASSOCIATED (ptr%next_list_element)) EXIT
        IF (ptr%next_list_element%field%info%name == name) THEN
          CALL remove_list_element (this_list, ptr%next_list_element)
          EXIT
        ENDIF
        ptr => ptr%next_list_element
      END DO
    ENDIF
    !
  END SUBROUTINE remove_var_list_element
  !------------------------------------------------------------------------------------------------
  !
  ! Get element of a list, specific routines for pointers of different rank
  !
  ! obtain pointer to 4d-field
  !
  SUBROUTINE get_var_list_element_r4d (this_list, name, ptr)
    TYPE(t_var_list),   INTENT(in) :: this_list ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    REAL(dp),         POINTER    :: ptr(:,:,:,:) ! reference to allocated field

    TYPE(t_list_element), POINTER :: element

    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%r_ptr
    
  END SUBROUTINE get_var_list_element_r4d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 3d-field
  !
  SUBROUTINE get_var_list_element_r3d (this_list, name, ptr)

    TYPE(t_var_list) ,INTENT(in) :: this_list       ! list
    CHARACTER (*)   ,INTENT(in) :: name         ! name of variable
    REAL(dp)        ,POINTER    :: ptr(:,:,:)   ! reference to allocated field

    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,:,:,1)

  END SUBROUTINE get_var_list_element_r3d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 2d-field
  !
  SUBROUTINE get_var_list_element_r2d (this_list, name, ptr)

    TYPE(t_var_list)   ,INTENT(in) :: this_list     ! list
    CHARACTER (*)     ,INTENT(in) :: name       ! name of variable
    REAL(dp)          ,POINTER    :: ptr(:,:)   ! reference to allocated field

    TYPE(t_list_element), POINTER :: element

    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,:,1,1)

  END SUBROUTINE get_var_list_element_r2d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 1d-field
  !
  SUBROUTINE get_var_list_element_r1d (this_list, name, ptr)

    TYPE(t_var_list) ,INTENT(in) :: this_list    ! list
    CHARACTER (*)   ,INTENT(in) :: name      ! name of variable
    REAL(dp)        ,POINTER    :: ptr(:)    ! reference to allocated field

    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,1,1,1)

  END SUBROUTINE get_var_list_element_r1d
  !------------------------------------------------------------------------------------------------
  !
  ! Print routines for control output and debuggung
  !
  SUBROUTINE print_memory_use (this_list)

    TYPE(t_var_list) ,INTENT(in) :: this_list ! list

    WRITE (message_text,'(a16,a,a,i10,a,i4,a)')                 &
         TRIM(this_list%name), '-buffer: ',                        &
         'Memory in use: ', this_list%memory_used/1024, ' kb in ', &
         this_list%list_elements, ' fields.'
    CALL message('',message_text)

  END SUBROUTINE print_memory_use
  !------------------------------------------------------------------------------------------------
  !
  ! print current memory table 
  !
  SUBROUTINE print_memory_table (this_list)
    TYPE(t_var_list),  INTENT(in) :: this_list ! list
    !
    CALL message('','')
    CALL message('','')
    CALL message('','Status of base memory:')    
    CALL message('','')
    !
!LK    CALL print_linked_list (this_list)
    !    
  END SUBROUTINE print_memory_table
  !------------------------------------------------------------------------------------------------
  !
  ! print current stat table 
  !
  SUBROUTINE print_sinfo (this_list)
    TYPE(t_var_list),  INTENT(in) :: this_list
    WRITE (message_text,'(a16,a)') TRIM(this_list%name), '-buffer: '
    CALL message('',message_text)
    CALL message('','')    
    CALL message('','')
    CALL message('','Statistic of base memory:')
    CALL message('','')
    !
!LK    CALL print_sinfo_list (this_list)
    !
  END SUBROUTINE print_sinfo
  !------------------------------------------------------------------------------------------------
  !
  ! add supplementary fields (eg. geopotential, surface pressure, gridbox area)
  !
  SUBROUTINE add_var_list_reference (to_var_list, name, from_var_list, lpost, kprec)
    TYPE(t_var_list)  ,INTENT(inout)            :: to_var_list
    CHARACTER(len=*) ,INTENT(in)              :: name
    CHARACTER(len=*) ,INTENT(in)    ,OPTIONAL :: from_var_list
    LOGICAL          ,INTENT(in)    ,OPTIONAL :: lpost
    INTEGER          ,INTENT(in)    ,OPTIONAL :: kprec
    !
    TYPE(t_var_list_element)  ,POINTER :: source
    TYPE(t_list_element) ,POINTER :: new_list_element
    !
    CALL locate (source ,name ,from_var_list)
    IF(ASSOCIATED(source)) THEN
!      IF (source%info%gribcode > 0) THEN
!        CALL append_list_element (to_var_list, new_list_element, &
!             code=source%info%gribcode)
!      ELSE
        CALL append_list_element (to_var_list, new_list_element)
!      ENDIF
      new_list_element%field                = source
      new_list_element%field%info%allocated = .FALSE.
      new_list_element%field%info%lrestart  = .FALSE.
      IF (PRESENT(lpost)) &
           new_list_element%field%info%lpost = lpost
!      IF (PRESENT(kprec)) &
!           new_list_element%field%info%gribbits = kprec
    ENDIF
    !
  CONTAINS
    !----------------------------------------------------------------------------------------------
    !
    ! find an entry
    !
    SUBROUTINE locate (element, name, in_var_list)
      TYPE(t_var_list_element), POINTER        :: element
      CHARACTER(len=*), INTENT(in)           :: name
      CHARACTER(len=*), INTENT(in), OPTIONAL :: in_var_list
      !
      INTEGER                     :: i
      TYPE(t_list_element), POINTER :: link
      !
      NULLIFY (element)
      !
      DO i = 1, nvar_lists
        IF (PRESENT(in_var_list)) THEN
          IF (in_var_list /= var_lists(i)%name) CYCLE
        ENDIF
        link => find_list_element (var_lists(i), name)
        IF (ASSOCIATED(link)) THEN
          element => link%field
          EXIT
        ENDIF
      END DO
      !
    END SUBROUTINE locate
    !
  END SUBROUTINE add_var_list_reference
  !------------------------------------------------------------------------------------------------
  !
  ! private routines to assign values if actual parameters are present
  !
  SUBROUTINE assign_if_present_character (y,x)
    CHARACTER(len=*), INTENT(inout)        :: y
    CHARACTER(len=*), INTENT(in) ,OPTIONAL :: x
    IF (.NOT. PRESENT(x)) RETURN
    IF ( x == ' ' )       RETURN      
    y = x
  END SUBROUTINE assign_if_present_character
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_logical (y,x)
    LOGICAL, INTENT(inout)        :: y
    LOGICAL, INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_logical
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_integer (y,x)
    INTEGER, INTENT(inout)        :: y
    INTEGER, INTENT(in) ,OPTIONAL :: x
    IF (.NOT. PRESENT(x)) RETURN
    IF ( x == -HUGE(x)  ) RETURN
    y = x
  END SUBROUTINE assign_if_present_integer
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_integers (y,x)
    INTEGER, INTENT(inout)        :: y (:)
    INTEGER, INTENT(in) ,OPTIONAL :: x (:)
    INTEGER :: n
    IF (PRESENT(x)) THEN
      n = MIN(SIZE(x), SIZE(y))
      y(1:n) = x(1:n)
    ENDIF
  END SUBROUTINE assign_if_present_integers
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_real (y,x)
    REAL(dp), INTENT(inout)        :: y
    REAL(dp), INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    IF ( x == -HUGE(x) ) RETURN
    y = x
  END SUBROUTINE assign_if_present_real
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_cf (y,x)
    TYPE(t_cf_var), INTENT(inout)        :: y
    TYPE(t_cf_var), INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    y = x
  END SUBROUTINE assign_if_present_cf
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_grib2 (y,x)
    TYPE(t_grib2_var), INTENT(inout)        :: y
    TYPE(t_grib2_var) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    y = x
  END SUBROUTINE assign_if_present_grib2
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_union_vals (y,x)
    TYPE(t_union_vals), INTENT(inout)        :: y
    TYPE(t_union_vals) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    y = x
  END SUBROUTINE assign_if_present_union_vals
  !------------------------------------------------------------------------------------------------
END MODULE mo_var_list
