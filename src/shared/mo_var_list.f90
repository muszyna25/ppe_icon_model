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
       &                         delete_list_element 
  USE mo_exception,        ONLY: message, message_text, finish

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: new_var_list              ! get a pointer to a new output var_list
  PUBLIC :: delete_var_list           ! delete an output var_list
  PUBLIC :: delete_var_lists          ! delete all output var_lists
  PUBLIC :: get_var_list              ! get a pointer to an existing output var_list
  PUBLIC :: set_var_list              ! set default parameters of an output var_list
  PUBLIC :: print_var_list
  PUBLIC :: print_memory_use

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
  END INTERFACE add_var
  
  INTERFACE get_var  ! obtain reference to a list entry
    MODULE PROCEDURE get_var_list_element_r4d
    MODULE PROCEDURE get_var_list_element_r3d
    MODULE PROCEDURE get_var_list_element_r2d
    MODULE PROCEDURE get_var_list_element_r1d
    MODULE PROCEDURE get_var_list_element_i4d
    MODULE PROCEDURE get_var_list_element_i3d
    MODULE PROCEDURE get_var_list_element_i2d
    MODULE PROCEDURE get_var_list_element_i1d
    MODULE PROCEDURE get_var_list_element_l4d
    MODULE PROCEDURE get_var_list_element_l3d
    MODULE PROCEDURE get_var_list_element_l2d
    MODULE PROCEDURE get_var_list_element_l1d
  END INTERFACE get_var
  
  INTERFACE assign_if_present  ! purely internal
    MODULE PROCEDURE assign_if_present_character
    MODULE PROCEDURE assign_if_present_logical
    MODULE PROCEDURE assign_if_present_integer
    MODULE PROCEDURE assign_if_present_integers
    MODULE PROCEDURE assign_if_present_real
    MODULE PROCEDURE assign_if_present_cf
    MODULE PROCEDURE assign_if_present_grib2
    MODULE PROCEDURE assign_if_present_union
  END INTERFACE assign_if_present
  
  INTEGER, PARAMETER             :: max_var_lists  = 128      ! max number of output var_lists
  INTEGER,                  SAVE :: nvar_lists =  0           ! var_lists allocated so far
  TYPE(t_var_list), TARGET, SAVE :: var_lists(max_var_lists)  ! memory buffer array
  
CONTAINS
  !------------------------------------------------------------------------------------------------
  !
  ! Create a new memory buffer / output var_list
  ! Get a pointer to the new var_list
  !
  SUBROUTINE new_var_list (this_list, name, output_type, restart_type,    &
       &                   post_suf, rest_suf, init_suf, lpost, lrestart, &
       &                   linitial)
    !
    TYPE(t_var_list), POINTER              :: this_list    ! anchor
    CHARACTER(len=*), INTENT(in)           :: name         ! name of output var_list
    INTEGER,          INTENT(in), OPTIONAL :: output_type  ! 'GRIB2' or 'NetCDF'
    INTEGER,          INTENT(in), OPTIONAL :: restart_type ! 'GRIB2' or 'NetCDF'
    CHARACTER(len=*), INTENT(in), OPTIONAL :: post_suf     ! suffix of output file
    CHARACTER(len=*), INTENT(in), OPTIONAL :: rest_suf     ! suffix of restart file
    CHARACTER(len=*), INTENT(in), OPTIONAL :: init_suf     ! suffix of initial file
    LOGICAL,          INTENT(in), OPTIONAL :: lpost        ! write to  output file
    LOGICAL,          INTENT(in), OPTIONAL :: lrestart     ! write to restart file
    LOGICAL,          INTENT(in), OPTIONAL :: linitial     ! read from initial file
    !
    INTEGER :: i
    !
    CALL message('new_var_list','adding new var_list '//TRIM(name))
    !
    ! name must be unique
    !
    NULLIFY (this_list)
    !
    ! look if name exists already in list
    !
    IF (ANY(var_lists(1:nvar_lists)%name == name)) THEN
      CALL finish('new_list', 'output var_list '//TRIM(name)//' already used.')
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
        CALL finish('new_list', 'var_lists container overflow, increase "max_var_lists" in mo_var_list.f90')
      ENDIF
      this_list => var_lists(nvar_lists)
    ENDIF
    !
    CALL new_list (this_list)
    !
    ! set default list characteristics
    !
    this_list%name     = name
    this_list%post_suf = '_'//TRIM(name)
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
    TYPE(t_var_list), POINTER    :: this_list ! pointer
    CHARACTER(len=*), INTENT(in) :: name      ! name of output var_list
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
       &                   post_suf, rest_suf, init_suf, lpost,  &
       &                   lrestart, linitial)
    !
    TYPE(t_var_list), INTENT(inout)        :: this_list   ! output var_list to change
    INTEGER,          INTENT(in), OPTIONAL :: output_type    ! 'GRIB' or 'NetCDF'
    INTEGER,          INTENT(in), OPTIONAL :: restart_type   ! 'GRIB' or 'NetCDF'
    CHARACTER(len=*), INTENT(in), OPTIONAL :: post_suf       ! suffix of output  file
    CHARACTER(len=*), INTENT(in), OPTIONAL :: rest_suf       ! suffix of restart file
    CHARACTER(len=*), INTENT(in), OPTIONAL :: init_suf       ! suffix of initial file
    LOGICAL,          INTENT(in), OPTIONAL :: lpost          ! in standard output file
    LOGICAL,          INTENT(in), OPTIONAL :: lrestart       ! in standard restartfile
    LOGICAL,          INTENT(in), OPTIONAL :: linitial       ! in standard initialfile
    !
    IF (PRESENT(output_type)) CALL assign_if_present(this_list%output_type,  output_type)
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
    TYPE(t_var_list),     INTENT(in)  :: this_list    ! list
    CHARACTER(len=*),     INTENT(in)  :: name         ! name of variable
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
  SUBROUTINE default_var_list_settings (this_list,                                  &
       &                                filename,                                   &
       &                                lpost, lrestart, linitial,                  &
       &                                post_suf, rest_suf, init_suf,               &
       &                                output_type, restart_type, compression_type) 
    !
    TYPE(t_var_list),   INTENT(inout)        :: this_list        ! output var_list
    LOGICAL,            INTENT(in), OPTIONAL :: lpost            ! to output 
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart         ! from/to restart
    LOGICAL,            INTENT(in), OPTIONAL :: linitial         ! from initial
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: filename         ! name of output file
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: post_suf         ! suffix of output  file         
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: rest_suf         ! suffix of restart file
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: init_suf         ! suffix of initial file
    INTEGER,            INTENT(in), OPTIONAL :: output_type      ! output file type
    INTEGER,            INTENT(in), OPTIONAL :: restart_type     ! restart file type 
    INTEGER,            INTENT(in), OPTIONAL :: compression_type ! compression type
    !
    CALL assign_if_present (this_list%lpost,            lpost)
    CALL assign_if_present (this_list%lrestart,         lrestart)
    CALL assign_if_present (this_list%linitial,         linitial)
    CALL assign_if_present (this_list%filename,         filename)
    CALL assign_if_present (this_list%post_suf,         post_suf)
    CALL assign_if_present (this_list%rest_suf,         rest_suf)
    CALL assign_if_present (this_list%init_suf,         init_suf)
    CALL assign_if_present (this_list%output_type,      output_type)
    CALL assign_if_present (this_list%restart_type,     restart_type)
    CALL assign_if_present (this_list%compression_type, compression_type)
    !
  END SUBROUTINE default_var_list_settings
  !------------------------------------------------------------------------------------------------
  FUNCTION default_var_list_metadata(this_list) RESULT(this_info)
    !
    TYPE(t_var_metadata)         :: this_info        ! memory info structure
    !
    TYPE(t_var_list), INTENT(in) :: this_list        ! output var_list    
    !
    this_info%key                 = 0    
    this_info%name                = ''
    !
    this_info%cf                  = t_cf_var('', '', '')
    this_info%grib2               = t_grib2_var(-1, -1, -1, -1, -1, -1, -1)
    !
    this_info%allocated           = .FALSE.
    this_info%ndims               = 0
    this_info%used_dimensions(:)  = 0
    !
    this_info%lpost               = this_list%lpost
    this_info%laccu               = this_list%laccu
    this_info%resetval            = t_union_vals( 0.0_dp, 0, .FALSE.)
    this_info%lrestart            = this_list%lrestart
    this_info%lrestart_cont       = .FALSE.
    this_info%lrestart_read       = .FALSE.
    this_info%lmiss               = this_list%lmiss
    this_info%missval             = t_union_vals( 0.0_dp, 0, .FALSE.)
    this_info%initval             = t_union_vals( 0.0_dp, 0, .FALSE.)
    !
    this_info%hgrid               = -1
    this_info%vgrid               = -1
    !
  END FUNCTION default_var_list_metadata
  !------------------------------------------------------------------------------------------------
  !
  ! Set parameters of list element already created
  ! (private routine within this module)
  !
  ! Set each parameter in data type var_metadata if the respective 
  ! optional parameter is present.
  !
  SUBROUTINE set_var_metadata (info,                                    &
         &                     name, hgrid, vgrid, cf, grib2, ldims,    &
         &                     lpost, lrestart, lrestart_cont, initval, &
         &                     laccu, resetval, lmiss, missval, verbose)
    !
    TYPE(t_var_metadata), INTENT(inout)        :: info          ! memory info struct.
    CHARACTER(len=*),     INTENT(in), OPTIONAL :: name          ! variable name
    INTEGER,              INTENT(in), OPTIONAL :: hgrid         ! horizontal grid type used 
    INTEGER,              INTENT(in), OPTIONAL :: vgrid         ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in), OPTIONAL :: cf            ! CF convention
    TYPE(t_grib2_var),    INTENT(in), OPTIONAL :: grib2         ! GRIB2
    INTEGER,              INTENT(in), OPTIONAL :: ldims(:)      ! used dimensions 
    LOGICAL,              INTENT(in), OPTIONAL :: lpost         ! into output var_list
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart      ! restart file flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont ! continue on restart
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: initval       ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu         ! accumulation flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: resetval      ! reset value
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss         ! missing value flag
    TYPE(t_union_vals),   INTENT(in), OPTIONAL :: missval       ! missing value
    LOGICAL,              INTENT(in), OPTIONAL :: verbose
    !
    LOGICAL :: lverbose
    !
    ! set flags from optional parameters
    !
    lverbose = .FALSE. 
    CALL assign_if_present (lverbose, verbose)
    !
    ! set components describing the 'Content of the field'
    !
    CALL assign_if_present (info%name,  name)
    CALL assign_if_present (info%cf,    cf)
    CALL assign_if_present (info%grib2, grib2)
    !
    CALL assign_if_present (info%used_dimensions(1:SIZE(ldims)), ldims)
    !
    ! set grid type
    !
    CALL assign_if_present (info%hgrid, hgrid)
    CALL assign_if_present (info%vgrid, vgrid)
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
  SUBROUTINE add_var_list_element_r4d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, lpost,           &
       lrestart, lrestart_cont, initval_r, laccu, resetval_r,  &
       lmiss, missval_r, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(dp),             POINTER              :: ptr(:,:,:,:)        ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(4)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(dp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    REAL(dp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(dp),             POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    new_list_element%field%info = default_var_list_metadata(this_list)
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
    CALL assign_if_present(missval%rval, missval_r)
    CALL assign_if_present(initval%rval, initval_r)
    CALL assign_if_present(resetval%rval, resetval_r)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:4), ldims(1:4))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = new_list_element%field%info%used_dimensions(4)
      new_list_element%field%info%ndims = 4
      new_list_element%field%var_base_size = 8
      ALLOCATE(new_list_element%field%r_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r4d',      &
                   &'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
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
  SUBROUTINE add_var_list_element_r3d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, lpost,           &
       lrestart, lrestart_cont, initval_r, laccu, resetval_r,  &
       lmiss, missval_r, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(dp),             POINTER              :: ptr(:,:,:)          ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(dp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    REAL(dp),             INTENT(in), OPTIONAL :: resetval_r            ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r             ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(dp),             POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    new_list_element%field%info = default_var_list_metadata(this_list)
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
    CALL assign_if_present(missval%rval, missval_r)
    CALL assign_if_present(initval%rval, initval_r)
    CALL assign_if_present(resetval%rval, resetval_r)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:3), ldims(1:3))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = 1
      new_list_element%field%info%ndims = 3
      new_list_element%field%var_base_size = 8
      ALLOCATE (new_list_element%field%r_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r3d',      &
                   &'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
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
  SUBROUTINE add_var_list_element_r2d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, lpost,           &
       lrestart, lrestart_cont, initval_r, laccu, resetval_r,  &
       lmiss, missval_r, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(dp),             POINTER              :: ptr(:,:)            ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(2)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(dp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    REAL(dp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(dp),             POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    new_list_element%field%info = default_var_list_metadata(this_list)
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
    CALL assign_if_present(missval%rval, missval_r)
    CALL assign_if_present(initval%rval, initval_r)
    CALL assign_if_present(resetval%rval, resetval_r)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:2), ldims(1:2))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = 1
      idims(4) = 1
      new_list_element%field%info%ndims = 2
      new_list_element%field%var_base_size = 8
      ALLOCATE (new_list_element%field%r_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r2d',      &
                   &'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
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
  SUBROUTINE add_var_list_element_r1d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, lpost,           &
       lrestart, lrestart_cont, initval_r, laccu, resetval_r,  &
       lmiss, missval_r, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(dp),             POINTER              :: ptr(:)              ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(1)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(dp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    REAL(dp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(dp),             POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    new_list_element%field%info = default_var_list_metadata(this_list)
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
    CALL assign_if_present(missval%rval, missval_r)
    CALL assign_if_present(initval%rval, initval_r)
    CALL assign_if_present(resetval%rval, resetval_r)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:1), ldims(1:1))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = 1
      idims(3) = 1
      idims(4) = 1
      new_list_element%field%info%ndims = 1
      new_list_element%field%var_base_size = 8
      ALLOCATE (new_list_element%field%r_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r1d',      &
                   &'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
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
  SUBROUTINE add_var_list_element_i4d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, lpost,           &
       lrestart, lrestart_cont, initval_i, laccu, resetval_i,  &
       lmiss, missval_i, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    INTEGER,              POINTER              :: ptr(:,:,:,:)        ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(4)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    INTEGER,              INTENT(in), OPTIONAL :: initval_i           ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    INTEGER,              INTENT(in), OPTIONAL :: resetval_i          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    INTEGER,              INTENT(in), OPTIONAL :: missval_i           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    INTEGER,              POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    new_list_element%field%info = default_var_list_metadata(this_list)
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
    CALL assign_if_present(missval%ival, missval_i)
    CALL assign_if_present(initval%ival, initval_i)
    CALL assign_if_present(resetval%ival, resetval_i)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:4), ldims(1:4))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = new_list_element%field%info%used_dimensions(4)
      new_list_element%field%info%ndims = 4
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%i_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_i4d',      &
                   &'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
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
  SUBROUTINE add_var_list_element_i3d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, lpost,           &
       lrestart, lrestart_cont, initval_i, laccu, resetval_i,  &
       lmiss, missval_i, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    INTEGER,              POINTER              :: ptr(:,:,:)          ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    INTEGER,              INTENT(in), OPTIONAL :: initval_i           ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    INTEGER,              INTENT(in), OPTIONAL :: resetval_i          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    INTEGER,              INTENT(in), OPTIONAL :: missval_i           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    INTEGER,              POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    new_list_element%field%info = default_var_list_metadata(this_list)
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
    CALL assign_if_present(missval%ival, missval_i)
    CALL assign_if_present(initval%ival, initval_i)
    CALL assign_if_present(resetval%ival, resetval_i)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:3), ldims(1:3))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = 1
      new_list_element%field%info%ndims = 3
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%i_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r3d',      &
                   &'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
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
  SUBROUTINE add_var_list_element_i2d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, lpost,           &
       lrestart, lrestart_cont, initval_i, laccu, resetval_i,  &
       lmiss, missval_i, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    INTEGER,              POINTER              :: ptr(:,:)            ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(2)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    INTEGER,              INTENT(in), OPTIONAL :: initval_i           ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    INTEGER,              INTENT(in), OPTIONAL :: resetval_i          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    INTEGER,              INTENT(in), OPTIONAL :: missval_i           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    INTEGER,              POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    new_list_element%field%info = default_var_list_metadata(this_list)
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
    CALL assign_if_present(missval%ival, missval_i)
    CALL assign_if_present(initval%ival, initval_i)
    CALL assign_if_present(resetval%ival, resetval_i)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:2), ldims(1:2))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = 1
      idims(4) = 1
      new_list_element%field%info%ndims = 2
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%i_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r2d',      &
                   &'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
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
  SUBROUTINE add_var_list_element_i1d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, lpost,           &
       lrestart, lrestart_cont, initval_i, laccu, resetval_i,  &
       lmiss, missval_i, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    INTEGER,              POINTER              :: ptr(:)              ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(1)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    INTEGER,              INTENT(in), OPTIONAL :: initval_i           ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    INTEGER,              INTENT(in), OPTIONAL :: resetval_i          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    INTEGER,              INTENT(in), OPTIONAL :: missval_i           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    INTEGER,              POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    new_list_element%field%info = default_var_list_metadata(this_list)
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
    CALL assign_if_present(missval%ival, missval_i)
    CALL assign_if_present(initval%ival, initval_i)
    CALL assign_if_present(resetval%ival, resetval_i)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:1), ldims(1:1))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = 1
      idims(3) = 1
      idims(4) = 1
      new_list_element%field%info%ndims = 1
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%i_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r1d',      &
                   &'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
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
  SUBROUTINE add_var_list_element_l4d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, lpost,           &
       lrestart, lrestart_cont, initval_l, laccu, resetval_l,  &
       lmiss, missval_l, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    LOGICAL,              POINTER              :: ptr(:,:,:,:)        ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(4)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: initval_l           ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    LOGICAL,              INTENT(in), OPTIONAL :: resetval_l          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    LOGICAL,              INTENT(in), OPTIONAL :: missval_l           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    LOGICAL,              POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    new_list_element%field%info = default_var_list_metadata(this_list)
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
    CALL assign_if_present(missval%lval, missval_l)
    CALL assign_if_present(initval%lval, initval_l)
    CALL assign_if_present(resetval%lval, resetval_l)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:4), ldims(1:4))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = new_list_element%field%info%used_dimensions(4)
      new_list_element%field%info%ndims = 4
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%l_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r4d',      &
                   &'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
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
  SUBROUTINE add_var_list_element_l3d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, lpost,           &
       lrestart, lrestart_cont, initval_l, laccu, resetval_l,  &
       lmiss, missval_l, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    LOGICAL,              POINTER              :: ptr(:,:,:)          ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: initval_l           ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    LOGICAL,              INTENT(in), OPTIONAL :: resetval_l          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    LOGICAL,              INTENT(in), OPTIONAL :: missval_l           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    LOGICAL,              POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    new_list_element%field%info = default_var_list_metadata(this_list)
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
    CALL assign_if_present(missval%lval, missval_l)
    CALL assign_if_present(initval%lval, initval_l)
    CALL assign_if_present(resetval%lval, resetval_l)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:3), ldims(1:3))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = 1
      new_list_element%field%info%ndims = 3
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%l_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r3d',      &
                   &'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
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
       hgrid, vgrid, cf, grib2, ldims, lpost,       &
       lrestart, lrestart_cont, initval_l, laccu, resetval_l,  &
       lmiss, missval_l, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    LOGICAL,              POINTER              :: ptr(:,:)            ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(2)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: initval_l           ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    LOGICAL,              INTENT(in), OPTIONAL :: resetval_l          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    LOGICAL,              INTENT(in), OPTIONAL :: missval_l           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    LOGICAL,              POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    new_list_element%field%info = default_var_list_metadata(this_list)
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
    CALL assign_if_present(missval%lval, missval_l)
    CALL assign_if_present(initval%lval, initval_l)
    CALL assign_if_present(resetval%lval, resetval_l)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:2), ldims(1:2))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = 1
      idims(4) = 1
      new_list_element%field%info%ndims = 2
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%l_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r2d',      &
                   &'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
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
  SUBROUTINE add_var_list_element_l1d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, lpost,           &
       lrestart, lrestart_cont, initval_l, laccu, resetval_l,  &
       lmiss, missval_l, info, p4, verbose)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    LOGICAL,              POINTER              :: ptr(:)              ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(1)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: lpost               ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: initval_l           ! value if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: laccu               ! accumulation flag
    LOGICAL,              INTENT(in), OPTIONAL :: resetval_l          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    LOGICAL,              INTENT(in), OPTIONAL :: missval_l           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    LOGICAL,              POINTER,    OPTIONAL :: p4(:,:,:,:)         ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(4)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    new_list_element%field%info = default_var_list_metadata(this_list)
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
    CALL assign_if_present(missval%lval, missval_l)
    CALL assign_if_present(initval%lval, initval_l)
    CALL assign_if_present(resetval%lval, resetval_l)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, lpost=lpost,                    &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         laccu=laccu, resetval=resetval, lmiss=lmiss, missval=missval,    & 
         verbose=verbose)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:1), ldims(1:1))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = 1
      idims(3) = 1
      idims(4) = 1
      new_list_element%field%info%ndims = 1
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%l_ptr(idims(1), idims(2), idims(3), idims(4)), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r1d',      &
                   &'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
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
  SUBROUTINE delete_var_list_element (this_list, name)
    TYPE(t_var_list), INTENT(inout) :: this_list
    CHARACTER(len=*), INTENT(in)    :: name   
    !
    TYPE(t_list_element), POINTER :: ptr
    !
    IF (this_list%first_list_element%field%info%name == name) THEN
      CALL delete_list_element (this_list, this_list%first_list_element)
      RETURN
    ELSE
      ptr => this_list%first_list_element
      DO
        IF (.NOT.ASSOCIATED (ptr%next_list_element)) EXIT
        IF (ptr%next_list_element%field%info%name == name) THEN
          CALL delete_list_element (this_list, ptr%next_list_element)
          EXIT
        ENDIF
        ptr => ptr%next_list_element
      END DO
    ENDIF
    !
  END SUBROUTINE delete_var_list_element
  !------------------------------------------------------------------------------------------------
  !
  ! Get element of a list
  !
  ! Specific routines for pointers of different rank
  !
  !================================================================================================ 
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  ! obtain pointer to 4d-field
  !
  SUBROUTINE get_var_list_element_r4d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    REAL(dp),         POINTER    :: ptr(:,:,:,:) ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%r_ptr
    !
  END SUBROUTINE get_var_list_element_r4d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 3d-field
  !
  SUBROUTINE get_var_list_element_r3d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    REAL(dp),         POINTER    :: ptr(:,:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,:,:,1)
    !
  END SUBROUTINE get_var_list_element_r3d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 2d-field
  !
  SUBROUTINE get_var_list_element_r2d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list  ! list
    CHARACTER(len=*), INTENT(in) :: name       ! name of variable
    REAL(dp),         POINTER    :: ptr(:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,:,1,1)
    !
  END SUBROUTINE get_var_list_element_r2d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 1d-field
  !
  SUBROUTINE get_var_list_element_r1d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list ! list
    CHARACTER(len=*), INTENT(in) :: name      ! name of variable
    REAL(dp),         POINTER    :: ptr(:)    ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,1,1,1)
    !
  END SUBROUTINE get_var_list_element_r1d
  !================================================================================================ 
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  ! obtain pointer to 4d-field
  !
  SUBROUTINE get_var_list_element_i4d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    INTEGER,          POINTER    :: ptr(:,:,:,:) ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%i_ptr
    !
  END SUBROUTINE get_var_list_element_i4d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 3d-field
  !
  SUBROUTINE get_var_list_element_i3d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    INTEGER,          POINTER    :: ptr(:,:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%i_ptr(:,:,:,1)
    !
  END SUBROUTINE get_var_list_element_i3d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 2d-field
  !
  SUBROUTINE get_var_list_element_i2d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list  ! list
    CHARACTER(len=*), INTENT(in) :: name       ! name of variable
    INTEGER,          POINTER    :: ptr(:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%i_ptr(:,:,1,1)
    !
  END SUBROUTINE get_var_list_element_i2d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 1d-field
  !
  SUBROUTINE get_var_list_element_i1d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list ! list
    CHARACTER(len=*), INTENT(in) :: name      ! name of variable
    INTEGER,          POINTER    :: ptr(:)    ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%i_ptr(:,1,1,1)
    !
  END SUBROUTINE get_var_list_element_i1d
  !================================================================================================ 
  ! LOGICAL SECTION -------------------------------------------------------------------------------
  !
  ! obtain pointer to 4d-field
  !
  SUBROUTINE get_var_list_element_l4d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    LOGICAL,          POINTER    :: ptr(:,:,:,:) ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%l_ptr
    !
  END SUBROUTINE get_var_list_element_l4d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 3d-field
  !
  SUBROUTINE get_var_list_element_l3d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    LOGICAL,          POINTER    :: ptr(:,:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%l_ptr(:,:,:,1)
    !
  END SUBROUTINE get_var_list_element_l3d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 2d-field
  !
  SUBROUTINE get_var_list_element_l2d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list  ! list
    CHARACTER(len=*), INTENT(in) :: name       ! name of variable
    LOGICAL,          POINTER    :: ptr(:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%l_ptr(:,:,1,1)
    !
  END SUBROUTINE get_var_list_element_l2d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 1d-field
  !
  SUBROUTINE get_var_list_element_l1d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list ! list
    CHARACTER(len=*), INTENT(in) :: name      ! name of variable
    LOGICAL,          POINTER    :: ptr(:)    ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%l_ptr(:,1,1,1)
    !
  END SUBROUTINE get_var_list_element_l1d
  !------------------------------------------------------------------------------------------------
  !
  ! Print routines for control output and debuggung
  !
  SUBROUTINE print_memory_use (this_list, ldetailed)
    TYPE(t_var_list) ,INTENT(in) :: this_list ! list
    LOGICAL, INTENT(in), OPTIONAL :: ldetailed
    !
    IF (PRESENT(ldetailed)) THEN
      WRITE (message_text,'(a32,a,a,i10,a,i4,a)')                  &
           TRIM(this_list%name), '-buffer: ',                      &
           'Memory in use: ', this_list%memory_used, ' bytes in ', &
           this_list%list_elements, ' fields.'
    ELSE
      WRITE (message_text,'(a32,a,a,i10,a,i6,a)')                    &
           TRIM(this_list%name), '-buffer: ',                        &
           'Memory in use: ', this_list%memory_used/1024, ' kb in ', &
           this_list%list_elements, ' fields.'
    ENDIF
    CALL message('',message_text)
    !
  END SUBROUTINE print_memory_use
  !------------------------------------------------------------------------------------------------
  !
  ! print current memory table 
  !
  SUBROUTINE print_var_list (this_list)
    TYPE(t_var_list),  INTENT(in) :: this_list ! list
    !
    TYPE(t_list_element), POINTER :: this_list_element
    CHARACTER(len=32) :: dimension_text, dtext
    INTEGER :: i
    !
    CALL message('','')
    CALL message('','')
    CALL message('','Status of variable list '//TRIM(this_list%name)//':')    
    CALL message('','')
    !
    this_list_element => this_list%first_list_element
    !
    DO WHILE (ASSOCIATED(this_list_element))
      ! 
      IF (this_list_element%field%info%name /= '') THEN
        !
        WRITE (message_text,'(a,a)')       &
             'Table entry name                            : ', &
             TRIM(this_list_element%field%info%name)
        CALL message('', message_text)
        !
        IF (ASSOCIATED(this_list_element%field%r_ptr) .OR. &
          & ASSOCIATED(this_list_element%field%i_ptr) .OR. &
          & ASSOCIATED(this_list_element%field%l_ptr)) THEN
          CALL message ('','Pointer status                              : in use.')
          dimension_text = '('
          DO i = 1, this_list_element%field%info%ndims
            WRITE(dtext,'(i0)') this_list_element%field%info%used_dimensions(i)
            IF (this_list_element%field%info%ndims == i) THEN
              dimension_text = TRIM(dimension_text)//TRIM(dtext)//')'
            ELSE
              dimension_text = TRIM(dimension_text)//TRIM(dtext)//','
            ENDIF
          ENDDO
          WRITE (message_text,'(a,a)') &
               'Local field dimensions                      : ', TRIM(dimension_text)
          CALL message('', message_text)
        ELSE
          CALL message('', 'Pointer status                              : not in use.')
        ENDIF
        !        
        WRITE (message_text,'(a,3i3)') &
             'Assigned GRIB discipline/category/parameter : ', &
             this_list_element%field%info%grib2%discipline,    &
             this_list_element%field%info%grib2%category,      &
             this_list_element%field%info%grib2%parameter
        CALL message('', message_text)
        !
        WRITE (message_text,'(a,a,a,a)') &
             ' CF convention standard name/unit            : ',    &
             TRIM(this_list_element%field%info%cf%standard_name), &
             '     ',                                             &
             TRIM(this_list_element%field%info%cf%units)
        CALL message('', message_text)
        !
        IF (this_list_element%field%info%laccu) THEN
          CALL message('', 'Accumulation                                : on.')
        ELSE
          CALL message('', 'Accumulation                                : off.')
        ENDIF
        !
        IF (this_list_element%field%info%lmiss) THEN
          IF (ASSOCIATED(this_list_element%field%r_ptr)) THEN
            WRITE (message_text,'(a,e20.12)')      &
                 'Missing value                               : ', &
                 this_list_element%field%info%missval%rval
          ELSE IF (ASSOCIATED(this_list_element%field%i_ptr)) THEN
            WRITE (message_text,'(a,i8)')      &
                 'Missing value                               : ', &
                 this_list_element%field%info%missval%ival
          ELSE IF (ASSOCIATED(this_list_element%field%l_ptr)) THEN
            WRITE (message_text,'(a,l8)')      &
                 'Missing value                               : ', &
                 this_list_element%field%info%missval%lval
          ENDIF
          CALL message('', message_text)
        ELSE
          CALL message('', 'Missing values                              : off.')
        ENDIF
        !
        IF (this_list_element%field%info%lrestart) THEN
          CALL message('', 'Added to restart                            : yes.')
        ELSE
          CALL message('', 'Added to Restart                            : no.')
        ENDIF
        CALL message('', '')
      ENDIF
      !
      ! select next element in linked list 
      !
      this_list_element => this_list_element%next_list_element
    ENDDO

    !    
  END SUBROUTINE print_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! print current stat table 
  !
  SUBROUTINE print_sinfo (this_list)
    TYPE(t_var_list),  INTENT(in) :: this_list
    !
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
  SUBROUTINE add_var_list_reference (to_var_list, name, from_var_list, lpost, bit_precision)
    TYPE(t_var_list), INTENT(inout)          :: to_var_list
    CHARACTER(len=*), INTENT(in)             :: name
    CHARACTER(len=*), INTENT(in)             :: from_var_list
    LOGICAL,          INTENT(in),   OPTIONAL :: lpost
    INTEGER,          INTENT(in),   OPTIONAL :: bit_precision
    !
    TYPE(t_var_list_element), POINTER :: source
    TYPE(t_list_element),     POINTER :: new_list_element
    !
    CALL locate (source, name, from_var_list)
    IF (ASSOCIATED(source)) THEN
      CALL append_list_element (to_var_list, new_list_element)
      new_list_element%field                = source
      new_list_element%field%info%allocated = .FALSE.
      new_list_element%field%info%lrestart  = .FALSE.
      CALL assign_if_present(new_list_element%field%info%lpost, lpost)
      CALL assign_if_present(new_list_element%field%info%grib2%bits, bit_precision)
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
  SUBROUTINE assign_if_present_union (y,x)
    TYPE(t_union_vals), INTENT(inout)        :: y
    TYPE(t_union_vals) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    y = x
  END SUBROUTINE assign_if_present_union
  !------------------------------------------------------------------------------------------------
END MODULE mo_var_list
