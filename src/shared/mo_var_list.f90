MODULE mo_var_list

  USE mo_kind,             ONLY: dp
  USE mo_cf_convention,    ONLY: cf_var
  USE mo_grib1,            ONLY: grib1_var
  USE mo_grib2,            ONLY: grib2_var
  USE mo_var_metadata,     ONLY: var_metadata
  USE mo_var_list_element, ONLY: var_list_element
  USE mo_linked_list,      ONLY: var_list, new_list, list_element, find_list_element, &
                                 append_list_element 
  USE mo_exception,        ONLY: message, message_text, finish

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: new_var_list             ! get a pointer to a new output var_list
  PUBLIC :: delete_var_list          ! delete an output var_list
  PUBLIC :: delete_var_lists         ! delete all output var_lists
  PUBLIC :: get_var_list             ! get a pointer to an existing output var_list
  PUBLIC :: set_var_list             ! set default parameters of an output var_list

  PUBLIC :: var_lists                ! vector of output var_lists
  PUBLIC :: nvar_lists               ! number of output var_lists defined so far

  PUBLIC :: add_var                  ! create/allocate a new var_list list entry
  PUBLIC :: get_var                  ! obtain reference to existing list entry

  INTERFACE add_var
     MODULE PROCEDURE add_var_list_element_4d ! create a new list entry
     MODULE PROCEDURE add_var_list_element_3d 
     MODULE PROCEDURE add_var_list_element_2d 
     MODULE PROCEDURE add_var_list_element_1d 
  END INTERFACE
  
  INTERFACE get_var
     MODULE PROCEDURE get_var_list_element_4d ! obtain reference to a list entry
     MODULE PROCEDURE get_var_list_element_3d
     MODULE PROCEDURE get_var_list_element_2d
     MODULE PROCEDURE get_var_list_element_1d
  END INTERFACE

  INTERFACE assign_if_present
     MODULE PROCEDURE assign_if_present_character
     MODULE PROCEDURE assign_if_present_logical
     MODULE PROCEDURE assign_if_present_integer
     MODULE PROCEDURE assign_if_present_integers
     MODULE PROCEDURE assign_if_present_real
     MODULE PROCEDURE assign_if_present_cf
     MODULE PROCEDURE assign_if_present_grib1
     MODULE PROCEDURE assign_if_present_grib2
  END INTERFACE

  INTEGER, PARAMETER           :: max_var_lists  = 64    ! max number of output var_lists
  LOGICAL,                SAVE :: lvar_list_init = .FALSE.
  INTEGER,                SAVE :: nvar_lists =  0    ! var_lists allocated so far
  TYPE(var_list), TARGET, SAVE :: var_lists(max_var_lists) ! memory buffer array

CONTAINS
  !------------------------------------------------------------------------------------------------
  !
  ! Create a new memory buffer / output var_list
  ! Get a pointer to the new var_list
  !
  SUBROUTINE new_var_list (ptr_var_list, name, output_type, restart_type, &
       post_suf, rest_suf, init_suf, lpost, lrestart, linitial)
    TYPE(var_list),   POINTER, INTENT(out)          :: ptr_var_list ! anchor
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
    NULLIFY (ptr_var_list)
    !
    IF (lvar_list_init) THEN
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
          ptr_var_list => var_lists(i)
          EXIT
        ENDIF
      END DO
      IF(.NOT. ASSOCIATED(ptr_var_list)) THEN
        nvar_lists = nvar_lists + 1
        IF (nvar_lists > max_var_lists) THEN
          CALL finish('new_list', 'var_lists container overflow, increase "max_var_lists"')
        ENDIF
        ptr_var_list => var_lists(nvar_lists)
      ENDIF
      CALL new_list (ptr_var_list)
      !
      ! set default list characteristics
      !
      ptr_var_list%name     = name
      ptr_var_list%post_suf = '_'//name
      ptr_var_list%rest_suf = ptr_var_list%post_suf
      ptr_var_list%init_suf = ptr_var_list%post_suf
      !
      ! set non-default list characteristics
      !
      CALL assign_if_present(ptr_var_list%output_type,  output_type)
      CALL assign_if_present(ptr_var_list%restart_type, restart_type)
      CALL assign_if_present(ptr_var_list%post_suf,     post_suf) 
      CALL assign_if_present(ptr_var_list%rest_suf,     rest_suf) 
      CALL assign_if_present(ptr_var_list%init_suf,     init_suf) 
      CALL assign_if_present(ptr_var_list%lpost,        lpost)
      CALL assign_if_present(ptr_var_list%lrestart,     lrestart)
      CALL assign_if_present(ptr_var_list%linitial,     linitial)
      !
    ELSE
      ! definition of var_list structure is available
      CALL get_var_list(ptr_var_list, name)
      !
    END IF
    !
  END SUBROUTINE new_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! Get a reference to a memory buffer / output var_list
  !
  SUBROUTINE get_var_list (ptr_var_list, name)
    TYPE(var_list),  POINTER, INTENT(out) :: ptr_var_list ! pointer
    CHARACTER(len=*),          INTENT(in)  :: name         ! name of output var_list
    !
    INTEGER :: i
    !
    NULLIFY (ptr_var_list)
    !
    DO i = 1, nvar_lists
      IF (var_lists(i)%name == name) THEN
        ptr_var_list => var_lists(i)
        EXIT
      ENDIF
    END DO
    !
  END SUBROUTINE get_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! Change parameters of an already existent output var_list
  !
  SUBROUTINE set_var_list (ptr_var_list, output_type, restart_type, &
       post_suf, rest_suf, init_suf, lpost, lrestart, linitial)
    TYPE(var_list),   INTENT(inout)        :: ptr_var_list   ! output var_list to change
    INTEGER,          INTENT(in), OPTIONAL :: output_type    ! 'GRIB' or 'NetCDF'
    INTEGER,          INTENT(in), OPTIONAL :: restart_type   ! 'GRIB' or 'NetCDF'
    CHARACTER(len=*), INTENT(in), OPTIONAL :: post_suf       ! suffix of output  file
    CHARACTER(len=*), INTENT(in), OPTIONAL :: rest_suf       ! suffix of restart file
    CHARACTER(len=*), INTENT(in), OPTIONAL :: init_suf       ! suffix of initial file
    LOGICAL,          INTENT(in), OPTIONAL :: lpost          ! in standard output file
    LOGICAL,          INTENT(in), OPTIONAL :: lrestart       ! in standard restartfile
    LOGICAL,          INTENT(in), OPTIONAL :: linitial       ! in standard initialfile
    !
    CALL assign_if_present(ptr_var_list%output_type,  output_type)
    CALL assign_if_present(ptr_var_list%restart_type, restart_type)
    CALL assign_if_present(ptr_var_list%post_suf,     post_suf) 
    CALL assign_if_present(ptr_var_list%rest_suf,     rest_suf) 
    CALL assign_if_present(ptr_var_list%init_suf,     init_suf) 
    CALL assign_if_present(ptr_var_list%lpost,        lpost)
    CALL assign_if_present(ptr_var_list%lrestart,     lrestart)
    CALL assign_if_present(ptr_var_list%linitial,     linitial)
    !
  END SUBROUTINE set_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! Delete an output var_list, nullify the associated pointer
  !
  SUBROUTINE delete_var_list(ptr_var_list)
    TYPE(var_list), POINTER :: ptr_var_list
    !
    IF (ASSOCIATED(ptr_var_list)) THEN
      CALL destruct_list(ptr_var_list)
      NULLIFY(ptr_var_list)
    ENDIF
    !
  END SUBROUTINE delete_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! Delete all output var_lists 
  !
  SUBROUTINE delete_var_lists
    TYPE(var_list), POINTER :: ptr_var_list
    !
    INTEGER :: i
    !
    DO i = 1, nvar_lists
      ptr_var_list => var_lists(i)
      CALL delete_var_list (ptr_var_list)
    END DO
    !
    lvar_list_init = .FALSE.
    !
  END SUBROUTINE delete_var_lists
  !------------------------------------------------------------------------------------------------
  !
  ! Get a copy of the metadata concerning a var_list element
  !
  SUBROUTINE get_var_list_element_info (ptr_var_list, name, info)
    TYPE(var_list),     INTENT(in)  :: ptr_var_list ! list
    CHARACTER(len=*),   INTENT(in)  :: name         ! name of variable
    TYPE(var_metadata), INTENT(out) :: info         ! variable meta data
    !
    TYPE(list_element), POINTER :: ptr
    ! 
    ptr => find_list_element (ptr_var_list, name)
    IF (ASSOCIATED (ptr)) THEN
      info = ptr%field%info
    ENDIF
    !
  END SUBROUTINE get_var_list_element_info
  !------------------------------------------------------------------------------------------------
  !
  ! Set default meta data of output var_list
  !
  SUBROUTINE default_var_list_setting (ptr_var_list, &
       name, cf, grib1, grib2, lrestart, lpost, laccu, reset, lmiss, missval, &
       lrestart_cont, grid_representation)
    TYPE(var_list)  ,INTENT(inout)         :: ptr_var_list        ! output var_list
    CHARACTER(len=*), INTENT(in), OPTIONAL :: name          ! variable name
    TYPE(cf_var),     INTENT(in), OPTIONAL :: cf         
    TYPE(grib1_var),  INTENT(in), OPTIONAL :: grib1
    TYPE(grib2_var),  INTENT(in), OPTIONAL :: grib2
    LOGICAL,          INTENT(in), OPTIONAL :: lrestart      ! restart file flag
    LOGICAL,          INTENT(in), OPTIONAL :: lpost         ! into output var_list
    LOGICAL,          INTENT(in), OPTIONAL :: laccu         ! accumulation flag
    REAL(dp),         INTENT(in), OPTIONAL :: reset         ! reset value
    LOGICAL,          INTENT(in), OPTIONAL :: lmiss         ! missing value flag
    REAL(dp),         INTENT(in), OPTIONAL :: missval       ! missing value
    LOGICAL,          INTENT(in), OPTIONAL :: lrestart_cont ! continue on restart
    INTEGER,          INTENT(in), OPTIONAL :: grid_representation

! more complicated ... haven't got default_info anymore ...
!    CALL set_var_list_metadata (ptr_var_list%default_info                                &
!       name, cf, grib1, grib2, lrestart, lpost, laccu, reset, lmiss, missval, &
!       lrestart_cont, grid_representation, verbose=.FALSE.)

  END SUBROUTINE default_var_list_setting
  !------------------------------------------------------------------------------------------------
  !
  ! Change metadata of a var_list element
  !
  SUBROUTINE set_var_list_element_metadata (ptr_var_list, &
       name, cf, grib1, grib2, lrestart, lpost, laccu, reset, lmiss, missval, &
       lrestart_cont, grid_representation, verbose)
    TYPE(var_list),   INTENT(inout)        :: ptr_var_list  ! memory info struct.
    CHARACTER(len=*), INTENT(in), OPTIONAL :: name          ! variable name
    TYPE(cf_var),     INTENT(in), OPTIONAL :: cf         
    TYPE(grib1_var),  INTENT(in), OPTIONAL :: grib1
    TYPE(grib2_var),  INTENT(in), OPTIONAL :: grib2
    LOGICAL,          INTENT(in), OPTIONAL :: lrestart      ! restart file flag
    LOGICAL,          INTENT(in), OPTIONAL :: lpost         ! into output var_list
    LOGICAL,          INTENT(in), OPTIONAL :: laccu         ! accumulation flag
    REAL(dp),         INTENT(in), OPTIONAL :: reset         ! reset value
    LOGICAL,          INTENT(in), OPTIONAL :: lmiss         ! missing value flag
    REAL(dp),         INTENT(in), OPTIONAL :: missval       ! missing value
    LOGICAL,          INTENT(in), OPTIONAL :: lrestart_cont ! continue on restart
    INTEGER,          INTENT(in), OPTIONAL :: grid_representation
    LOGICAL,          INTENT(in), OPTIONAL :: verbose
    !
    TYPE(list_element), POINTER :: element
    !
    element => find_list_element (ptr_var_list, name)
    IF (ASSOCIATED (element)) THEN
      CALL set_var_metadata (element%field%info,                              &
       name, cf, grib1, grib2, lrestart, lpost, laccu, reset, lmiss, missval, &
       lrestart_cont, grid_representation, verbose)
    ENDIF
  END SUBROUTINE set_var_list_element_metadata
  !------------------------------------------------------------------------------------------------
  !
  ! Set parameters of list element already created
  ! (private routine within this module)
  !
  ! Set each parameter in data type var_metadata if the respective 
  ! optional parameter is present.
  !
  SUBROUTINE set_var_metadata (info, &
       name, cf, grib1, grib2, lrestart, lpost, laccu, reset, lmiss, missval, &
       lrestart_cont, grid_representation, verbose)
    TYPE(var_metadata), INTENT(inout)        :: info          ! memory info struct.
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: name          ! variable name
    TYPE(cf_var),       INTENT(in), OPTIONAL :: cf         
    TYPE(grib1_var),    INTENT(in), OPTIONAL :: grib1
    TYPE(grib2_var),    INTENT(in), OPTIONAL :: grib2
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart      ! restart file flag
    LOGICAL,            INTENT(in), OPTIONAL :: lpost         ! into output var_list
    LOGICAL,            INTENT(in), OPTIONAL :: laccu         ! accumulation flag
    REAL(dp),           INTENT(in), OPTIONAL :: reset         ! reset value
    LOGICAL,            INTENT(in), OPTIONAL :: lmiss         ! missing value flag
    REAL(dp),           INTENT(in), OPTIONAL :: missval       ! missing value
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart_cont ! continue on restart
    INTEGER,            INTENT(in), OPTIONAL :: grid_representation
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
    CALL assign_if_present (info%grib1, grib1)
    CALL assign_if_present (info%grib2, grib2)
    !
    ! set grid type
    !
    CALL assign_if_present (info%grid_representation, grid_representation)
    !
    ! set flags concerning I/O
    !
    CALL assign_if_present (info%lpost,         lpost)
    CALL assign_if_present (info%reset,         reset)
    CALL assign_if_present (info%laccu,         laccu)
    CALL assign_if_present (info%lmiss,         lmiss)
    CALL assign_if_present (info%missval,       missval)
    CALL assign_if_present (info%lrestart,      lrestart)
    CALL assign_if_present (info%lrestart_cont, lrestart_cont)
    !
    ! printout (optional)
    !
    IF (lverbose) CALL print_var_metadata (info)

  END SUBROUTINE set_var_metadata
  !------------------------------------------------------------------------------------------------
  !
  ! Create a list new entry
  !
  ! Specific routines for pointers of different rank
  !
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 4d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_4d (ptr_var_list, &
       name, ptr, ldims, gdims,       &
       klev, ktrac, units, longname, repr,                           &
       lpost, laccu, lmiss, missval, reset, lrestart, contnorest,      &
       table, code, bits, leveltype,                                 &
       dimnames, mem_info, p4, no_default, verbose)
    TYPE(var_list), INTENT(inout)        :: ptr_var_list        ! list
    !
    CHARACTER (len=*)    ,INTENT(in)           :: name          ! name of variable
    REAL(dp)         ,POINTER              :: ptr(:,:,:,:)  ! reference to field
    INTEGER          ,INTENT(in) ,OPTIONAL :: ldims (4)     ! local dimensions
    INTEGER          ,INTENT(in) ,OPTIONAL :: gdims (4)     ! global dimensions
    INTEGER          ,INTENT(in) ,OPTIONAL :: klev          ! number of levels
    INTEGER          ,INTENT(in) ,OPTIONAL :: ktrac         ! number of 'tracers'
    CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: units         ! units
    CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: longname      ! long name
    INTEGER          ,INTENT(in) ,OPTIONAL :: repr          ! representation
    LOGICAL          ,INTENT(in) ,OPTIONAL :: lpost         ! into output var_list
    LOGICAL          ,INTENT(in) ,OPTIONAL :: laccu         ! accumulation flag
    LOGICAL          ,INTENT(in) ,OPTIONAL :: lmiss         ! missing value flag
    REAL(dp)         ,INTENT(in) ,OPTIONAL :: missval       ! missing value
    REAL(dp)         ,INTENT(in) ,OPTIONAL :: reset         ! reset value
    LOGICAL          ,INTENT(in) ,OPTIONAL :: lrestart        ! restart file flag
    LOGICAL          ,INTENT(in) ,OPTIONAL :: contnorest    ! continue
    INTEGER          ,INTENT(in) ,OPTIONAL :: table         ! gribtable number
    INTEGER          ,INTENT(in) ,OPTIONAL :: code          ! gribcode number
    INTEGER          ,INTENT(in) ,OPTIONAL :: bits          ! bits used for GRIB
    INTEGER          ,INTENT(in) ,OPTIONAL :: leveltype     ! grib level type
    CHARACTER (*)    ,INTENT(in) ,OPTIONAL :: dimnames(4)   ! dimension names
    TYPE(var_metadata),POINTER    ,OPTIONAL :: mem_info      ! returned reference
    REAL(dp)         ,POINTER    ,OPTIONAL :: p4(:,:,:,:)   ! provided pointer
    LOGICAL          ,INTENT(in) ,OPTIONAL :: no_default    ! use no defaults
    LOGICAL          ,INTENT(in) ,OPTIONAL :: verbose       ! produce printout
    !
    TYPE(list_element), POINTER :: new_list_element
    LOGICAL                      :: alloc
    !
    ! add list entry
    !
    CALL append_list_element (ptr_var_list, new_list_element, code=code)
    !
    ! and set meta data
    !
    alloc = .NOT. PRESENT(p4)
    new_list_element%field%info%allocated = alloc

!!$    CALL set_var_metadata (new_list_element%field%info                       &
!!$     ,name=name  ,longname=longname ,units=units                              &
!!$     ,ldims=ldims ,gdims=gdims ,ndim=4 ,klev=klev, ktrac=ktrac ,alloc=alloc   &
!!$     ,repr=repr ,lpost=lpost ,laccu=laccu, lmiss=lmiss, missval=missval       &
!!$     ,reset=reset, lrestart=lrestart                                              &
!!$     ,contnorest=contnorest ,table=table ,code=code ,bits=bits                &
!!$     ,leveltype=leveltype ,dimnames=dimnames                                  &
!!$     ,no_default=no_default, verbose=verbose)

    IF (alloc) THEN
      ALLOCATE (new_list_element%field%     ptr(     &
                new_list_element%field%info%used_dimensions(1), &
                new_list_element%field%info%used_dimensions(2), &
                new_list_element%field%info%used_dimensions(3), &
                new_list_element%field%info%used_dimensions(4)))
      ptr_var_list%memory_used = ptr_var_list%memory_used &
                 +8*SIZE(new_list_element%field%ptr)
    ELSE
      new_list_element%field%ptr => p4
    ENDIF
    ptr => new_list_element%field%ptr(:,:,:,:)
    IF(PRESENT(mem_info)) mem_info => new_list_element%field%info

    IF ( PRESENT(lmiss) ) THEN
      new_list_element%field%ptr=new_list_element%field%info%missval
    ELSE
      new_list_element%field%ptr=0._dp
    END IF

  END SUBROUTINE add_var_list_element_4d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 3d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_3d (ptr_var_list, name, ptr, ldims, gdims,       &
             klev, units, longname, repr,                                  &
             lpost, laccu, lmiss, missval, reset, lrestart, contnorest,      &
             table, code, bits, leveltype,                                 &
             dimnames, mem_info, p4, no_default, verbose, tracidx)
  TYPE(var_list)  ,INTENT(inout)        :: ptr_var_list        ! list

  CHARACTER (*)    ,INTENT(in)           :: name          ! name of variable
  REAL(dp)         ,POINTER              :: ptr(:,:,:)    ! reference to field
  INTEGER          ,INTENT(in) ,OPTIONAL :: ldims (3)     ! local dimensions
  INTEGER          ,INTENT(in) ,OPTIONAL :: gdims (3)     ! global dimensions
  INTEGER          ,INTENT(in) ,OPTIONAL :: klev          ! number of levels
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: units         ! units
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: longname      ! long name
  INTEGER          ,INTENT(in) ,OPTIONAL :: repr          ! representation
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lpost         ! into output var_list
  LOGICAL          ,INTENT(in) ,OPTIONAL :: laccu         ! accumulation flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lmiss         ! missing value flag
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: missval       ! missing value
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: reset         ! reset value
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lrestart        ! restart file flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: contnorest    ! continue
  INTEGER          ,INTENT(in) ,OPTIONAL :: table         ! gribtable number
  INTEGER          ,INTENT(in) ,OPTIONAL :: code          ! gribcode number
  INTEGER          ,INTENT(in) ,OPTIONAL :: bits          ! bits used for GRIB
  INTEGER          ,INTENT(in) ,OPTIONAL :: leveltype     ! grib level type
  CHARACTER (*)    ,INTENT(in) ,OPTIONAL :: dimnames(3)   ! dimension names
  TYPE(var_metadata),POINTER    ,OPTIONAL :: mem_info      ! returned reference
  REAL(dp)         ,POINTER    ,OPTIONAL :: p4(:,:,:,:)   ! provided pointer
  LOGICAL          ,INTENT(in) ,OPTIONAL :: no_default    ! use no defaults
  LOGICAL          ,INTENT(in) ,OPTIONAL :: verbose       ! produce printout
  INTEGER          ,INTENT(in) ,OPTIONAL :: tracidx

    TYPE(list_element), POINTER :: new_list_element
    LOGICAL                      :: alloc
    !
    ! add list entry
    !
    CALL append_list_element (ptr_var_list, new_list_element, code=code)
    !
    ! and set meta data
    !

    alloc = .NOT. PRESENT(p4)
    new_list_element%field%info%allocated = alloc

!!$    CALL set_var_metadata (new_list_element%field%info                       &
!!$     ,name=name  ,longname=longname ,units=units                              &
!!$     ,ldims=ldims ,gdims=gdims ,ndim=3 ,klev=klev ,alloc=alloc ,repr=repr     &
!!$     ,lpost=lpost ,laccu=laccu, lmiss=lmiss, missval=missval                  &
!!$     ,reset=reset, lrestart=lrestart                                              &
!!$     ,contnorest=contnorest ,table=table ,code=code ,bits=bits                &
!!$     ,leveltype=leveltype ,dimnames=dimnames                                  &
!!$     ,no_default=no_default, verbose=verbose, tracidx=tracidx)

    IF (alloc) THEN
      ALLOCATE (new_list_element%field%     ptr(      &
                new_list_element%field%info%used_dimensions(1),  &
                new_list_element%field%info%used_dimensions(2),  &
                new_list_element%field%info%used_dimensions(3),1))
      ptr_var_list%memory_used = ptr_var_list%memory_used &
                 +8*SIZE(new_list_element%field%ptr)
    ELSE
      new_list_element%field%ptr => p4
    ENDIF
    ptr => new_list_element%field%ptr(:,:,:,1)
    IF(PRESENT(mem_info)) mem_info => new_list_element%field%info

    IF ( PRESENT(lmiss) ) THEN
      new_list_element%field%ptr=new_list_element%field%info%missval
    ELSE
      new_list_element%field%ptr=0._dp
    END IF

  END SUBROUTINE add_var_list_element_3d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 2d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_2d (ptr_var_list, name, ptr, ldims, gdims,       &
             klev, units, longname, repr,                                  &
             lpost, laccu, lmiss, missval, reset, lrestart, contnorest,      &
             table, code, bits, leveltype,                                 &
             dimnames, no_default, verbose)
  TYPE(var_list)  ,INTENT(inout)        :: ptr_var_list        ! list

  CHARACTER (*)    ,INTENT(in)           :: name          ! name of variable
  REAL(dp)         ,POINTER              :: ptr(:,:)      ! reference to field
  INTEGER          ,INTENT(in) ,OPTIONAL :: ldims (2)     ! local dimensions
  INTEGER          ,INTENT(in) ,OPTIONAL :: gdims (2)     ! global dimensions

  INTEGER          ,INTENT(in) ,OPTIONAL :: klev          ! number of levels
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: units         ! units
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: longname      ! long name

  INTEGER          ,INTENT(in) ,OPTIONAL :: repr          ! representation

  LOGICAL          ,INTENT(in) ,OPTIONAL :: lpost         ! into output var_list
  LOGICAL          ,INTENT(in) ,OPTIONAL :: laccu         ! accumulation flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lmiss         ! missing value flag
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: missval       ! missing value
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: reset         ! reset value
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lrestart        ! restart file flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: contnorest    ! continue

  INTEGER          ,INTENT(in) ,OPTIONAL :: table         ! gribtable number
  INTEGER          ,INTENT(in) ,OPTIONAL :: code          ! gribcode number
  INTEGER          ,INTENT(in) ,OPTIONAL :: bits          ! bits used for GRIB
  INTEGER          ,INTENT(in) ,OPTIONAL :: leveltype     ! grib level type

  CHARACTER (*)    ,INTENT(in) ,OPTIONAL :: dimnames(2)   ! dimension names
  LOGICAL          ,INTENT(in) ,OPTIONAL :: no_default    ! use no defaults

  LOGICAL          ,INTENT(in) ,OPTIONAL :: verbose       ! produce printout

    LOGICAL :: alloc

    TYPE(list_element), POINTER :: element

    CALL append_list_element (ptr_var_list, element, code=code)

    alloc = .TRUE.

!!$    CALL set_var_metadata (element%field%info                                &
!!$     ,name=name  ,longname=longname ,units=units                              &
!!$     ,ldims=ldims ,gdims=gdims ,ndim=2 ,klev=klev, alloc=alloc ,repr=repr     &
!!$     ,lpost=lpost ,laccu=laccu, lmiss=lmiss, missval=missval                  &
!!$     ,reset=reset, lrestart=lrestart                                              &
!!$     ,contnorest=contnorest ,table=table ,code=code ,bits=bits                &
!!$     ,leveltype=leveltype ,dimnames=dimnames                                  &
!!$     ,no_default=no_default, verbose=verbose)

    ALLOCATE (element%field%ptr (             &
              element%field%info%used_dimensions(1),    &
              element%field%info%used_dimensions(2),1,1))
    ptr_var_list%memory_used = ptr_var_list%memory_used &
               +8*SIZE(element%field%ptr)
    ptr => element%field%ptr(:,:,1,1)

    IF ( PRESENT(lmiss) ) THEN
      element%field%ptr = element%field%info%missval
    ELSE
      element%field%ptr=0._dp
    END IF

  END SUBROUTINE add_var_list_element_2d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 1d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_1d (ptr_var_list, name, ptr, ldims, gdims, &
             klev, units, longname, repr,                                  &
             lpost, laccu, lmiss, missval, reset, lrestart, contnorest,      &
             table, code, bits, leveltype,                                 &
             dimnames, no_default, verbose)

    TYPE(var_list)       ,INTENT(inout):: ptr_var_list      ! list
    CHARACTER (*)         ,INTENT(in)   :: name        ! name of variable
    REAL(dp)              ,POINTER      :: ptr(:)      ! reference to field
    INTEGER      ,OPTIONAL,INTENT(in)   :: ldims(1)    ! shape of array 
    INTEGER      ,OPTIONAL,INTENT(in)   :: gdims(1)    ! global size of field

    INTEGER      ,OPTIONAL,INTENT(in)   :: klev        ! number of levels
    CHARACTER(*) ,OPTIONAL,INTENT(in)   :: units       ! units
    CHARACTER(*) ,OPTIONAL,INTENT(in)   :: longname    ! long name

    INTEGER      ,OPTIONAL,INTENT(in)   :: repr        ! representation

    LOGICAL      ,OPTIONAL,INTENT(in)   :: lpost       ! into output var_list
    LOGICAL      ,OPTIONAL,INTENT(in)   :: lmiss       ! missing value flag
    REAL(dp)     ,OPTIONAL,INTENT(in)   :: missval     ! missing value
    REAL(dp)     ,OPTIONAL,INTENT(in)   :: reset       ! reset value

    INTEGER      ,OPTIONAL,INTENT(in)   :: bits        ! bits used for GRIB
    INTEGER      ,OPTIONAL,INTENT(in)   :: leveltype   ! grib level type

    LOGICAL      ,OPTIONAL,INTENT(in)   :: no_default  ! use no defaults

    LOGICAL      ,OPTIONAL,INTENT(in)   :: verbose     ! produce printout
    CHARACTER(*) ,OPTIONAL,INTENT(in)   :: dimnames(1) ! dimension names
    INTEGER      ,OPTIONAL,INTENT(in)   :: code        ! gribcode number
    INTEGER      ,OPTIONAL,INTENT(in)   :: table       ! gribcode table number
    LOGICAL      ,OPTIONAL,INTENT(in)   :: laccu       ! accumulation flag
    LOGICAL      ,OPTIONAL,INTENT(in)   :: lrestart      ! restart file flag
    LOGICAL      ,OPTIONAL,INTENT(in)   :: contnorest  ! continue on restart
!    TYPE(list_element), POINTER :: new_list_element
!    !
!    ! add list entry
!    !
!    CALL append_list_element (var_list, new_list_element, code=code)
!    !
!    ! and set meta data
!    !
!    new_list_element%field%info%name = name
!    ALLOCATE (new_list_element%field%ptr(ldims(1),1,1,1))
!
!    new_list_element%field%ptr=0._dp
!
!    var_list%memory_used = var_list%memory_used &
!               +8*SIZE(new_list_element%field%ptr)
!
!    new_list_element%field%info%alloc = .TRUE.
!
!    new_list_element%field%info%dim(1) = ldims(1)
!    new_list_element%field%info%dim(2) = 1
!    new_list_element%field%info%dim(3) = 1
!    new_list_element%field%info%dim(4) = 1
!    new_list_element%field%info%dima   = new_list_element%field%info%dim
!    new_list_element%field%info%lreg   = .true.
!
!    new_list_element%field%info%gdim(1) = gdims(1)
!    new_list_element%field%info%gdim(2) = 1
!    new_list_element%field%info%gdim(3) = 1
!    new_list_element%field%info%gdim(4) = 1
!
!    new_list_element%field%info%ndim = 1
!
!    ptr => new_list_element%field%ptr(:,1,1,1)
!    !
!    ! pass optional arguments
!    !
!    IF (PRESENT(dimnames)) THEN
!       new_list_element%field%info%IO_var_indx(1) = IO_get_varindx(dimnames(1))
!    END IF
!
!    IF (PRESENT(table))     new_list_element%field%info%gribtable  = table 
!    IF (PRESENT(code))      new_list_element%field%info%gribcode   = code
!    IF (PRESENT(laccu))     new_list_element%field%info%laccu      = laccu
!    IF (PRESENT(lrestart))    new_list_element%field%info%lrestart     = lrestart
!    IF (PRESENT(contnorest))new_list_element%field%info%contnorest = contnorest

    LOGICAL :: alloc

    TYPE(list_element), POINTER :: element

    CALL append_list_element (ptr_var_list, element, code=code)

    alloc = .TRUE.
!!$
!!$    CALL set_var_metadata (element%field%info                                &
!!$     ,name=name  ,longname=longname ,units=units                              &
!!$     ,ldims=ldims ,gdims=gdims ,ndim=1 ,klev=klev, alloc=alloc ,repr=repr     &
!!$     ,lpost=lpost ,laccu=laccu, lmiss=lmiss, missval=missval                  &
!!$     ,reset=reset, lrestart=lrestart                                              &
!!$     ,contnorest=contnorest ,table=table ,code=code ,bits=bits                &
!!$     ,leveltype=leveltype ,dimnames=dimnames                                  &
!!$     ,no_default=no_default, verbose=verbose)

    ALLOCATE (element%field%ptr (             &
              element%field%info%used_dimensions(1),    &
              1,1,1))
    ptr_var_list%memory_used = ptr_var_list%memory_used &
               +8*SIZE(element%field%ptr)
    ptr => element%field%ptr(:,1,1,1)

    IF ( PRESENT(lmiss) ) THEN
      element%field%ptr=element%field%info%missval
    ELSE
      element%field%ptr=0._dp
    END IF

  END SUBROUTINE add_var_list_element_1d
  !------------------------------------------------------------------------------------------------
  !
  ! remove one element from the list
  ! the element is identified by its name
  !
  SUBROUTINE remove_var_list_element (this_list, name)
    TYPE(var_list),   INTENT(inout) :: this_list
    CHARACTER(len=*), INTENT(in)    :: name   
    TYPE(list_element), POINTER :: ptr

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
  SUBROUTINE get_var_list_element_4d (ptr_var_list, name, ptr)
    TYPE(var_list),   INTENT(in) :: ptr_var_list ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    REAL(dp),         POINTER    :: ptr(:,:,:,:) ! reference to allocated field

    TYPE(list_element), POINTER :: element

    element => find_list_element (ptr_var_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%ptr
    
  END SUBROUTINE get_var_list_element_4d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 3d-field
  !
  SUBROUTINE get_var_list_element_3d (ptr_var_list, name, ptr)

    TYPE(var_list) ,INTENT(in) :: ptr_var_list       ! list
    CHARACTER (*)   ,INTENT(in) :: name         ! name of variable
    REAL(dp)        ,POINTER    :: ptr(:,:,:)   ! reference to allocated field

    TYPE(list_element), POINTER :: element
    element => find_list_element (ptr_var_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%ptr(:,:,:,1)

  END SUBROUTINE get_var_list_element_3d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 2d-field
  !
  SUBROUTINE get_var_list_element_2d (ptr_var_list, name, ptr)

    TYPE(var_list)   ,INTENT(in) :: ptr_var_list     ! list
    CHARACTER (*)     ,INTENT(in) :: name       ! name of variable
    REAL(dp)          ,POINTER    :: ptr(:,:)   ! reference to allocated field

    TYPE(list_element), POINTER :: element

    element => find_list_element (ptr_var_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%ptr(:,:,1,1)

  END SUBROUTINE get_var_list_element_2d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 1d-field
  !
  SUBROUTINE get_var_list_element_1d (ptr_var_list, name, ptr)

    TYPE(var_list) ,INTENT(in) :: ptr_var_list    ! list
    CHARACTER (*)   ,INTENT(in) :: name      ! name of variable
    REAL(dp)        ,POINTER    :: ptr(:)    ! reference to allocated field

    TYPE(list_element), POINTER :: element
    element => find_list_element (ptr_var_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%ptr(:,1,1,1)

  END SUBROUTINE get_var_list_element_1d
  !------------------------------------------------------------------------------------------------
  !
  ! Print routines for control output and debuggung
  !
  SUBROUTINE print_memory_use (ptr_var_list)

    TYPE(var_list) ,INTENT(in) :: ptr_var_list ! list

    WRITE (message_text,'(a16,a,a,i10,a,i4,a)')                 &
         TRIM(ptr_var_list%name), '-buffer: ',                        &
         'Memory in use: ', ptr_var_list%memory_used/1024, ' kb in ', &
         ptr_var_list%list_elements, ' fields.'
    CALL message('',message_text)

  END SUBROUTINE print_memory_use
  !------------------------------------------------------------------------------------------------
  !
  ! print current memory table 
  !
  SUBROUTINE print_memory_table (ptr_var_list)
    TYPE(var_list),  INTENT(in) :: ptr_var_list ! list
    !
    CALL message('','')
    CALL message('','')
    CALL message('','Status of base memory:')    
    CALL message('','')
    !
    CALL print_linked_list (ptr_var_list)
    !    
  END SUBROUTINE print_memory_table
  !------------------------------------------------------------------------------------------------
  !
  ! print current stat table 
  !
  SUBROUTINE print_sinfo (ptr_var_list)
    TYPE(var_list),  INTENT(in) :: ptr_var_list
    WRITE (message_text,'(a16,a)') TRIM(ptr_var_list%name), '-buffer: '
    CALL message('',message_text)
    CALL message('','')    
    CALL message('','')
    CALL message('','Statistic of base memory:')
    CALL message('','')
    !
    CALL print_sinfo_list (ptr_var_list)
    !
  END SUBROUTINE print_sinfo
  !------------------------------------------------------------------------------------------------
  !
  ! add supplementary fields (eg. geopotential, surface pressure, gridbox area)
  !
  SUBROUTINE add_var_list_reference (to_var_list, name, from_var_list, lpost, kprec)
    TYPE(var_list)  ,INTENT(inout)            :: to_var_list
    CHARACTER(len=*) ,INTENT(in)              :: name
    CHARACTER(len=*) ,INTENT(in)    ,OPTIONAL :: from_var_list
    LOGICAL          ,INTENT(in)    ,OPTIONAL :: lpost
    INTEGER          ,INTENT(in)    ,OPTIONAL :: kprec
    !
    TYPE(var_list_element)  ,POINTER :: source
    TYPE(list_element) ,POINTER :: new_list_element
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
      TYPE(var_list_element), POINTER        :: element
      CHARACTER(len=*), INTENT(in)           :: name
      CHARACTER(len=*), INTENT(in), OPTIONAL :: in_var_list
      !
      INTEGER                     :: i
      TYPE(list_element), POINTER :: link
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
    TYPE(cf_var), INTENT(inout)        :: y
    TYPE(cf_var), INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    y = x
  END SUBROUTINE assign_if_present_cf
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_grib1 (y,x)
    TYPE(grib1_var), INTENT(inout)        :: y
    TYPE(grib1_var) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    y = x
  END SUBROUTINE assign_if_present_grib1
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_grib2 (y,x)
    TYPE(grib2_var), INTENT(inout)        :: y
    TYPE(grib2_var) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    y = x
  END SUBROUTINE assign_if_present_grib2
  !------------------------------------------------------------------------------------------------
END MODULE mo_var_list
