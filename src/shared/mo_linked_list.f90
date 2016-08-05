!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
#if  (defined(__SX__) || defined(__SUNPRO_F95) || defined (__GNUC__)) 
#define HAVE_F95
#endif
MODULE mo_linked_list
  !
  ! This is a specific linked list implementation for handling ICON output.
  ! When Fortran 2003 is available on almost all production machines this 
  ! should be replaced by a proper generic version.
  !
  !----------------------------------------------------------------------------
  !
  USE mo_kind,             ONLY: i8
  USE mo_exception,        ONLY: finish, message
  USE mo_var_list_element, ONLY: t_var_list_element, level_type_ml
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: t_var_list          ! anchor for a whole list
  PUBLIC :: t_list_element
  !
  PUBLIC :: new_list            ! construct an (empty) list
  PUBLIC :: delete_list         ! clean up the list
  !
  PUBLIC :: append_list_element ! add an element to the list
  PUBLIC :: delete_list_element ! remove one element from the list
  !
#ifdef HAVE_F95
  PUBLIC :: t_var_list_intrinsic
#endif
  !
  ! t_list_element provides the entry to the actual information 
  ! and a reference to the next element in the list
  !  
  TYPE t_list_element
    TYPE(t_var_list_element)      :: field
    TYPE(t_list_element), POINTER :: next_list_element
  END TYPE t_list_element
  !
  TYPE t_var_list_intrinsic
    INTEGER                       :: key                ! hash value of name   
    CHARACTER(len=128)            :: name               ! stream name
    TYPE(t_list_element), POINTER :: first_list_element ! reference to first
    INTEGER(i8)                   :: memory_used        ! memory allocated
    INTEGER                       :: list_elements      ! allocated elements
    LOGICAL                       :: loutput            ! output stream
    LOGICAL                       :: lrestart           ! restart stream
    LOGICAL                       :: linitial           ! initial stream
    CHARACTER(len=256)            :: filename           ! name of file
    CHARACTER(len=8)              :: post_suf           ! suffix of output  file
    CHARACTER(len=8)              :: rest_suf           ! suffix of restart file
    CHARACTER(len=8)              :: init_suf           ! suffix of initial file
    LOGICAL                       :: first              ! first var_list in file
    INTEGER                       :: output_type        ! CDI format
    INTEGER                       :: restart_type       ! CDI format
    INTEGER                       :: compression_type   ! CDI compression type
    LOGICAL                       :: restart_opened     ! true, if restart file opened
    LOGICAL                       :: output_opened      ! true, if output file opened
    CHARACTER(len=8)              :: model_type         ! store model type
    INTEGER                       :: patch_id           ! ID of patch to which list variables belong
    INTEGER                       :: vlevel_type        ! 1: model levels, 2: pressure levels, 3: height levels
    !--------------------------------------------------------------------------------------------
    ! Internal used handler for CDI setup of synchronous restart
    !
    ! Todo: This metadata should not be placed in this location ?!
    INTEGER                       :: cdiTimeIndex
    !
    INTEGER                       :: nvars

    ! Metadata for missing value masking

    LOGICAL                    :: lmiss          ! flag: true, if variables should be initialized with missval
    LOGICAL                    :: lmask_boundary ! flag: true, if interpolation zone should be masked *in output*
  END TYPE t_var_list_intrinsic
  !
  TYPE t_var_list
    TYPE(t_var_list_intrinsic), POINTER :: p
  END type t_var_list

  !
CONTAINS
  !
  !-----------------------------------------------------------------------------
  !
  ! initialize a variable of type var_list with default values and
  ! nullify anchor to linked list
  !
  SUBROUTINE new_list(this_list)
    TYPE(t_var_list), INTENT(inout) :: this_list
    !
    ALLOCATE(this_list%p)
    !
    this_list%p%key                = 0
    this_list%p%name               = ''
    this_list%p%first_list_element => NULL()
    this_list%p%memory_used        = 0_i8
    this_list%p%list_elements      = 0
    !
    this_list%p%loutput            = .FALSE.
    this_list%p%lmiss              = .FALSE.
    this_list%p%lmask_boundary     = .TRUE.
    this_list%p%lrestart           = .FALSE.
    this_list%p%linitial           = .FALSE.
    !
    this_list%p%filename           = ''
    !
    this_list%p%post_suf           = ''
    this_list%p%rest_suf           = ''
    this_list%p%init_suf           = ''
    !
    this_list%p%first              = .FALSE.
    !
    this_list%p%output_type        = -1
    this_list%p%restart_type       = -1
    this_list%p%compression_type   = -1
    !
    this_list%p%restart_opened     = .FALSE.
    this_list%p%output_opened      = .FALSE.
    this_list%p%model_type         = 'atm'
    this_list%p%patch_id           = -1
    this_list%p%vlevel_type        =  level_type_ml ! Default is model levels
    !
    this_list%p%cdiTimeIndex       = -1
    !
    this_list%p%nvars              = 0
  END SUBROUTINE new_list
  !-----------------------------------------------------------------------------
  !
  ! remove all elements of a linked list
  ! check if all elements are removed
  !
  SUBROUTINE delete_list(this_list)
    TYPE(t_var_list), INTENT(inout) :: this_list
    !
    CALL delete_list_elements(this_list, this_list%p%first_list_element)
    !
    this_list%p%first_list_element => NULL()
    !
    IF (this_list%p%memory_used /= 0_i8) THEN
      CALL finish ('delete_list', 'List delete didnt work proper (memory counter)')
    ENDIF
    !
    IF (this_list%p%list_elements /= 0) THEN
      CALL finish ('delete_list', 'List delete didnt work proper (element counter)')
    ENDIF
    !
  END SUBROUTINE delete_list
  !-----------------------------------------------------------------------------
  !
  ! deallocate a list element and all its sucessors, private routine
  !
  SUBROUTINE delete_list_elements(this_list, this_list_element)
    !
    TYPE(t_var_list),     INTENT(inout) :: this_list
    TYPE(t_list_element), POINTER       :: this_list_element
    !
    TYPE(t_list_element), POINTER       :: this, next
    !
    next => this_list_element
    this_list_element => NULL()
    !
    DO
      IF (.NOT. ASSOCIATED(next)) EXIT
      this => next
      next => this%next_list_element
      !
      IF (this%field%info%allocated) THEN
        IF (ASSOCIATED(this%field%r_ptr)) THEN
          this_list%p%memory_used = this_list%p%memory_used &
               &                   -INT(this%field%var_base_size*SIZE(this%field%r_ptr),i8)
          DEALLOCATE (this%field%r_ptr)
        ELSE IF (ASSOCIATED(this%field%i_ptr)) THEN
          this_list%p%memory_used = this_list%p%memory_used &
               &                   -INT(this%field%var_base_size*SIZE(this%field%i_ptr),i8)
          DEALLOCATE (this%field%i_ptr)
        ELSE IF (ASSOCIATED(this%field%l_ptr)) THEN
          this_list%p%memory_used = this_list%p%memory_used &
               &                   -INT(this%field%var_base_size*SIZE(this%field%l_ptr),i8)
          DEALLOCATE (this%field%l_ptr)
        ENDIF
        this%field%info%allocated = .FALSE.
      ENDIF
      this_list%p%list_elements = this_list%p%list_elements-1
      DEALLOCATE (this)
     END DO
    !
  END SUBROUTINE delete_list_elements
  !-----------------------------------------------------------------------------
  SUBROUTINE create_list_element (this_list, current_list_element)
    !
    TYPE(t_var_list),     INTENT(inout) :: this_list
    TYPE(t_list_element), POINTER       :: current_list_element
    !
    INTEGER :: ist
    !
    ALLOCATE (current_list_element, STAT=ist)
    IF (ist /= 0) THEN
      CALL finish('create_list_element','Cannot add element to linked list ...')
    ENDIF
    this_list%p%list_elements = this_list%p%list_elements+1
    !
    current_list_element%next_list_element => NULL()
    current_list_element%field%r_ptr       => NULL()
    current_list_element%field%i_ptr       => NULL()
    current_list_element%field%l_ptr       => NULL()
    !
  END SUBROUTINE create_list_element
  !-----------------------------------------------------------------------------
  !
  ! add a list element to the linked list
  !
  SUBROUTINE append_list_element (this_list, new_list_element)
    !
    TYPE(t_var_list),     INTENT(inout) :: this_list
    TYPE(t_list_element), POINTER       :: new_list_element
    !
    TYPE(t_list_element), POINTER :: current_list_element
    !
    ! insert as first element if list is empty
    !
    IF (.NOT. ASSOCIATED (this_list%p%first_list_element)) THEN
      CALL create_list_element (this_list, this_list%p%first_list_element)
      new_list_element => this_list%p%first_list_element
      RETURN
    ENDIF
    !
    ! loop over list elements to find position
    !
    current_list_element => this_list%p%first_list_element
    DO WHILE (ASSOCIATED(current_list_element%next_list_element)) 
      current_list_element => current_list_element%next_list_element
    ENDDO
    !
    ! insert element
    !
    CALL create_list_element (this_list, new_list_element)
    new_list_element%next_list_element => current_list_element%next_list_element
    current_list_element%next_list_element => new_list_element
    !
    this_list%p%nvars = this_list%p%nvars + 1
  END SUBROUTINE append_list_element
  !-----------------------------------------------------------------------------
  SUBROUTINE delete_list_element (this_list, delete_this_list_element)
    !
    TYPE(t_var_list),     INTENT(inout) :: this_list
    TYPE(t_list_element), POINTER       :: delete_this_list_element
    !
    TYPE(t_list_element), POINTER :: current_list_element
    !
    IF (ASSOCIATED(delete_this_list_element, this_list%p%first_list_element)) THEN
      this_list%p%first_list_element => delete_this_list_element%next_list_element
    ELSE
      current_list_element => this_list%p%first_list_element
      DO WHILE ((ASSOCIATED(current_list_element)) &
           &           .AND. (.NOT. ASSOCIATED(current_list_element%next_list_element, &
           &           delete_this_list_element)))
        current_list_element => current_list_element%next_list_element
      ENDDO
      IF (.NOT. ASSOCIATED(current_list_element)) THEN
        CALL message('', 'Cannot find element to be deleted ...')
        RETURN
      ENDIF
      current_list_element%next_list_element &
           &          => current_list_element%next_list_element%next_list_element
    ENDIF
    !
    IF (delete_this_list_element%field%info%allocated) THEN
      IF (ASSOCIATED(delete_this_list_element%field%r_ptr)) THEN
        this_list%p%memory_used = this_list%p%memory_used                          &
             &                   -INT(delete_this_list_element%field%var_base_size &
             &                   *SIZE(delete_this_list_element%field%r_ptr),i8)
        DEALLOCATE (delete_this_list_element%field%r_ptr)
      ELSE IF (ASSOCIATED(delete_this_list_element%field%i_ptr)) THEN
        this_list%p%memory_used = this_list%p%memory_used                          &
             &                   -INT(delete_this_list_element%field%var_base_size &
             &                   *SIZE(delete_this_list_element%field%i_ptr),i8)
        DEALLOCATE (delete_this_list_element%field%i_ptr)
      ELSE IF (ASSOCIATED(delete_this_list_element%field%l_ptr)) THEN
        this_list%p%memory_used = this_list%p%memory_used                          &
             &                   -INT(delete_this_list_element%field%var_base_size &
             &                   *SIZE(delete_this_list_element%field%l_ptr),i8)
        DEALLOCATE (delete_this_list_element%field%l_ptr)
      ENDIF
      delete_this_list_element%field%info%allocated = .FALSE.
    ENDIF
    !
    this_list%p%list_elements = this_list%p%list_elements-1
    DEALLOCATE (delete_this_list_element)
    !
    this_list%p%nvars = this_list%p%nvars - 1
  END SUBROUTINE delete_list_element
  !-----------------------------------------------------------------------------
END MODULE mo_linked_list
