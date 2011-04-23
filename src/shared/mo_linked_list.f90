MODULE mo_linked_list
  !
  ! This is a specific linked list implementation for handling ICON output.
  ! When Fortran 2003 is availbale on almost all production machines this 
  ! should be replaced by a proper generic version.
  !
  ! Authors:
  !
  ! Luis Kornblueh, MPI,             original code
  ! Andreas Rhodin, MPI, April 2001, extended and documented
  ! Luis Kornblueh, MPI, April 2011, adapted for ICON
  !----------------------------------------------------------------------------
  !
  USE mo_kind,             ONLY: dp, i8
  USE mo_exception,        ONLY: message_text, finish, message
  USE mo_util_string,      ONLY: separator
  USE mo_var_list_element, ONLY: var_list_element
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: var_list            ! anchor for a whole list
  PUBLIC :: list_element
  !
  PUBLIC :: new_list            ! construct an (empty) list
  PUBLIC :: delete_list         ! clean up the list
  !
  PUBLIC :: append_list_element ! add an element to the list
  PUBLIC :: remove_list_element ! remove one element from the list
  PUBLIC :: find_list_element   ! find an element in the list
  !
  ! Type list_element provides the entry to the actual information 
  ! and a reference to the next element in the list.
  !  
  TYPE list_element
    TYPE(var_list_element)      :: field
    TYPE(list_element), POINTER :: next_list_element
  END TYPE list_element
  !
  TYPE var_list
    INTEGER                     :: key                ! hash value of name   
    CHARACTER(len=16)           :: name               ! stream name
    TYPE(list_element), POINTER :: first_list_element ! reference to first
    INTEGER(i8)                 :: memory_used        ! memory allocated
    INTEGER                     :: list_elements      ! allocated elements
    !
    ! Postprocessing information
    !
    LOGICAL                     :: lpost           ! postprocessing stream
    LOGICAL                     :: lrestart        ! restart stream
    LOGICAL                     :: linitial        ! initial stream
    !
    CHARACTER(len=128)          :: filename        ! name of file
    !
    CHARACTER(len=8)            :: post_suf        ! suffix of output  file
    CHARACTER(len=8)            :: rest_suf        ! suffix of restart file
    CHARACTER(len=8)            :: init_suf        ! suffix of initial file
    !
    ! CDI support and handler
    !
    LOGICAL                     :: first           ! first stream in file
    !
    INTEGER                     :: output_type      ! CDI format
    INTEGER                     :: restart_type     ! CDI format
    INTEGER                     :: compression_type ! CDI comp. type
    !
    INTEGER                     :: fileID           ! CDI file ID
    INTEGER                     :: vlistID          ! CDI vlist ID
    !
    INTEGER                     :: timestep         ! Output timestep
    !
  END TYPE var_list
  !
CONTAINS
  !
  !-----------------------------------------------------------------------------
  !
  ! initialize a variable of type var_list with default values and
  ! nullify anchor to linked lisat
  !
  SUBROUTINE new_list(this_list)
    TYPE(var_list), INTENT(out) :: this_list
    !
    this_list%key                = 0
    this_list%name               = ''
    this_list%first_list_element => NULL()
    this_list%memory_used        = 0
    this_list%list_elements      = 0
    !
    this_list%lpost              = .TRUE.
    this_list%lrestart           = .TRUE.
    this_list%linitial           = .TRUE.
    !
    this_list%filename           = ''
    !
    this_list%post_suf           = ''
    this_list%rest_suf           = ''
    this_list%init_suf           = ''
    !
    this_list%first              = .FALSE.
    !
    this_list%output_type        = -1
    this_list%restart_type       = -1
    this_list%compression_type   = -1
    !
    this_list%fileID             = -1
    this_list%vlistID            = -1
    !
    this_list%timestep           = -1
    !
  END SUBROUTINE new_list
  !-----------------------------------------------------------------------------
  !
  ! remove all elements of a linked list
  ! check if all elements are removed
  !
  SUBROUTINE delete_list(this_list)
    TYPE(var_list), INTENT(inout) :: this_list

    CALL delete_list_element(this_list, this_list%first_list_element)

    this_list%first_list_element => NULL()
    
    IF (this_list%memory_used /= 0_i8) THEN
      CALL message('', &
           'List delete didn''t work properly (memory counter) ...')
      CALL finish ('delete_list',&
           'List delete didnt work proper (memory counter)')
    ENDIF

    IF (this_list%list_elements /= 0) THEN
      CALL message('', &
           'List delele didn''t work proper (element counter) ...')
      CALL finish ('delete_list',&
           'List delete didnt work proper (element counter)')
    ENDIF

  END SUBROUTINE delete_list
  !-----------------------------------------------------------------------------
  !
  ! deallocate a list element and all its sucessors
  !
  SUBROUTINE delete_list_element(this_list, this_list_element)
    TYPE(var_list),     INTENT(inout) :: this_list
    TYPE(list_element), POINTER       :: this_list_element

    TYPE(list_element), POINTER       :: this, next

    next => this_list_element
    this_list_element => NULL()
    DO
      IF (.NOT. ASSOCIATED(next)) EXIT
      this => next
      next => this%next_list_element
      !
      ! 8 as constant has to be adjusted with information from mo_machine
      ! the variable to be used is mp_real8
      !
      IF (this%field%info%allocated) THEN
        this_list%memory_used = this_list%memory_used-8*SIZE(this%field%ptr)
        DEALLOCATE (this%field%ptr)
        this%field%info%allocated = .FALSE.
      ENDIF
      this_list%list_elements = this_list%list_elements-1
      DEALLOCATE (this)
    END DO
    
  END SUBROUTINE delete_list_element
  !-----------------------------------------------------------------------------
  !
  ! remove an element from the list which is identified by its pointer
  !
  SUBROUTINE remove_list_element(this_list, this_list_element)
    TYPE(var_list),     INTENT(inout) :: this_list
    TYPE(list_element), POINTER       :: this_list_element

    TYPE(list_element) ,POINTER :: p

    p => this_list_element
    this_list_element => this_list_element%next_list_element

    IF (p%field%info%allocated) THEN
      this_list%memory_used = this_list%memory_used-8*SIZE(p%field%ptr)
      DEALLOCATE (p%field%ptr)
      p%field%info%allocated = .FALSE.
    ENDIF
    this_list%list_elements = this_list%list_elements-1
    DEALLOCATE (p) ! LK necessary?
    
  END SUBROUTINE remove_list_element
  !-----------------------------------------------------------------------------
  SUBROUTINE create_list_element (this_list, current_list_element)
    TYPE(var_list),     INTENT(inout) :: this_list
    TYPE(list_element), POINTER       :: current_list_element

    INTEGER :: ist

    ALLOCATE (current_list_element, STAT=ist)
    IF (ist /= 0) THEN
      CALL message('', 'Cannot add element to linked list ...')
      CALL finish('create_list_element','Cannot add element to linked list')
    ENDIF
    this_list%list_elements = this_list%list_elements+1
    
    current_list_element%next_list_element => NULL()
    current_list_element%field%ptr         => NULL()
    
  END SUBROUTINE create_list_element
  !-----------------------------------------------------------------------------
  !
  ! Add a list element to the linked list, if code is present, order with 
  ! respect to the grib code else search for end of the list to append
  !
  SUBROUTINE append_list_element (this_list, new_list_element, code)
    TYPE(var_list),      INTENT(inout) :: this_list
    TYPE(list_element),  POINTER       :: new_list_element
    INTEGER, OPTIONAL,   INTENT(in)    :: code

    TYPE(list_element), POINTER :: current_list_element
    INTEGER                      :: ic

    ic = 1000
    IF (PRESENT(code)) ic = code    
    !
    ! insert as first element if list is empty
    !
    IF (.NOT. ASSOCIATED (this_list%first_list_element)) THEN
      CALL create_list_element (this_list, this_list%first_list_element)
      new_list_element => this_list%first_list_element
      RETURN
    ENDIF
    !
    ! insert before first element
    !
    IF (this_list%first_list_element%field%info%grib2%parameter > ic) THEN
      CALL create_list_element (this_list, new_list_element)
      new_list_element%next_list_element => this_list%first_list_element
      this_list%      first_list_element => new_list_element
      RETURN
    ENDIF
    !
    ! loop over list elements to find position
    !
    current_list_element => this_list%first_list_element
    DO WHILE (ASSOCIATED(current_list_element%next_list_element)) 
      IF (current_list_element%next_list_element%field%info%grib2%parameter > ic) &
           EXIT
      current_list_element => current_list_element%next_list_element
    ENDDO
    !
    ! insert element
    !
    CALL create_list_element (this_list, new_list_element)
    new_list_element%    next_list_element => current_list_element%&
         next_list_element
    current_list_element%next_list_element => new_list_element

  END SUBROUTINE append_list_element
  !-----------------------------------------------------------------------------
  SUBROUTINE delete_one_list_element (this_list, delete_this_list_element)
    TYPE(var_list),     INTENT(inout) :: this_list
    TYPE(list_element), POINTER       :: delete_this_list_element

    TYPE(list_element), POINTER :: current_list_element

    IF (ASSOCIATED(delete_this_list_element, &
         this_list%first_list_element)) THEN
      this_list%first_list_element &
           => delete_this_list_element%next_list_element
    ELSE
      current_list_element => this_list%first_list_element
      DO WHILE ((ASSOCIATED(current_list_element)) &
           .AND. (.NOT. ASSOCIATED(current_list_element%next_list_element, &
           delete_this_list_element)))
        current_list_element => current_list_element%next_list_element
      ENDDO
      IF (.NOT. ASSOCIATED(current_list_element)) THEN
        CALL message('', 'Cannot find element to be deleted ...')
        RETURN
      ENDIF
      current_list_element%next_list_element &
           => current_list_element%next_list_element%next_list_element
    ENDIF
    
    IF (delete_this_list_element%field%info%allocated) THEN
      this_list%memory_used = this_list%memory_used &
           -8*SIZE(delete_this_list_element%field%ptr)
      DEALLOCATE (delete_this_list_element%field%ptr)
      delete_this_list_element%field%info%allocated = .FALSE.
    ENDIF
    
    this_list%list_elements = this_list%list_elements-1
    DEALLOCATE (delete_this_list_element)
    
  END SUBROUTINE delete_one_list_element
  !-----------------------------------------------------------------------------
  !
  ! Should be overloaded to be able to search for the different information 
  ! In the proposed structure for the linked list, in the example only
  ! A character string is used so it is straight forward only one find
  !
  FUNCTION find_list_element (this_list, name) RESULT(this_list_element)

    TYPE(var_list),   INTENT(in) :: this_list
    CHARACTER(len=*), INTENT(in) :: name
    
    TYPE(list_element), POINTER :: this_list_element
    INTEGER :: KEY
    
    this_list_element => this_list%first_list_element
    DO WHILE (ASSOCIATED(this_list_element))
      IF (name == this_list_element%field%info%name) THEN
        RETURN
      ENDIF
      this_list_element => this_list_element%next_list_element
    ENDDO
    
    NULLIFY (this_list_element)
    
  END FUNCTION find_list_element
  !-----------------------------------------------------------------------------
END MODULE mo_linked_list
