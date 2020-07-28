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
  !----------------------------------------------------------------------------
  USE mo_kind,             ONLY: i8
  USE mo_var_list_element, ONLY: t_var_list_element, level_type_ml
  !
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC :: t_var_list          ! anchor for a whole list
  PUBLIC :: t_list_element
  PUBLIC :: delete_list         ! clean up the list
  PUBLIC :: append_list_element ! add an element to the list
#ifdef HAVE_F95
  PUBLIC :: t_var_list_intrinsic
#endif

  !
  ! t_list_element provides the entry to the actual information 
  ! and a reference to the next element in the list
  TYPE t_list_element
    TYPE(t_var_list_element)      :: field
    TYPE(t_list_element), POINTER :: next_list_element => NULL()
  END TYPE t_list_element
  !
  TYPE t_var_list_intrinsic
    INTEGER                       :: key = 0            ! hash value of name   
    CHARACTER(len=128)            :: name = ''          ! stream name
    TYPE(t_list_element), POINTER :: first_list_element => NULL() ! reference to first
    INTEGER(i8)                   :: memory_used = 0_i8 ! memory allocated
    INTEGER                       :: list_elements = 0  ! allocated elements
    LOGICAL                       :: loutput = .FALSE.  ! output stream
    LOGICAL                       :: lrestart = .FALSE. ! restart stream
    LOGICAL                       :: linitial = .FALSE. ! initial stream
    CHARACTER(len=256)            :: filename = ''      ! name of file
    CHARACTER(len=8)              :: post_suf = ''      ! suffix of output  file
    CHARACTER(len=8)              :: rest_suf = ''      ! suffix of restart file
    CHARACTER(len=8)              :: init_suf = ''      ! suffix of initial file
    LOGICAL                       :: first              ! first var_list in file
    INTEGER                       :: output_type = -1   ! CDI format
    INTEGER                       :: restart_type = -1  ! CDI format
    INTEGER                       :: compression_type = -1 ! CDI compression type
    LOGICAL                       :: restart_opened = .FALSE. ! true, if restart file opened
    LOGICAL                       :: output_opened = .FALSE. ! true, if output file opened
    CHARACTER(len=8)              :: model_type = ''    ! store model type
    INTEGER                       :: patch_id = -1      ! ID of patch to which list variables belong
    INTEGER                       :: vlevel_type = level_type_ml ! 1: model levels, 2: pressure levels, 3: height levels
    INTEGER                       :: nvars = 0
    ! Metadata for missing value masking
    LOGICAL                    :: lmiss = .FALSE.         ! flag: true, if variables should be initialized with missval
    LOGICAL                    :: lmask_boundary =.FALSE. ! flag: true, if interpolation zone should be masked *in output*
  END TYPE t_var_list_intrinsic
  !
  TYPE t_var_list
    TYPE(t_var_list_intrinsic), POINTER :: p => NULL()
  END type t_var_list
  !
CONTAINS
  !-----------------------------------------------------------------------------
  ! remove all elements of a linked list
  ! check if all elements are removed
  SUBROUTINE delete_list(this_list)
    TYPE(t_var_list), INTENT(inout) :: this_list
    TYPE(t_list_element), POINTER       :: this, next

    next => this_list%p%first_list_element
    DO WHILE (ASSOCIATED(next))
      this => next
      next => this%next_list_element
      IF (this%field%info%allocated) THEN
        IF (ASSOCIATED(this%field%r_ptr)) THEN
          !$ACC EXIT DATA DELETE( this%field%r_ptr ) IF( this%field%info%lopenacc )
          DEALLOCATE (this%field%r_ptr)
        ELSE IF (ASSOCIATED(this%field%s_ptr)) THEN
          !$ACC EXIT DATA DELETE( this%field%s_ptr ) IF( this%field%info%lopenacc )
          DEALLOCATE (this%field%s_ptr)
        ELSE IF (ASSOCIATED(this%field%i_ptr)) THEN
          !$ACC EXIT DATA DELETE( this%field%i_ptr ) IF( this%field%info%lopenacc )
          DEALLOCATE (this%field%i_ptr)
        ELSE IF (ASSOCIATED(this%field%l_ptr)) THEN
          !$ACC EXIT DATA DELETE( this%field%l_ptr ) IF( this%field%info%lopenacc )
          DEALLOCATE (this%field%l_ptr)
        ENDIF
        this%field%info%allocated = .FALSE.
        CALL this%field%info%finalize()
        CALL this%field%info_dyn%finalize()
      ENDIF
      DEALLOCATE (this)
    END DO
    DEALLOCATE(this_list%p)
    NULLIFY(this_list%p)
  END SUBROUTINE delete_list
  !-----------------------------------------------------------------------------
  ! add a list element to the linked list
  SUBROUTINE append_list_element(this_list, new_list_element)
    TYPE(t_var_list),     INTENT(inout) :: this_list
    TYPE(t_list_element), POINTER       :: new_list_element
    TYPE(t_list_element), POINTER :: current_list_element
    !
    ! insert as first element if list is empty
    IF (.NOT. ASSOCIATED (this_list%p%first_list_element)) THEN
      ALLOCATE(this_list%p%first_list_element)
      new_list_element => this_list%p%first_list_element
      this_list%p%nvars = this_list%p%nvars + 1
      RETURN
    ENDIF
    !
    ! loop over list elements to find position
    current_list_element => this_list%p%first_list_element
    DO WHILE (ASSOCIATED(current_list_element%next_list_element)) 
      current_list_element => current_list_element%next_list_element
    ENDDO
    !
    ! insert element
    ALLOCATE(new_list_element)
    new_list_element%next_list_element => current_list_element%next_list_element
    current_list_element%next_list_element => new_list_element
    this_list%p%list_elements = this_list%p%list_elements+1
    this_list%p%nvars = this_list%p%nvars + 1
  END SUBROUTINE append_list_element

END MODULE mo_linked_list
