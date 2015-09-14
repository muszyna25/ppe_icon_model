!>
!! Routines for handling proxy variables e.g. accumulation buffers
!!
MODULE mo_derived_variable_handling

  USE self_vector
  USE self_map
  USE self_assert

  USE mo_kind, ONLY: wp
  USE mo_model_domain, ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma
  USE mo_ocean_nml,           ONLY: n_zlev
  USE mo_var_metadata_types, ONLY: varname_len
  USE mo_impl_constants, ONLY: vname_len, success, max_char_length
  USE mo_name_list_output_types, ONLY: t_output_name_list
  USE mo_mpi, ONLY: my_process_is_stdio
  USE mo_var_list_element, ONLY: level_type_ml, level_type_pl, level_type_hl, level_type_il
  USE mo_name_list_output_config, ONLY: first_output_name_list
  USE mo_var_list, ONLY: nvar_lists, max_var_lists, var_lists, new_var_list,&
       total_number_of_variables, collect_group, get_var_timelevel,&
       get_var_name, default_var_list_settings, add_var
  USE mo_linked_list, ONLY: find_list_element, t_var_list, t_list_element
  USE mo_exception, ONLY: finish, message, message_text

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_derived_variable_handling'

  TYPE(map)       , SAVE :: meanMap
  TYPE(vector)    , SAVE :: meanVariables(10)
  TYPE(t_var_list)   :: mean_stream_list
  INTEGER, PARAMETER :: ntotal = 1024

  PUBLIC :: init_mean_stream
  PUBLIC :: finish_mean_stream
  PUBLIC :: collect_meanstream_variables
  PUBLIC :: mean_stream_list
  PUBLIC :: copy_var_to_list

  TYPE :: t_accumulation_pair
    TYPE(t_list_element), POINTER :: source, destination
  END TYPE t_accumulation_pair

  interface copy_var_to_list
    module procedure copy_var3d_to_list
    module procedure copy_var2d_to_list
    module procedure copy_var1d_to_list
  end interface copy_var_to_list

!!!  SUBROUTINE collect_target_variables()
!!!  END SUBROUTINE collect_target_variables
  
CONTAINS

  !>
  !! Print contents of a vector, just giving the name for t_list_elements
  !!
  !! Optional label is printed first, on a line by its own
  !!
  SUBROUTINE var_print(this, label)
    TYPE(vector) , INTENT(in) :: this
    CHARACTER(*), INTENT(in), OPTIONAL :: label
    
    TYPE(vector_iterator) :: my_iter
    CLASS(*), POINTER :: my_buffer
    INTEGER :: i

    IF (PRESENT(label)) PRINT *, label

    my_iter = this%each()
    DO WHILE(my_iter%next(my_buffer))
      SELECT TYPE(my_buffer)
      TYPE is (t_list_element)
        PRINT *,'t_list_element:varname:',         trim(my_buffer%field%info%name)
      TYPE is (t_accumulation_pair)
        PRINT *,'t_accumulation_pair:source     :',trim(my_buffer%source%field%info%name)
        PRINT *,'t_accumulation_pair:destination:',trim(my_buffer%destination%field%info%name)
      CLASS default
        PRINT *,' default class print  :'
        CALL class_print(my_buffer)
      END SELECT
    END DO
  END SUBROUTINE var_print

  !>
  !! Create a variable list
  !!
  SUBROUTINE init_mean_stream(patch_2d)
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2d
    
    CHARACTER(LEN=max_char_length) :: listname
    
    integer :: i
    
    meanMap = map()
    do i=1,size(meanVariables,1)
    meanVariables(i)= vector(debug=.true.)
    enddo

    WRITE(listname,'(a)')  'mean_stream_list'
    CALL new_var_list(mean_stream_list, listname, patch_id=patch_2d%id)
    CALL default_var_list_settings( mean_stream_list,lrestart=.FALSE.,loutput=.TRUE., model_type='oce' )
  END SUBROUTINE init_mean_stream

  !>
  !!
  !!
  SUBROUTINE finish_mean_stream()
    CALL print_green(&
         '===================================================================')
    CALL print_green('FINISH MAP:')
    CALL meanMap%PRINT()
!!!    CALL print_green('FINISH VECTOR:')
!!!    CALL meanVariables%PRINT()
    CALL print_green(&
         '===================================================================')
!!!    CALL print_green('FINISH BUFFERS:')
!!!    PRINT *,varlist_buffer
!!!    PRINT *,periods_buffer
!!!    CALL print_green(&
!!!         '===================================================================')
  END SUBROUTINE finish_mean_stream

  !>
  !!
  !!
  SUBROUTINE collect_meanstream_variables(src_varlist1, src_varlist2,patch)
    TYPE(t_var_list)   :: src_varlist1
    TYPE(t_var_list)   :: src_varlist2
    type(t_patch)      :: patch

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::collect_meanStream_variables"
    CHARACTER(LEN=VARNAME_LEN) :: varname, mean_varname, message_text
    INTEGER :: nvars, i_typ, ierrstat, i, ntotal_vars, j, varlist_length
    CHARACTER(LEN=vname_len), POINTER :: in_varlist(:)
    TYPE (t_output_name_list), POINTER :: p_onl
    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: varlist(:)
    TYPE(t_list_element), POINTER :: src_element, dest_element
    REAL(wp), POINTER :: ptr(:,:,:)  !< reference to field
    TYPE(vector) :: keys 
    integer :: inml 
    type(vector) :: vector_buffer, value_buffer
    class(*), pointer :: buf
    type(vector_iterator) :: iter
    integer :: shape3d(3)
    type(t_accumulation_pair) :: accumulation_pair

    nvars = 1
    ntotal_vars = total_number_of_variables()
    ! temporary variables needed for variable group parsing
    ALLOCATE(varlist(ntotal_vars), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! -- loop over all output namelists
    p_onl => first_output_name_list
    inml = 1
    
    DO
      IF (.NOT.ASSOCIATED(p_onl)) EXIT
      IF ("mean" .EQ. TRIM(p_onl%operation)) THEN
!!!        WRITE(message_text,'(3a)') 'outputInterval: ',TRIM(p_onl%output_interval(1))
!!!        CALL message('',message_text)

        DO i_typ = 1, 4
   
          IF (i_typ == level_type_ml) in_varlist => p_onl%ml_varlist
          IF (i_typ == level_type_pl) in_varlist => p_onl%pl_varlist
          IF (i_typ == level_type_hl) in_varlist => p_onl%hl_varlist
          IF (i_typ == level_type_il) in_varlist => p_onl%il_varlist
   
          varlist_length = SIZE(in_varlist)

          DO
            IF (in_varlist(nvars) == ' ') EXIT
            nvars = nvars + 1
          END DO
          nvars = nvars - 1

!!!          WRITE(message_text,FMT=*) 'nvars: ',nvars
!!!          CALL message('',message_text)
   
          IF (nvars > 0)  varlist(1:nvars) = in_varlist(1:nvars)
          varlist((nvars+1):ntotal_vars) = " "
          
          shape3d = (/nproma,n_zlev,patch%alloc_cell_blocks/)
          IF (i_typ == level_type_ml) THEN
          write (0,*)'INML:',inml
          IF ( meanMap%has_key(TRIM(p_onl%output_interval(1))) ) THEN
            CALL meanMap%get(TRIM(p_onl%output_interval(1)),vector_buffer)
            call meanVariables(inml)%add(vector_buffer)
          !ELSE
            !call meanVariables(inml)%clear()
            !meanVariables = vector()
          END IF
            DO i=1,nvars
              ! collect data variables only, there variables names like
              ! 'grid:clon' which should be excluded
              IF ( INDEX(varlist(i),':') < 1 ) THEN
!!!                j = (periods_counter-1)*nvars + i
     
                ! find existing variable
                src_element => find_list_element (src_varlist1, TRIM(varlist(i)))
                IF (.NOT. ASSOCIATED (src_element)) src_element => &
                     find_list_element (src_varlist2, TRIM(varlist(i)))
                IF (.NOT. ASSOCIATED (src_element)) CALL finish( "collect_meanStream_variables", "Variable not found!")
                ! add new variable, copy the meta-data from the existing variable

              ! copy the source variable to destination pointer
              CALL copy_var_to_list(mean_stream_list,get_accumulation_varname(varlist(i),p_onl),src_element, ptr)
              dest_element => find_list_element(mean_stream_list,get_accumulation_varname(varlist(i),p_onl))
              CALL print_green('var:'//TRIM(src_element%field%info%name)//'---')
              CALL meanVariables(inml)%add(src_element) ! source element comes first
              CALL meanVariables(inml)%add(dest_element)
              end if
            end do

            call meanMap%add(trim(p_onl%output_interval(1)),meanVariables(inml),copy=.true.)
          END IF
        END DO
      inml = inml + 1
      END IF
      p_onl => p_onl%next
    END DO
    IF ( my_process_is_stdio() ) THEN
      call print_aqua('collected map {{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{')
      keys = meanMap%get_keys()
      call keys%print()
      value_buffer = meanMap%get_values()
      iter = value_buffer%each()

        DO WHILE(iter%next(buf))
          SELECT TYPE(buf)
          type is (vector)
            call var_print(buf)
          CLASS default
          call class_print(buf)
        end select
        end do
      call print_aqua('}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}')
    END IF
    !
    !1. Collect variables uniq by the output interval
    !   this will allow collective events for all variables in this group
    !2. for each variable, store
    !      source pointer [got with find_element]
    !      target pointer: copy meta info the the source, but choose new name with '_m'
    !
    !
    !write(0,*)'varlist:',varlist
  END SUBROUTINE collect_meanStream_variables
  SUBROUTINE copy_var3d_to_list(list,name,source_element,ptr)
    TYPE(t_var_list) :: list
    CHARACTER(LEN=VARNAME_LEN) :: name
    TYPE(t_list_element),POINTER :: source_element
    REAL(wp), POINTER           :: ptr(:,:,:)
    CALL add_var( list, &
      & name, &
      & ptr, source_element%field%info%hgrid, source_element%field%info%vgrid, &
      & source_element%field%info%cf, source_element%field%info%grib2,  &
      & ldims=source_element%field%info%used_dimensions(1:source_element%field%info%ndims), &
      & post_op=source_element%field%info%post_op, &
      & loutput=.TRUE., lrestart=.FALSE., &
      & var_class=source_element%field%info%var_class )
  END SUBROUTINE copy_var3d_to_list
  SUBROUTINE copy_var2d_to_list(list,name,source_element,ptr)
    TYPE(t_var_list) :: list
    CHARACTER(LEN=VARNAME_LEN) :: name
    TYPE(t_list_element),POINTER :: source_element
    REAL(wp), POINTER           :: ptr(:,:)
    CALL add_var( list, &
      & name, &
      & ptr, source_element%field%info%hgrid, source_element%field%info%vgrid, &
      & source_element%field%info%cf, source_element%field%info%grib2,  &
      & ldims=source_element%field%info%used_dimensions(1:source_element%field%info%ndims), &
      & post_op=source_element%field%info%post_op, &
      & loutput=.TRUE., lrestart=.FALSE., &
      & var_class=source_element%field%info%var_class )
  END SUBROUTINE copy_var2d_to_list
  SUBROUTINE copy_var1d_to_list(list,name,source_element,ptr)
    TYPE(t_var_list) :: list
    CHARACTER(LEN=VARNAME_LEN) :: name
    TYPE(t_list_element),POINTER :: source_element
    REAL(wp), POINTER           :: ptr(:)
    CALL add_var( list, &
      & name, &
      & ptr, source_element%field%info%hgrid, source_element%field%info%vgrid, &
      & source_element%field%info%cf, source_element%field%info%grib2,  &
      & ldims=source_element%field%info%used_dimensions(1:source_element%field%info%ndims), &
      & post_op=source_element%field%info%post_op, &
      & loutput=.TRUE., lrestart=.FALSE., &
      & var_class=source_element%field%info%var_class )
  END SUBROUTINE copy_var1d_to_list
  FUNCTION get_accumulation_varname(varname,output_setup)
    CHARACTER(LEN=VARNAME_LEN)  :: varname
    type(t_output_name_list) :: output_setup

    CHARACTER(LEN=VARNAME_LEN)  :: get_accumulation_varname
    CHARACTER(LEN=1)            :: separator

    separator = '_'
    get_accumulation_varname = &
      &TRIM(varname)//separator//&
      &trim(output_setup%operation)//separator//&
      &TRIM(output_setup%output_interval(1))//separator//&
      &trim(output_setup%output_start(1))

  END FUNCTION get_accumulation_varname
END MODULE mo_derived_variable_handling
