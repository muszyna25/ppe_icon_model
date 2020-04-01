MODULE mo_ser_all

#ifdef SERIALIZE
  USE m_serialize, ONLY: fs_write_field, &
                         fs_read_field, &
                         fs_add_savepoint_metainfo, &
                         fs_create_savepoint
  USE utils_ppser, ONLY: ppser_set_mode, &
                         ppser_get_mode, &
                         ppser_savepoint, &
                         ppser_serializer, &
                         ppser_serializer_ref, &
                         ppser_intlength, &
                         ppser_reallength, &
                         ppser_realtype, &
                         ppser_zrperturb
#endif

  USE mo_impl_constants,     ONLY: REAL_T, SINGLE_T, INT_T, BOOL_T, VARNAME_LEN
  USE mo_kind,               ONLY: vp, wp, sp
  USE mo_exception,          ONLY: warning, message
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_dynamics_config,    ONLY: nnow, nnew, nnow_rcf, nnew_rcf
  USE mo_linked_list,        ONLY: t_var_list, t_list_element
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_model_domain,       ONLY: t_patch
  USE mo_var_list,           ONLY: get_var_list
  USE mo_run_config,         ONLY: iforcing
  USE mo_impl_constants,     ONLY: inwp
  USE mo_ser_nml,            ONLY: ser_output_diag, ser_latbc_data, ser_dynamics, ser_diffusion, ser_step_advection, &
                                   ser_physics, ser_nudging, ser_all_debug

  IMPLICIT NONE

  PUBLIC :: serialize_all

  PRIVATE
  INTEGER :: output_diag_cnt = 0
  INTEGER :: latbc_data_cnt = 0
  INTEGER :: dynamics_cnt = 0
  INTEGER :: diffusion_cnt = 0
  INTEGER :: step_advection_cnt = 0
  INTEGER :: physics_cnt = 0
  INTEGER :: nudging_cnt = 0
  INTEGER :: debug_cnt = 0

  CONTAINS

  SUBROUTINE char_to_hash(c, a)
      CHARACTER(len=*), INTENT(in) :: c
      INTEGER, INTENT(out) :: a
      INTEGER :: i
  
      INTEGER :: p = 31
      INTEGER :: m = 1e5 + 5
      INTEGER :: p_pow
  
      a = 0
      p_pow = 1
      DO i=1,len_trim(c)
          a = MOD(a + (ichar(c(i:i)) + 1) * p_pow, m)
          p_pow = MOD(p_pow * p, m)
      END DO
  
  END SUBROUTINE


  SUBROUTINE ser_var_list(name, lupdate_cpu, domain, substr, timelev)
    CHARACTER(len=*), INTENT(IN)            :: name    ! name of output var_list
    LOGICAL, INTENT(IN)                     :: lupdate_cpu
    INTEGER, INTENT(IN), OPTIONAL           :: domain  ! domain index to append
    CHARACTER(len=*), INTENT(IN), OPTIONAL  :: substr  ! String after domain, before timelev
    INTEGER, INTENT(IN), OPTIONAL           :: timelev ! timelev index to append

    TYPE(t_var_list),     POINTER :: list
    TYPE(t_var_metadata), POINTER :: info
    TYPE(t_list_element), POINTER :: element
    CHARACTER(len=VARNAME_LEN)    :: listname

    INTEGER :: dims(5)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: ser_name, listhashchar
    INTEGER :: listhash

#ifdef SERIALIZE
    listname = TRIM(name)
    ! Append domain index
    IF( PRESENT(domain) ) THEN
      WRITE(listname, '(a,i2.2)') TRIM(listname), domain
    END IF

    ! Append substr
    IF( PRESENT(substr) ) THEN
      listname = TRIM(listname)//TRIM(substr)
    END IF

    ! Append timelev index
    IF( PRESENT(timelev) ) THEN
      WRITE(listname, '(a,i2.2)') TRIM(listname), timelev
    END IF

    CALL get_var_list(list, listname)

    element => list%p%first_list_element
    for_all_list_elements: DO WHILE (ASSOCIATED(element))
      info    => element%field%info
      IF(info%lopenacc) THEN
        call char_to_hash(listname, listhash)
        write(ser_name, '(a,i5.5)') TRIM(info%name)//'_', listhash
        dims = info%used_dimensions
        !print*, "writing out field ", ser_name, " dims: ", dims, " ndims: ", info%ndims
        SELECT CASE(info%data_type)
        CASE (REAL_T)
          SELECT CASE ( ppser_get_mode() )
            CASE(0)
              !$ACC UPDATE HOST( element%field%r_ptr ) IF( lupdate_cpu )
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,1,1,1,1))
                CASE(2)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,1,1,1))
                CASE(3)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,:,1,1))
                CASE(4)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,:,:,1))
              END SELECT
            CASE(1)
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,1,1,1,1))
                CASE(2)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,1,1,1))
                CASE(3)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,:,1,1))
                CASE(4)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,:,:,1))
              END SELECT
              !$ACC UPDATE DEVICE( element%field%r_ptr )
            CASE(2)
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,1,1,1,1), ppser_zrperturb)
                CASE(2)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,1,1,1), ppser_zrperturb)
                CASE(3)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,:,1,1), ppser_zrperturb)
                CASE(4)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,:,:,1), ppser_zrperturb)
              END SELECT
              !$ACC UPDATE DEVICE( element%field%r_ptr )
          END SELECT
        CASE (SINGLE_T)
          SELECT CASE ( ppser_get_mode() )
            CASE(0)
              !$ACC UPDATE HOST( element%field%s_ptr ) IF( lupdate_cpu )
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,1,1,1,1))
                CASE(2)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,1,1,1))
                CASE(3)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,:,1,1))
                CASE(4)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,:,:,1))
              END SELECT
            CASE(1)
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,1,1,1,1))
                CASE(2)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,1,1,1))
                CASE(3)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,:,1,1))
                CASE(4)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,:,:,1))
              END SELECT
              !$ACC UPDATE DEVICE( element%field%s_ptr )
            CASE(2)
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,1,1,1,1), ppser_zrperturb)
                CASE(2)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,1,1,1), ppser_zrperturb)
                CASE(3)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,:,1,1), ppser_zrperturb)
                CASE(4)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,:,:,1), ppser_zrperturb)
              END SELECT
              !$ACC UPDATE DEVICE( element%field%s_ptr )
          END SELECT
        CASE (INT_T)
          SELECT CASE ( ppser_get_mode() )
            CASE(0)
              !$ACC UPDATE HOST( element%field%i_ptr ) IF( lupdate_cpu )
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,1,1,1,1))
                CASE(2)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,1,1,1))
                CASE(3)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,:,1,1))
                CASE(4)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,:,:,1))
              END SELECT
            CASE(1)
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,1,1,1,1))
                CASE(2)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,1,1,1))
                CASE(3)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,:,1,1))
                CASE(4)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,:,:,1))
              END SELECT
              !$ACC UPDATE DEVICE( element%field%i_ptr )
            CASE(2)
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,1,1,1,1), ppser_zrperturb)
                CASE(2)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,1,1,1), ppser_zrperturb)
                CASE(3)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,:,1,1), ppser_zrperturb)
                CASE(4)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,:,:,1), ppser_zrperturb)
              END SELECT
              !$ACC UPDATE DEVICE( element%field%i_ptr )
          END SELECT
        CASE (BOOL_T)
        END SELECT
      END IF ! info%name /= 'tracer'
      element => element%next_list_element
    END DO for_all_list_elements
#endif

  END SUBROUTINE ser_var_list

  SUBROUTINE serialize_all(nproma, jg, savepoint_base, is_input, opt_lupdate_cpu, opt_id)

    INTEGER, INTENT(IN) :: nproma, jg
    CHARACTER(LEN=*), INTENT(IN) :: savepoint_base
    LOGICAL, INTENT(IN) :: is_input

    LOGICAL, INTENT(IN), OPTIONAL :: opt_lupdate_cpu
    INTEGER, INTENT(IN), OPTIONAL :: opt_id

    LOGICAL :: lupdate_cpu
    INTEGER :: id, ser_mode

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    CHARACTER(LEN=LEN(savepoint_base)+4) :: savepoint_name
    LOGICAL :: do_serialization

    do_serialization = .TRUE.

   SELECT CASE(savepoint_base)
      CASE("output_diag")
         output_diag_cnt = output_diag_cnt + 1
         IF(output_diag_cnt > ser_output_diag*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE("latbc_data")
         latbc_data_cnt = latbc_data_cnt + 1
         IF(latbc_data_cnt > ser_latbc_data*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE("dynamics")
         dynamics_cnt = dynamics_cnt + 1
         IF(dynamics_cnt > ser_dynamics*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE("diffusion")
         diffusion_cnt = diffusion_cnt + 1
         IF(diffusion_cnt > ser_diffusion*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE("step_advection")
         step_advection_cnt = step_advection_cnt + 1
         IF(step_advection_cnt > ser_step_advection*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE("physics")
         physics_cnt = physics_cnt + 1
         IF(physics_cnt > ser_physics*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE("nudging")
         nudging_cnt = nudging_cnt + 1
         IF(nudging_cnt > ser_nudging*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE DEFAULT
         debug_cnt = debug_cnt + 1
         IF(debug_cnt > ser_all_debug) THEN
             do_serialization = .FALSE.
         ENDIF
   END SELECT

    IF(do_serialization) THEN

       IF(is_input) THEN
#ifdef SERIALIZE_CREATE_REFERENCE
          ser_mode = 0 
#elif SERIALIZE_READ_REFERENCE
          ser_mode = 1
#elif SERIALIZE_PERTURB_REFERENCE
          ser_mode = 2
#endif
          savepoint_name = savepoint_base//"-in"
       ELSE
          savepoint_name = savepoint_base//"-out"
          ser_mode = 0
       ENDIF


       ! only update CPU data if requested and for write mode
       IF(PRESENT(opt_lupdate_cpu) .AND. ser_mode==0) THEN
           lupdate_cpu = opt_lupdate_cpu
       ELSE
           lupdate_cpu = .FALSE.
       ENDIF

       IF(PRESENT(opt_id)) THEN
           id = opt_id
       ELSE
           id = 0
       ENDIF

#ifdef SERIALIZE
       CALL warning('SER:'//TRIM(savepoint_name),'Serialization is active!')
       
#if defined(_OPENACC)
       IF(lupdate_cpu) THEN
         CALL warning('GPU:'//TRIM(savepoint_name),'GPU HOST synchronization forced by serialization!')
       ENDIF
#endif
       
       CALL datetimeToString(time_config%tc_current_date, date)
       CALL init('icon')
       call fs_create_savepoint(TRIM(savepoint_name), ppser_savepoint)
       call fs_add_savepoint_metainfo(ppser_savepoint, 'nproma', nproma)
       call fs_add_savepoint_metainfo(ppser_savepoint, 'date', TRIM(date))
       call fs_add_savepoint_metainfo(ppser_savepoint, 'id', id)
       call ppser_set_mode(ser_mode) ! 0=write, 1=read, 2=readperturb

       IF(iforcing == inwp) THEN
           ! Serialize NWP fields
           CALL ser_var_list('prm_diag_of_domain_', lupdate_cpu, domain=jg)
           CALL ser_var_list('prm_tend_of_domain_', lupdate_cpu, domain=jg)
           CALL ser_var_list('lnd_prog_of_domain_', lupdate_cpu, domain=jg, substr='_and_timelev_', timelev=nnow(jg))
           CALL ser_var_list('lnd_prog_of_domain_', lupdate_cpu, domain=jg, substr='_and_timelev_', timelev=nnew(jg))
           CALL ser_var_list('lnd_diag_of_domain_', lupdate_cpu, domain=jg)
           CALL ser_var_list('wtr_prog_of_domain_', lupdate_cpu, domain=jg, substr='_and_timelev_', timelev=nnow(jg))
           CALL ser_var_list('wtr_prog_of_domain_', lupdate_cpu, domain=jg, substr='_and_timelev_', timelev=nnew(jg))
           CALL ser_var_list('ext_data_atm_D', lupdate_cpu, domain=jg)
       ENDIF

       ! Serialize dynamics fields
       CALL ser_var_list('nh_state_metrics_of_domain_', lupdate_cpu, domain=jg)
       CALL ser_var_list('nh_state_diag_of_domain_', lupdate_cpu, domain=jg)
       CALL ser_var_list('nh_state_prog_of_domain_', lupdate_cpu, domain=jg, substr='_and_timelev_', timelev=nnew(jg)) !p_prog
       IF(nnow_rcf(jg) /= nnew(jg)) THEN
         CALL ser_var_list('nh_state_prog_of_domain_', lupdate_cpu, domain=jg, substr='_and_timelev_', timelev=nnow_rcf(jg)) !p_prog_now_rcf
       ENDIF
       IF(nnew_rcf(jg) /= nnew(jg) .AND. nnew_rcf(jg) /= nnow_rcf(jg)) THEN
         CALL ser_var_list('nh_state_prog_of_domain_', lupdate_cpu, domain=jg, substr='_and_timelev_', timelev=nnew_rcf(jg)) !p_prog_rcf
       ENDIF
       
#endif

    ENDIF !do_serialization
    
  END SUBROUTINE serialize_all

END MODULE mo_ser_all

