MODULE mo_ser_all

#ifdef SERIALIZE
  USE m_serialize, ONLY: fs_write_field, &
                         fs_read_field, &
                         fs_add_savepoint_metainfo, &
                         fs_create_savepoint
  USE utils_ppser, ONLY: ppser_savepoint, &
                         ppser_serializer, &
                         ppser_serializer_ref, &
                         ppser_zrperturb
#endif

  USE mo_impl_constants,     ONLY: REAL_T, SINGLE_T, INT_T, BOOL_T, VNAME_LEN
  USE mo_kind,               ONLY: vp, wp, sp
  USE mo_exception,          ONLY: warning, message
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_dynamics_config,    ONLY: nnow, nnew, nnow_rcf, nnew_rcf
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_model_domain,       ONLY: t_patch
  USE mo_var_list,           ONLY: t_var_list_ptr
  USE mo_var,                ONLY: t_var
  USE mo_var_list_register,  ONLY: vlr_get
  USE mo_run_config,         ONLY: iforcing, ldass_lhn
  USE mo_impl_constants,     ONLY: inwp
  USE mo_ser_nml,            ONLY: ser_output_diag, ser_latbc_data, ser_dynamics, ser_diffusion, ser_step_advection, &
                                   ser_physics, ser_lhn, ser_nudging, ser_all_debug, ser_surface, &
                                   ser_microphysics, ser_convection, ser_cover, ser_radiation, &
                                   ser_radheat, ser_gwdrag, ser_nfail, ser_nreport
  USE mo_mpi,                ONLY: get_my_mpi_work_id

  IMPLICIT NONE

  PUBLIC :: serialize_all

  PRIVATE
  INTEGER :: output_diag_cnt = 0
  INTEGER :: latbc_data_cnt = 0
  INTEGER :: dynamics_cnt = 0
  INTEGER :: diffusion_cnt = 0
  INTEGER :: step_advection_cnt = 0
  INTEGER :: physics_cnt = 0
  INTEGER :: lhn_cnt = 0
  INTEGER :: nudging_cnt = 0
  INTEGER :: surface_cnt = 0
  INTEGER :: turbtrans_cnt = 0
  INTEGER :: turbdiff_cnt = 0
  INTEGER :: microphysics_cnt = 0
  INTEGER :: convection_cnt = 0
  INTEGER :: cover_cnt = 0
  INTEGER :: radiation_cnt = 0
  INTEGER :: radheat_cnt = 0
  INTEGER :: gwdrag_cnt = 0
  INTEGER :: debug_cnt = 0

  REAL(wp) :: eps_r = EPSILON(1._wp)
  REAL(sp) :: eps_s = EPSILON(1._sp)
  !$acc declare copyin(eps_r, eps_s)

  INTERFACE compare
    MODULE PROCEDURE compare_r_1d
    MODULE PROCEDURE compare_r_2d
    MODULE PROCEDURE compare_r_3d
    MODULE PROCEDURE compare_r_4d
    MODULE PROCEDURE compare_s_1d
    MODULE PROCEDURE compare_s_2d
    MODULE PROCEDURE compare_s_3d
    MODULE PROCEDURE compare_s_4d
    MODULE PROCEDURE compare_i_1d
    MODULE PROCEDURE compare_i_2d
    MODULE PROCEDURE compare_i_3d
    MODULE PROCEDURE compare_i_4d
  END INTERFACE compare

  INTERFACE is_close
    MODULE PROCEDURE is_close_r
    MODULE PROCEDURE is_close_s
    MODULE PROCEDURE is_close_i
  END INTERFACE is_close

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


  SUBROUTINE ser_var_list(name, abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain, substr, timelev)
    CHARACTER(len=*), INTENT(IN)            :: name    ! name of output var_list
    INTEGER, INTENT(IN)                     :: abs_threshold
    INTEGER, INTENT(IN)                     :: rel_threshold
    LOGICAL, INTENT(IN)                     :: lupdate_cpu
    INTEGER, INTENT(IN)                     :: ser_mode
    INTEGER, INTENT(IN), OPTIONAL           :: domain  ! domain index to append
    CHARACTER(len=*), INTENT(IN), OPTIONAL  :: substr  ! String after domain, before timelev
    INTEGER, INTENT(IN), OPTIONAL           :: timelev ! timelev index to append

    TYPE(t_var_list_ptr)          :: list
    TYPE(t_var_metadata), POINTER :: info
    TYPE(t_var),          POINTER :: element
    CHARACTER(len=VNAME_LEN)      :: listname

    INTEGER :: dims(5), ii

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

    CALL vlr_get(list, listname)



!DR Test
    IF ( ASSOCIATED(list%p) ) THEN
      for_all_list_elements: DO ii = 1, list%p%nvars
        element => list%p%vl(ii)%p
        info    => element%info
        call char_to_hash(listname, listhash)
        write(ser_name, '(a,i5.5)') TRIM(info%name)//'_', listhash
        dims = info%used_dimensions

        SELECT CASE(info%data_type)
        CASE (REAL_T)
          SELECT CASE ( ser_mode )
            CASE(0) ! write
              !$ACC UPDATE HOST( element%r_ptr ) IF( lupdate_cpu .and.  info%lopenacc )
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,1,1,1,1))
                CASE(2)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,:,1,1,1))
                CASE(3)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,:,:,1,1))
                CASE(4)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,:,:,:,1))
              END SELECT
            CASE(1) ! read
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,1,1,1,1))
                CASE(2)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,:,1,1,1))
                CASE(3)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,:,:,1,1))
                CASE(4)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,:,:,:,1))
              END SELECT
              !$ACC UPDATE DEVICE( element%r_ptr ) if ( info%lopenacc )
            CASE(2) ! read perturb
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,1,1,1,1), ppser_zrperturb)
                CASE(2)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,:,1,1,1), ppser_zrperturb)
                CASE(3)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,:,:,1,1), ppser_zrperturb)
                CASE(4)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,:,:,:,1), ppser_zrperturb)
              END SELECT
              !$ACC UPDATE DEVICE( element%r_ptr ) if ( info%lopenacc )
            CASE(3) ! compare
              !$ACC UPDATE HOST( element%r_ptr ) IF( lupdate_cpu .and.  info%lopenacc )
              SELECT CASE ( info%ndims )
                CASE(1)
                  call compare(ser_name, element%r_ptr(:,1,1,1,1), info%lopenacc, abs_threshold, rel_threshold)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,1,1,1,1))
                CASE(2)
                  call compare(ser_name, element%r_ptr(:,:,1,1,1), info%lopenacc, abs_threshold, rel_threshold)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,:,1,1,1))
                CASE(3)
                  call compare(ser_name, element%r_ptr(:,:,:,1,1), info%lopenacc, abs_threshold, rel_threshold)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,:,:,1,1))
                CASE(4)
                  call compare(ser_name, element%r_ptr(:,:,:,:,1), info%lopenacc, abs_threshold, rel_threshold)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%r_ptr(:,:,:,:,1))
              END SELECT
          END SELECT
        CASE (SINGLE_T)
          SELECT CASE ( ser_mode )
            CASE(0)
              !$ACC UPDATE HOST( element%s_ptr ) IF( lupdate_cpu .and.  info%lopenacc )
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,1,1,1,1))
                CASE(2)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,:,1,1,1))
                CASE(3)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,:,:,1,1))
                CASE(4)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,:,:,:,1))
              END SELECT
            CASE(1)
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,1,1,1,1))
                CASE(2)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,:,1,1,1))
                CASE(3)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,:,:,1,1))
                CASE(4)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,:,:,:,1))
              END SELECT
              !$ACC UPDATE DEVICE( element%s_ptr ) if ( info%lopenacc )
            CASE(2)
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,1,1,1,1), ppser_zrperturb)
                CASE(2)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,:,1,1,1), ppser_zrperturb)
                CASE(3)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,:,:,1,1), ppser_zrperturb)
                CASE(4)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,:,:,:,1), ppser_zrperturb)
              END SELECT
              !$ACC UPDATE DEVICE( element%s_ptr ) if ( info%lopenacc )
            CASE(3) ! compare
              !$ACC UPDATE HOST( element%s_ptr ) IF( lupdate_cpu .and.  info%lopenacc )
              SELECT CASE ( info%ndims )
                CASE(1)
                  call compare(ser_name, element%s_ptr(:,1,1,1,1), info%lopenacc, abs_threshold, rel_threshold)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,1,1,1,1))
                CASE(2)
                  call compare(ser_name, element%s_ptr(:,:,1,1,1), info%lopenacc, abs_threshold, rel_threshold)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,:,1,1,1))
                CASE(3)
                  call compare(ser_name, element%s_ptr(:,:,:,1,1), info%lopenacc, abs_threshold, rel_threshold)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,:,:,1,1))
                CASE(4)
                  call compare(ser_name, element%s_ptr(:,:,:,:,1), info%lopenacc, abs_threshold, rel_threshold)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%s_ptr(:,:,:,:,1))
              END SELECT
          END SELECT
        CASE (INT_T)
          SELECT CASE ( ser_mode )
            CASE(0)
              !$ACC UPDATE HOST( element%i_ptr ) IF( lupdate_cpu .and.  info%lopenacc )
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,1,1,1,1))
                CASE(2)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,:,1,1,1))
                CASE(3)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,:,:,1,1))
                CASE(4)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,:,:,:,1))
              END SELECT
            CASE(1)
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,1,1,1,1))
                CASE(2)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,:,1,1,1))
                CASE(3)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,:,:,1,1))
                CASE(4)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,:,:,:,1))
              END SELECT
              !$ACC UPDATE DEVICE( element%i_ptr ) if ( info%lopenacc )
            CASE(2)
              SELECT CASE ( info%ndims )
                CASE(1)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,1,1,1,1), ppser_zrperturb)
                CASE(2)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,:,1,1,1), ppser_zrperturb)
                CASE(3)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,:,:,1,1), ppser_zrperturb)
                CASE(4)
                  call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,:,:,:,1), ppser_zrperturb)
              END SELECT
              !$ACC UPDATE DEVICE( element%i_ptr ) if ( info%lopenacc )
            CASE(3) ! compare
              !$ACC UPDATE HOST( element%i_ptr ) IF( lupdate_cpu .and.  info%lopenacc )
              SELECT CASE ( info%ndims )
                CASE(1)
                  call compare(ser_name, element%i_ptr(:,1,1,1,1), info%lopenacc, abs_threshold, rel_threshold)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,1,1,1,1))
                CASE(2)
                  call compare(ser_name, element%i_ptr(:,:,1,1,1), info%lopenacc, abs_threshold, rel_threshold)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,:,1,1,1))
                CASE(3)
                  call compare(ser_name, element%i_ptr(:,:,:,1,1), info%lopenacc, abs_threshold, rel_threshold)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,:,:,1,1))
                CASE(4)
                  call compare(ser_name, element%i_ptr(:,:,:,:,1), info%lopenacc, abs_threshold, rel_threshold)
                  call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%i_ptr(:,:,:,:,1))
              END SELECT
          END SELECT
        CASE (BOOL_T)
        END SELECT
        ! END IF

      ENDDO for_all_list_elements
    END IF
!DR End Test

!!$    element => list%p%first_list_element
!!$    for_all_list_elements: DO WHILE (ASSOCIATED(element))
!!$      info    => element%field%info
!!$      call char_to_hash(listname, listhash)
!!$      write(ser_name, '(a,i5.5)') TRIM(info%name)//'_', listhash
!!$      dims = info%used_dimensions
!!$      ! IF(info%lopenacc) THEN
!!$      SELECT CASE(info%data_type)
!!$      CASE (REAL_T)
!!$        SELECT CASE ( ser_mode )
!!$          CASE(0) ! write
!!$            !$ACC UPDATE HOST( element%field%r_ptr ) IF( lupdate_cpu .and.  info%lopenacc )
!!$            SELECT CASE ( info%ndims )
!!$              CASE(1)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,1,1,1,1))
!!$              CASE(2)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,1,1,1))
!!$              CASE(3)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,:,1,1))
!!$              CASE(4)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,:,:,1))
!!$            END SELECT
!!$          CASE(1) ! read
!!$            SELECT CASE ( info%ndims )
!!$              CASE(1)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,1,1,1,1))
!!$              CASE(2)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,1,1,1))
!!$              CASE(3)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,:,1,1))
!!$              CASE(4)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,:,:,1))
!!$            END SELECT
!!$            !$ACC UPDATE DEVICE( element%field%r_ptr ) if ( info%lopenacc )
!!$          CASE(2) ! read perturb
!!$            SELECT CASE ( info%ndims )
!!$              CASE(1)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,1,1,1,1), ppser_zrperturb)
!!$              CASE(2)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,1,1,1), ppser_zrperturb)
!!$              CASE(3)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,:,1,1), ppser_zrperturb)
!!$              CASE(4)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,:,:,1), ppser_zrperturb)
!!$            END SELECT
!!$            !$ACC UPDATE DEVICE( element%field%r_ptr ) if ( info%lopenacc )
!!$          CASE(3) ! compare
!!$            !$ACC UPDATE HOST( element%field%r_ptr ) IF( lupdate_cpu .and.  info%lopenacc )
!!$            SELECT CASE ( info%ndims )
!!$              CASE(1)
!!$                call compare(ser_name, element%field%r_ptr(:,1,1,1,1), info%lopenacc, abs_threshold, rel_threshold)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,1,1,1,1))
!!$              CASE(2)
!!$                call compare(ser_name, element%field%r_ptr(:,:,1,1,1), info%lopenacc, abs_threshold, rel_threshold)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,1,1,1))
!!$              CASE(3)
!!$                call compare(ser_name, element%field%r_ptr(:,:,:,1,1), info%lopenacc, abs_threshold, rel_threshold)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,:,1,1))
!!$              CASE(4)
!!$                call compare(ser_name, element%field%r_ptr(:,:,:,:,1), info%lopenacc, abs_threshold, rel_threshold)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%r_ptr(:,:,:,:,1))
!!$            END SELECT
!!$        END SELECT
!!$      CASE (SINGLE_T)
!!$        SELECT CASE ( ser_mode )
!!$          CASE(0)
!!$            !$ACC UPDATE HOST( element%field%s_ptr ) IF( lupdate_cpu .and.  info%lopenacc )
!!$            SELECT CASE ( info%ndims )
!!$              CASE(1)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,1,1,1,1))
!!$              CASE(2)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,1,1,1))
!!$              CASE(3)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,:,1,1))
!!$              CASE(4)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,:,:,1))
!!$            END SELECT
!!$          CASE(1)
!!$            SELECT CASE ( info%ndims )
!!$              CASE(1)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,1,1,1,1))
!!$              CASE(2)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,1,1,1))
!!$              CASE(3)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,:,1,1))
!!$              CASE(4)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,:,:,1))
!!$            END SELECT
!!$            !$ACC UPDATE DEVICE( element%field%s_ptr ) if ( info%lopenacc )
!!$          CASE(2)
!!$            SELECT CASE ( info%ndims )
!!$              CASE(1)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,1,1,1,1), ppser_zrperturb)
!!$              CASE(2)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,1,1,1), ppser_zrperturb)
!!$              CASE(3)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,:,1,1), ppser_zrperturb)
!!$              CASE(4)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,:,:,1), ppser_zrperturb)
!!$            END SELECT
!!$            !$ACC UPDATE DEVICE( element%field%s_ptr ) if ( info%lopenacc )
!!$          CASE(3) ! compare
!!$            !$ACC UPDATE HOST( element%field%s_ptr ) IF( lupdate_cpu .and.  info%lopenacc )
!!$            SELECT CASE ( info%ndims )
!!$              CASE(1)
!!$                call compare(ser_name, element%field%s_ptr(:,1,1,1,1), info%lopenacc, abs_threshold, rel_threshold)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,1,1,1,1))
!!$              CASE(2)
!!$                call compare(ser_name, element%field%s_ptr(:,:,1,1,1), info%lopenacc, abs_threshold, rel_threshold)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,1,1,1))
!!$              CASE(3)
!!$                call compare(ser_name, element%field%s_ptr(:,:,:,1,1), info%lopenacc, abs_threshold, rel_threshold)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,:,1,1))
!!$              CASE(4)
!!$                call compare(ser_name, element%field%s_ptr(:,:,:,:,1), info%lopenacc, abs_threshold, rel_threshold)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%s_ptr(:,:,:,:,1))
!!$            END SELECT
!!$        END SELECT
!!$      CASE (INT_T)
!!$        SELECT CASE ( ser_mode )
!!$          CASE(0)
!!$            !$ACC UPDATE HOST( element%field%i_ptr ) IF( lupdate_cpu .and.  info%lopenacc )
!!$            SELECT CASE ( info%ndims )
!!$              CASE(1)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,1,1,1,1))
!!$              CASE(2)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,1,1,1))
!!$              CASE(3)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,:,1,1))
!!$              CASE(4)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,:,:,1))
!!$            END SELECT
!!$          CASE(1)
!!$            SELECT CASE ( info%ndims )
!!$              CASE(1)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,1,1,1,1))
!!$              CASE(2)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,1,1,1))
!!$              CASE(3)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,:,1,1))
!!$              CASE(4)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,:,:,1))
!!$            END SELECT
!!$            !$ACC UPDATE DEVICE( element%field%i_ptr ) if ( info%lopenacc )
!!$          CASE(2)
!!$            SELECT CASE ( info%ndims )
!!$              CASE(1)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,1,1,1,1), ppser_zrperturb)
!!$              CASE(2)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,1,1,1), ppser_zrperturb)
!!$              CASE(3)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,:,1,1), ppser_zrperturb)
!!$              CASE(4)
!!$                call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,:,:,1), ppser_zrperturb)
!!$            END SELECT
!!$            !$ACC UPDATE DEVICE( element%field%i_ptr ) if ( info%lopenacc )
!!$          CASE(3) ! compare
!!$            !$ACC UPDATE HOST( element%field%i_ptr ) IF( lupdate_cpu .and.  info%lopenacc )
!!$            SELECT CASE ( info%ndims )
!!$              CASE(1)
!!$                call compare(ser_name, element%field%i_ptr(:,1,1,1,1), info%lopenacc, abs_threshold, rel_threshold)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,1,1,1,1))
!!$              CASE(2)
!!$                call compare(ser_name, element%field%i_ptr(:,:,1,1,1), info%lopenacc, abs_threshold, rel_threshold)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,1,1,1))
!!$              CASE(3)
!!$                call compare(ser_name, element%field%i_ptr(:,:,:,1,1), info%lopenacc, abs_threshold, rel_threshold)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,:,1,1))
!!$              CASE(4)
!!$                call compare(ser_name, element%field%i_ptr(:,:,:,:,1), info%lopenacc, abs_threshold, rel_threshold)
!!$                call fs_write_field(ppser_serializer, ppser_savepoint, TRIM(ser_name), element%field%i_ptr(:,:,:,:,1))
!!$            END SELECT
!!$        END SELECT
!!$      CASE (BOOL_T)
!!$      END SELECT
!!$      ! END IF
!!$      element => element%next_list_element
!!$    END DO for_all_list_elements
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
    CHARACTER(LEN=100) :: compare_file_name
    LOGICAL :: do_serialization
    INTEGER :: rel_threshold, abs_threshold

    do_serialization = .TRUE.

  ! make sure CPU and GPU are synchronous
  !$acc wait

   SELECT CASE(savepoint_base)
      CASE("output_diag")
         output_diag_cnt = output_diag_cnt + 1
         rel_threshold = ser_output_diag(2)
         abs_threshold = ser_output_diag(3)
         IF(output_diag_cnt > ser_output_diag(1)*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE("latbc_data")
         latbc_data_cnt = latbc_data_cnt + 1
         rel_threshold = ser_latbc_data(2)
         abs_threshold = ser_latbc_data(3)
         IF(latbc_data_cnt > ser_latbc_data(1)*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE("dynamics")
         dynamics_cnt = dynamics_cnt + 1
         rel_threshold = ser_dynamics(2)
         abs_threshold = ser_dynamics(3)
         IF(dynamics_cnt > ser_dynamics(1)*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE("diffusion")
         diffusion_cnt = diffusion_cnt + 1
         rel_threshold = ser_diffusion(2)
         abs_threshold = ser_diffusion(3)
         IF(diffusion_cnt > ser_diffusion(1)*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE("step_advection")
         step_advection_cnt = step_advection_cnt + 1
         rel_threshold = ser_step_advection(2)
         abs_threshold = ser_step_advection(3)
         IF(step_advection_cnt > ser_step_advection(1)*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE("physics")
         physics_cnt = physics_cnt + 1
         rel_threshold = ser_physics(2)
         abs_threshold = ser_physics(3)
         IF(physics_cnt > ser_physics(1)*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE("lhn")
          lhn_cnt = lhn_cnt + 1
          rel_threshold = ser_lhn(2)
          abs_threshold = ser_lhn(3)
          IF(lhn_cnt > ser_lhn(1)*2) THEN
              do_serialization = .FALSE.
          ENDIF
      CASE("nudging")
         nudging_cnt = nudging_cnt + 1
         rel_threshold = ser_nudging(2)
         abs_threshold = ser_nudging(3)
         IF(nudging_cnt > ser_nudging(1)*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE("surface")
         surface_cnt = surface_cnt + 1
         rel_threshold = ser_surface(2)
         abs_threshold = ser_surface(3)
         IF(surface_cnt > ser_surface(1)*2) THEN
             do_serialization = .FALSE.
         ENDIF
      CASE("microphysics")
          microphysics_cnt = microphysics_cnt + 1
          rel_threshold = ser_microphysics(2)
          abs_threshold = ser_microphysics(3)
          IF(microphysics_cnt > ser_microphysics(1)*2) THEN
              do_serialization = .FALSE.
          ENDIF
      CASE("convection")
          convection_cnt = convection_cnt + 1
          rel_threshold = ser_convection(2)
          abs_threshold = ser_convection(3)
          IF(convection_cnt > ser_convection(1)*2) THEN
              do_serialization = .FALSE.
          ENDIF
      CASE("cover")
          cover_cnt = cover_cnt + 1
          rel_threshold = ser_cover(2)
          abs_threshold = ser_cover(3)
          IF(cover_cnt > ser_cover(1)*2) THEN
              do_serialization = .FALSE.
          ENDIF
      CASE("radiation")
          radiation_cnt = radiation_cnt + 1
          rel_threshold = ser_radiation(2)
          abs_threshold = ser_radiation(3)
          IF(radiation_cnt > ser_radiation(1)*2) THEN
              do_serialization = .FALSE.
          ENDIF
      CASE("radheat")
          radheat_cnt = radheat_cnt + 1
          rel_threshold = ser_radheat(2)
          abs_threshold = ser_radheat(3)
          IF(radheat_cnt > ser_radheat(1)*2) THEN
              do_serialization = .FALSE.
          ENDIF
      CASE("gwdrag")
          gwdrag_cnt = gwdrag_cnt + 1
          rel_threshold = ser_gwdrag(2)
          abs_threshold = ser_gwdrag(3)
          IF(gwdrag_cnt > ser_gwdrag(1)*2) THEN
              do_serialization = .FALSE.
          ENDIF
      CASE DEFAULT
         debug_cnt = debug_cnt + 1
         rel_threshold = ser_all_debug(2)
         abs_threshold = ser_all_debug(3)
         IF(debug_cnt > ser_all_debug(1)) THEN
             do_serialization = .FALSE.
         ENDIF
   END SELECT

    IF(do_serialization) THEN

       IF(is_input) THEN
#ifdef SERIALIZE_CREATE_REFERENCE
          ser_mode = 0 ! write
#elif SERIALIZE_READ_REFERENCE
          ser_mode = 1 ! read
#elif SERIALIZE_PERTURB_REFERENCE
          ser_mode = 2 ! read perturb
#endif
          savepoint_name = savepoint_base//"-in"
       ELSE
          savepoint_name = savepoint_base//"-out"
#ifdef SERIALIZE_CREATE_REFERENCE
          ser_mode = 0 ! write
#elif SERIALIZE_READ_REFERENCE
          ser_mode = 3 ! compare
#elif SERIALIZE_PERTURB_REFERENCE
          ser_mode = 0 ! write
#endif
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
       IF(ser_mode == 3) THEN
         WRITE(compare_file_name, "(A,A1,A,A1,I2.2,A5,I2.2,A8)") savepoint_name, "_", TRIM(date), "_", id, &
           "_rank", get_my_mpi_work_id(), "_sum.txt"
         OPEN( unit=123, file=compare_file_name, action="WRITE")
         WRITE(compare_file_name, "(A,A1,A,A1,I2.2,A5,I2.2,A4)") savepoint_name, "_", TRIM(date), "_", id, &
           "_rank", get_my_mpi_work_id(), ".txt"
         OPEN( unit=1234, file=compare_file_name, action="WRITE")
         WRITE(123, "(A, T40,A7, T50,A7, T60,A7, T70,A7, T80,A7)"), "field", "rel", "abs", "%", "nfail", "ntot"
       ELSE
      
       ENDIF

       CALL init('icon')
       call fs_create_savepoint(TRIM(savepoint_name), ppser_savepoint)
       call fs_add_savepoint_metainfo(ppser_savepoint, 'nproma', nproma)
       call fs_add_savepoint_metainfo(ppser_savepoint, 'date', TRIM(date))
       call fs_add_savepoint_metainfo(ppser_savepoint, 'id', id)

       IF(iforcing == inwp) THEN
           ! Serialize NWP fields
           CALL ser_var_list('prm_diag_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
           CALL ser_var_list('prm_tend_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
           CALL ser_var_list('lnd_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnow(jg))
           CALL ser_var_list('lnd_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnew(jg))
           CALL ser_var_list('lnd_diag_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
           CALL ser_var_list('wtr_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnow(jg))
           CALL ser_var_list('wtr_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnew(jg))
           CALL ser_var_list('ext_data_atm_D', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
       ENDIF

       IF(ldass_lhn) THEN
           ! Serialize radar data fields
           CALL ser_var_list('radar_data_ct_dom_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
           CALL ser_var_list('radar_data_td_dom_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
       ENDIF

       ! Serialize dynamics fields
       CALL ser_var_list('nh_state_metrics_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
       CALL ser_var_list('nh_state_diag_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
       CALL ser_var_list('nh_state_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnew(jg)) !p_prog
       IF(nnow_rcf(jg) /= nnew(jg)) THEN
         CALL ser_var_list('nh_state_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnow_rcf(jg)) !p_prog_now_rcf
       ENDIF
       IF(nnew_rcf(jg) /= nnew(jg) .AND. nnew_rcf(jg) /= nnow_rcf(jg)) THEN
         CALL ser_var_list('nh_state_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnew_rcf(jg)) !p_prog_rcf
       ENDIF

       CLOSE( unit=123 )
       
#endif

    ENDIF !do_serialization
    
  END SUBROUTINE serialize_all

  SUBROUTINE is_close_r(ref, cur, abs_threshold, rel_threshold, rel_diff, abs_diff, out)
    REAL(wp), INTENT(IN) :: ref, cur
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold
    REAL(wp), INTENT(OUT) :: rel_diff, abs_diff
    LOGICAL, INTENT(OUT) :: out

    REAL(wp) :: maxval, at, rt
    !$acc routine seq

    abs_diff = ABS(cur - ref)
    maxval = MAX(ABS(cur), ABS(ref))

    ! compute relative difference for report
    IF(maxval <  eps_r) THEN 
      rel_diff = 0
    ELSE 
      rel_diff = abs_diff / maxval
    END IF

    ! thresholds given as negative exponents
    rt = 10._wp**(-rel_threshold)
    at = 10._wp**(-abs_threshold)

    ! threshold computed as in python's math.isclose
    out = abs_diff <= MAX(rt*maxval, at)

  END SUBROUTINE is_close_r

  SUBROUTINE is_close_s(ref, cur, abs_threshold, rel_threshold, rel_diff, abs_diff, out)
    REAL(sp), INTENT(IN) :: ref, cur
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold
    REAL(sp), INTENT(OUT) :: rel_diff, abs_diff
    LOGICAL, INTENT(OUT) :: out

    REAL(sp) :: maxval, at, rt
    !$acc routine seq

    abs_diff = ABS(cur - ref)
    maxval = MAX(ABS(cur), ABS(ref))

    ! compute relative difference for report
    IF(maxval < eps_s) THEN 
      rel_diff = 0
    ELSE 
      rel_diff = abs_diff / maxval
    END IF

    ! thresholds given as negative exponents
    rt = 10._sp**(-rel_threshold)
    at = 10._sp**(-abs_threshold)

    ! threshold computed as in python's math.isclose
    out = abs_diff <= MAX(rt*maxval, at)

  END SUBROUTINE is_close_s

  SUBROUTINE is_close_i(ref, cur, abs_threshold, rel_threshold, rel_diff, abs_diff, out)
    INTEGER, INTENT(IN) :: ref, cur
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold
    INTEGER, INTENT(OUT) :: abs_diff
    REAL(sp), INTENT(OUT) :: rel_diff
    LOGICAL, INTENT(OUT) :: out

    INTEGER  :: maxval
    REAL(sp) :: at, rt
    !$acc routine seq

    abs_diff = ABS(cur - ref)
    maxval = MAX(ABS(cur), ABS(ref))

    ! compute relative difference for report
    IF(maxval ==  0) THEN 
      rel_diff = 0
    ELSE 
      rel_diff = REAL(abs_diff, sp) / REAL(maxval, sp)
    END IF

    ! thresholds given as negative exponents
    rt = 10._sp**(-rel_threshold)
    at = 10._sp**(-abs_threshold)

    ! threshold computed as in python's math.isclose
    out = abs_diff <= MAX(rt*maxval, at)

  END SUBROUTINE is_close_i

  SUBROUTINE report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, n_tot)
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(wp), DIMENSION(:), INTENT(IN) :: report_rel_diff, report_abs_diff, report_cur, report_ref
    CHARACTER(len=*), DIMENSION(:), INTENT(IN) :: report_idx
    INTEGER, INTENT(IN) :: n_fail, n_tot
    REAL(wp) :: q
    INTEGER :: z

    q = REAL(n_fail, wp) / REAL(n_tot, wp) * 100
    IF(q > ser_nfail) THEN 
      WRITE(123, "(A, T40,E7.1E2, T50,E7.1E2, T60,F7.3, T70,I7, T80,I7)"), TRIM(name), report_rel_diff(1), report_abs_diff(1), q, n_fail, n_tot
      WRITE(1234, "(A, A, I, A, I, A, F7.3, A)"), TRIM(name), ": ", n_fail, " out of ", n_tot, " elements are off (", q, " %)"
      WRITE(1234, "(T10,A, T30,A, T50,A, T70,A, T90,A)"), "rel diff", "abs diff", "current", "reference", "index"
      DO z=1,ser_nreport
        WRITE(1234, "(T10,E14.8E2, T30,E14.8E2, T50,E14.8E2, T70,E14.8E2, T90,A)") report_rel_diff(z), report_abs_diff(z), report_cur(z), report_ref(z), TRIM(report_idx(z))
      END DO
    ENDIF
  END SUBROUTINE report

  SUBROUTINE compare_r_1d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(wp), INTENT(IN) :: cur(:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(wp), DIMENSION(size(cur, 1)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(1)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, z

    n_fail = 0
    !$acc data create(ref, rel_diff, abs_diff) present(cur) if(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$acc update device(ref) if(lopenacc)
    !$acc parallel default(none) if(lopenacc)
    !$acc loop gang vector collapse(1) reduction(+: n_fail)
    DO i=1,size(cur, 1)
      call is_close(ref(i), cur(i), abs_threshold, rel_threshold, rel_diff(i), abs_diff(i), out)
      IF (.NOT. out) THEN 
        n_fail = n_fail + 1
      ENDIF
    END DO
    !$acc end parallel

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$acc update host(rel_diff, abs_diff) if(lopenacc)
      mask = .TRUE.
      !$acc kernels copyout(cur_cpy) if(lopenacc)
      cur_cpy = cur
      !$acc end kernels
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        WRITE(report_idx(z), "(A,I6,A)") "(", i, ")"
        report_abs_diff(z) = abs_diff(i)
        report_rel_diff(z) = rel_diff(i)
        report_cur(z) = cur_cpy(i)
        report_ref(z) = ref(i)
        mask(i) = .FALSE.
      END DO
    END IF

    !$acc end data

    ! call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_r_1d

  SUBROUTINE compare_r_2d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(wp), INTENT(IN) :: cur(:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(wp), DIMENSION(size(cur, 1), size(cur, 2)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(2)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, z

    n_fail = 0
    !$acc data create(ref, rel_diff, abs_diff) present(cur) if(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$acc update device(ref) if(lopenacc)
    !$acc parallel default(none) if(lopenacc)
    !$acc loop gang vector collapse(2) reduction(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        call is_close(ref(i,j), cur(i,j), abs_threshold, rel_threshold, rel_diff(i,j), abs_diff(i,j), out)
        IF (.NOT. out) THEN 
          n_fail = n_fail + 1
        ENDIF
      END DO
    END DO
    !$acc end parallel

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$acc update host(rel_diff, abs_diff) if(lopenacc)
      mask = .TRUE.
      !$acc kernels copyout(cur_cpy) if(lopenacc)
      cur_cpy = cur
      !$acc end kernels
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        WRITE(report_idx(z), "(A,I6,A,I6,A)") "(", i, ",", j, ")"
        report_abs_diff(z) = abs_diff(i,j)
        report_rel_diff(z) = rel_diff(i,j)
        report_cur(z) = cur_cpy(i,j)
        report_ref(z) = ref(i,j)
        mask(i,j) = .FALSE.
      END DO
    END IF

    !$acc end data

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_r_2d

  SUBROUTINE compare_r_3d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(wp), INTENT(IN) :: cur(:,:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(wp), DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2), size(cur, 3))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(3)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, k, z

    n_fail = 0
    !$acc data create(ref, rel_diff, abs_diff) present(cur) if(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$acc update device(ref) if(lopenacc)
    !$acc parallel default(none) if(lopenacc)
    !$acc loop gang vector collapse(3) reduction(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        DO k=1,size(cur, 3)
          call is_close(ref(i,j,k), cur(i,j,k), abs_threshold, rel_threshold, rel_diff(i,j,k), abs_diff(i,j,k), out)
          IF (.NOT. out) THEN 
            n_fail = n_fail + 1
          ENDIF
        END DO
      END DO
    END DO
    !$acc end parallel

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$acc update host(rel_diff, abs_diff) if(lopenacc)
      mask = .TRUE.
      !$acc kernels copyout(cur_cpy) if(lopenacc)
      cur_cpy = cur
      !$acc end kernels
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        k = idx(3)
        WRITE(report_idx(z), "(A,I6,A,I6,A,I6,A)") "(", i, ",", j, ",", k, ")"
        report_abs_diff(z) = abs_diff(i,j,k)
        report_rel_diff(z) = rel_diff(i,j,k)
        report_cur(z) = cur_cpy(i,j,k)
        report_ref(z) = ref(i,j,k)
        mask(i,j,k) = .FALSE.
      END DO
    END IF

    !$acc end data

    ! call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_r_3d

  SUBROUTINE compare_r_4d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(wp), INTENT(IN) :: cur(:,:,:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(wp), DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(4)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, k, l, z

    n_fail = 0
    !$acc data create(ref, rel_diff, abs_diff) present(cur) if(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$acc update device(ref) if(lopenacc)
    !$acc parallel default(none) if(lopenacc)
    !$acc loop gang vector collapse(4) reduction(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        DO k=1,size(cur, 3)
          DO l=1,size(cur, 4)
            call is_close(ref(i,j,k,l), cur(i,j,k,l), abs_threshold, rel_threshold, rel_diff(i,j,k,l), abs_diff(i,j,k,l), out)
            IF (.NOT. out) THEN 
              n_fail = n_fail + 1
            ENDIF
          END DO
        END DO
      END DO
    END DO
    !$acc end parallel

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$acc update host(rel_diff, abs_diff) if(lopenacc)
      mask = .TRUE.
      !$acc kernels copyout(cur_cpy) if(lopenacc)
      cur_cpy = cur
      !$acc end kernels
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        k = idx(3)
        l = idx(4)
        WRITE(report_idx(z), "(A,I6,A,I6,A,I6,A,I6,A)") "(", i, ",", j, ",", k, ",", l, ")"
        report_abs_diff(z) = abs_diff(i,j,k,l)
        report_rel_diff(z) = rel_diff(i,j,k,l)
        report_cur(z) = cur_cpy(i,j,k,l)
        report_ref(z) = ref(i,j,k,l)
        mask(i,j,k,l) = .FALSE.
      END DO
    END IF

    !$acc end data

    ! call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_r_4d

  SUBROUTINE compare_s_1d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(sp), INTENT(IN) :: cur(:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(sp), DIMENSION(size(cur, 1)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(1)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, z

    n_fail = 0
    !$acc data create(ref, rel_diff, abs_diff) present(cur) if(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$acc update device(ref) if(lopenacc)
    !$acc parallel default(none) if(lopenacc)
    !$acc loop gang vector collapse(1) reduction(+: n_fail)
    DO i=1,size(cur, 1)
      call is_close(ref(i), cur(i), abs_threshold, rel_threshold, rel_diff(i), abs_diff(i), out)
      IF (.NOT. out) THEN 
        n_fail = n_fail + 1
      ENDIF
    END DO
    !$acc end parallel

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$acc update host(rel_diff, abs_diff) if(lopenacc)
      mask = .TRUE.
      !$acc kernels copyout(cur_cpy) if(lopenacc)
      cur_cpy = cur
      !$acc end kernels
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        WRITE(report_idx(z), "(A,I6,A)") "(", i, ")"
        report_abs_diff(z) = abs_diff(i)
        report_rel_diff(z) = rel_diff(i)
        report_cur(z) = cur_cpy(i)
        report_ref(z) = ref(i)
        mask(i) = .FALSE.
      END DO
    END IF

    !$acc end data

    ! call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_s_1d

  SUBROUTINE compare_s_2d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(sp), INTENT(IN) :: cur(:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(sp), DIMENSION(size(cur, 1), size(cur, 2)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(2)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, z

    n_fail = 0
    !$acc data create(ref, rel_diff, abs_diff) present(cur) if(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$acc update device(ref) if(lopenacc)
    !$acc parallel default(none) if(lopenacc)
    !$acc loop gang vector collapse(2) reduction(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        call is_close(ref(i,j), cur(i,j), abs_threshold, rel_threshold, rel_diff(i,j), abs_diff(i,j), out)
        IF (.NOT. out) THEN 
          n_fail = n_fail + 1
        ENDIF
      END DO
    END DO
    !$acc end parallel

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$acc update host(rel_diff, abs_diff) if(lopenacc)
      mask = .TRUE.
      !$acc kernels copyout(cur_cpy) if(lopenacc)
      cur_cpy = cur
      !$acc end kernels
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        WRITE(report_idx(z), "(A,I6,A,I6,A)") "(", i, ",", j, ")"
        report_abs_diff(z) = abs_diff(i,j)
        report_rel_diff(z) = rel_diff(i,j)
        report_cur(z) = cur_cpy(i,j)
        report_ref(z) = ref(i,j)
        mask(i,j) = .FALSE.
      END DO
    END IF

    !$acc end data

    ! call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_s_2d

  SUBROUTINE compare_s_3d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(sp), INTENT(IN) :: cur(:,:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(sp), DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2), size(cur, 3))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(3)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, k, z

    n_fail = 0
    !$acc data create(ref, rel_diff, abs_diff) present(cur) if(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$acc update device(ref) if(lopenacc)
    !$acc parallel default(none) if(lopenacc)
    !$acc loop gang vector collapse(3) reduction(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        DO k=1,size(cur, 3)
          call is_close(ref(i,j,k), cur(i,j,k), abs_threshold, rel_threshold, rel_diff(i,j,k), abs_diff(i,j,k), out)
          IF (.NOT. out) THEN 
            n_fail = n_fail + 1
          ENDIF
        END DO
      END DO
    END DO
    !$acc end parallel

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$acc update host(rel_diff, abs_diff) if(lopenacc)
      mask = .TRUE.
      !$acc kernels copyout(cur_cpy) if(lopenacc)
      cur_cpy = cur
      !$acc end kernels
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        k = idx(3)
        WRITE(report_idx(z), "(A,I6,A,I6,A,I6,A)") "(", i, ",", j, ",", k, ")"
        report_abs_diff(z) = abs_diff(i,j,k)
        report_rel_diff(z) = rel_diff(i,j,k)
        report_cur(z) = cur_cpy(i,j,k)
        report_ref(z) = ref(i,j,k)
        mask(i,j,k) = .FALSE.
      END DO
    END IF

    !$acc end data

    ! call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_s_3d

  SUBROUTINE compare_s_4d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(sp), INTENT(IN) :: cur(:,:,:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(sp), DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(4)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, k, l, z

    n_fail = 0
    !$acc data create(ref, rel_diff, abs_diff) present(cur) if(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$acc update device(ref) if(lopenacc)
    !$acc parallel default(none) if(lopenacc)
    !$acc loop gang vector collapse(4) reduction(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        DO k=1,size(cur, 3)
          DO l=1,size(cur, 4)
            call is_close(ref(i,j,k,l), cur(i,j,k,l), abs_threshold, rel_threshold, rel_diff(i,j,k,l), abs_diff(i,j,k,l), out)
            IF (.NOT. out) THEN 
              n_fail = n_fail + 1
            ENDIF
          END DO
        END DO
      END DO
    END DO
    !$acc end parallel

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$acc update host(rel_diff, abs_diff) if(lopenacc)
      mask = .TRUE.
      !$acc kernels copyout(cur_cpy) if(lopenacc)
      cur_cpy = cur
      !$acc end kernels
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        k = idx(3)
        l = idx(4)
        WRITE(report_idx(z), "(A,I6,A,I6,A,I6,A,I6,A)") "(", i, ",", j, ",", k, ",", l, ")"
        report_abs_diff(z) = abs_diff(i,j,k,l)
        report_rel_diff(z) = rel_diff(i,j,k,l)
        report_cur(z) = cur_cpy(i,j,k,l)
        report_ref(z) = ref(i,j,k,l)
        mask(i,j,k,l) = .FALSE.
      END DO
    END IF

    !$acc end data

    ! call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_s_4d

  SUBROUTINE compare_i_1d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: cur(:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    INTEGER, DIMENSION(size(cur, 1)) :: ref, cur_cpy, abs_diff
    REAL(sp), DIMENSION(size(cur, 1)) :: rel_diff
    LOGICAL :: mask(size(cur, 1))
    INTEGER :: report_cur(ser_nreport), report_ref(ser_nreport), report_abs_diff(ser_nreport)
    REAL(wp) :: report_rel_diff(ser_nreport)
    INTEGER :: idx(1)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, z

    n_fail = 0
    !$acc data create(ref, rel_diff, abs_diff) present(cur) if(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$acc update device(ref) if(lopenacc)
    !$acc parallel default(none) if(lopenacc)
    !$acc loop gang vector collapse(1) reduction(+: n_fail)
    DO i=1,size(cur, 1)
      call is_close(ref(i), cur(i), abs_threshold, rel_threshold, rel_diff(i), abs_diff(i), out)
      IF (.NOT. out) THEN 
        n_fail = n_fail + 1
      ENDIF
    END DO
    !$acc end parallel

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$acc update host(rel_diff, abs_diff) if(lopenacc)
      mask = .TRUE.
      !$acc kernels copyout(cur_cpy) if(lopenacc)
      cur_cpy = cur
      !$acc end kernels
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        WRITE(report_idx(z), "(A,I6,A)") "(", i, ")"
        report_abs_diff(z) = abs_diff(i)
        report_rel_diff(z) = rel_diff(i)
        report_cur(z) = cur_cpy(i)
        report_ref(z) = ref(i)
        mask(i) = .FALSE.
      END DO
    END IF

    !$acc end data

    ! call report(name, report_rel_diff, REAL(report_abs_diff, wp), REAL(report_cur, wp), REAL(report_ref, wp), report_idx, n_fail, size(cur))

  END SUBROUTINE compare_i_1d

  SUBROUTINE compare_i_2d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: cur(:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    INTEGER, DIMENSION(size(cur, 1), size(cur, 2)) :: ref, cur_cpy, abs_diff
    REAL(sp), DIMENSION(size(cur, 1), size(cur, 2)) :: rel_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2))
    INTEGER :: report_cur(ser_nreport), report_ref(ser_nreport), report_abs_diff(ser_nreport)
    REAL(wp) :: report_rel_diff(ser_nreport)
    INTEGER :: idx(2)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, z

    n_fail = 0
    !$acc data create(ref, rel_diff, abs_diff) present(cur) if(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$acc update device(ref) if(lopenacc)
    !$acc parallel default(none) if(lopenacc)
    !$acc loop gang vector collapse(2) reduction(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        call is_close(ref(i,j), cur(i,j), abs_threshold, rel_threshold, rel_diff(i,j), abs_diff(i,j), out)
        IF (.NOT. out) THEN 
          n_fail = n_fail + 1
        ENDIF
      END DO
    END DO
    !$acc end parallel

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$acc update host(rel_diff, abs_diff) if(lopenacc)
      mask = .TRUE.
      !$acc kernels copyout(cur_cpy) if(lopenacc)
      cur_cpy = cur
      !$acc end kernels
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        WRITE(report_idx(z), "(A,I6,A,I6,A)") "(", i, ",", j, ")"
        report_abs_diff(z) = abs_diff(i,j)
        report_rel_diff(z) = rel_diff(i,j)
        report_cur(z) = cur_cpy(i,j)
        report_ref(z) = ref(i,j)
        mask(i,j) = .FALSE.
      END DO
    END IF

    !$acc end data

    ! call report(name, report_rel_diff, REAL(report_abs_diff, wp), REAL(report_cur, wp), REAL(report_ref, wp), report_idx, n_fail, size(cur))

  END SUBROUTINE compare_i_2d

  SUBROUTINE compare_i_3d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: cur(:,:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    INTEGER, DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3)) :: ref, cur_cpy, abs_diff
    REAL(sp), DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3)) :: rel_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2), size(cur, 3))
    INTEGER :: report_cur(ser_nreport), report_ref(ser_nreport), report_abs_diff(ser_nreport)
    REAL(wp) :: report_rel_diff(ser_nreport)
    INTEGER :: idx(3)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, k, z

    n_fail = 0
    !$acc data create(ref, rel_diff, abs_diff) present(cur) if(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$acc update device(ref) if(lopenacc)
    !$acc parallel default(none) if(lopenacc)
    !$acc loop gang vector collapse(3) reduction(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        DO k=1,size(cur, 3)
          call is_close(ref(i,j,k), cur(i,j,k), abs_threshold, rel_threshold, rel_diff(i,j,k), abs_diff(i,j,k), out)
          IF (.NOT. out) THEN 
            n_fail = n_fail + 1
          ENDIF
        END DO
      END DO
    END DO
    !$acc end parallel

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$acc update host(rel_diff, abs_diff) if(lopenacc)
      mask = .TRUE.
      !$acc kernels copyout(cur_cpy) if(lopenacc)
      cur_cpy = cur
      !$acc end kernels
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        k = idx(3)
        WRITE(report_idx(z), "(A,I6,A,I6,A,I6,A)") "(", i, ",", j, ",", k, ")"
        report_abs_diff(z) = abs_diff(i,j,k)
        report_rel_diff(z) = rel_diff(i,j,k)
        report_cur(z) = cur_cpy(i,j,k)
        report_ref(z) = ref(i,j,k)
        mask(i,j,k) = .FALSE.
      END DO
    END IF

    !$acc end data

    ! call report(name, report_rel_diff, REAL(report_abs_diff, wp), REAL(report_cur, wp), REAL(report_ref, wp), report_idx, n_fail, size(cur))

  END SUBROUTINE compare_i_3d

  SUBROUTINE compare_i_4d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: cur(:,:,:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    INTEGER, DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4)) :: ref, cur_cpy, abs_diff
    REAL(sp), DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4)) :: rel_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4))
    INTEGER :: report_cur(ser_nreport), report_ref(ser_nreport), report_abs_diff(ser_nreport)
    REAL(wp) :: report_rel_diff(ser_nreport)
    INTEGER :: idx(4)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, k, l, z

    n_fail = 0
    !$acc data create(ref, rel_diff, abs_diff) present(cur) if(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$acc update device(ref) if(lopenacc)
    !$acc parallel default(none) if(lopenacc)
    !$acc loop gang vector collapse(4) reduction(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        DO k=1,size(cur, 3)
          DO l=1,size(cur, 4)
            call is_close(ref(i,j,k,l), cur(i,j,k,l), abs_threshold, rel_threshold, rel_diff(i,j,k,l), abs_diff(i,j,k,l), out)
            IF (.NOT. out) THEN 
              n_fail = n_fail + 1
            ENDIF
          END DO
        END DO
      END DO
    END DO
    !$acc end parallel

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$acc update host(rel_diff, abs_diff) if(lopenacc)
      mask = .TRUE.
      !$acc kernels copyout(cur_cpy) if(lopenacc)
      cur_cpy = cur
      !$acc end kernels
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        k = idx(3)
        l = idx(4)
        WRITE(report_idx(z), "(A,I6,A,I6,A,I6,A,I6,A)") "(", i, ",", j, ",", k, ",", l, ")"
        report_abs_diff(z) = abs_diff(i,j,k,l)
        report_rel_diff(z) = rel_diff(i,j,k,l)
        report_cur(z) = cur_cpy(i,j,k,l)
        report_ref(z) = ref(i,j,k,l)
        mask(i,j,k,l) = .FALSE.
      END DO
    END IF

    !$acc end data

    ! call report(name, report_rel_diff, REAL(report_abs_diff, wp), REAL(report_cur, wp), REAL(report_ref, wp), report_idx, n_fail, size(cur))

  END SUBROUTINE compare_i_4d


END MODULE mo_ser_all

