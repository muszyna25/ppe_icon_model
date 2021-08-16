MODULE mo_ser_manually
  ! Serialize any scalar and array component in selected global derived type
  ! variables and that were not initialized with add_var or add_ref.
  !
  ! HOW TO:
  ! 1) Add a USE-ONLY statement for each derived type variable
  ! 2) Add ser_component() calls to ser_manually() for each component

#ifdef SERIALIZE
  USE mo_ser_common,  ONLY: ser_component, t_ser_options
  USE mo_kind,        ONLY: sp, dp

#ifdef _OPENACC
  USE openacc,               ONLY: acc_is_present
#endif

  ! Global variables with scalars that are going to be serialized
  USE mo_nonhydro_state,     ONLY: p_nh_state
  USE mo_prepadv_state,      ONLY: prep_adv

    IMPLICIT NONE

  PUBLIC :: ser_manually

  PRIVATE

  ! Test dataset
  TYPE :: t_test
    INTEGER :: int
    REAL(dp) :: double
    REAL(sp) :: float
    LOGICAL :: log_a, log_b
  END TYPE t_test
  TYPE(t_test), TARGET :: p_test
  LOGICAL :: lscalars_test_was_run = .FALSE.


  CONTAINS

  SUBROUTINE ser_manually(abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain)
    INTEGER, INTENT(IN) :: abs_threshold
    INTEGER, INTENT(IN) :: rel_threshold
    LOGICAL, INTENT(IN) :: lupdate_cpu
    INTEGER, INTENT(IN) :: ser_mode
    INTEGER, INTENT(IN) :: domain
    TYPE(t_ser_options) :: o

    o%abs_threshold = abs_threshold
    o%rel_threshold = rel_threshold
    o%lupdate_cpu = lupdate_cpu
    o%lopenacc = .TRUE. ! all variables tested here should be available on GPU
    o%ser_mode = ser_mode
    o%domain = domain

    ! Make sure that serialbox and openACC work with scalars
    CALL ser_test_scalars()

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! List all scalars that should be serialized here !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL ser_component(o, "p_nh_state(jg)%metrics%zd_listdim", p_nh_state(domain)%metrics%zd_listdim)
    CALL ser_component(o, "p_nh_state(jg)%metrics%pg_listdim", p_nh_state(domain)%metrics%pg_listdim)
    CALL ser_component(o, "p_nh_state(jg)%metrics%nudge_c_dim", p_nh_state(domain)%metrics%nudge_c_dim)
    CALL ser_component(o, "p_nh_state(jg)%metrics%nudge_e_dim", p_nh_state(domain)%metrics%nudge_e_dim)
    CALL ser_component(o, "p_nh_state(jg)%metrics%bdy_halo_c_dim", p_nh_state(domain)%metrics%bdy_halo_c_dim)
    CALL ser_component(o, "p_nh_state(jg)%metrics%bdy_mflx_e_dim", p_nh_state(domain)%metrics%bdy_mflx_e_dim)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! List all arrays that should be serialized here  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF(ALLOCATED(prep_adv)) THEN
#ifdef _OPENACC
      ! prep_adv can be tested only if variable is already present on device
      IF( acc_is_present(prep_adv(domain:domain)) .AND. &
        & acc_is_present(prep_adv(domain)%vn_traj) ) THEN
#else
      IF(.TRUE.) THEN
#endif
        CALL ser_component(o, "mass_flx_me", prep_adv(domain)%mass_flx_me)
        CALL ser_component(o, "mass_flx_ic", prep_adv(domain)%mass_flx_ic)
        CALL ser_component(o, "vn_traj",     prep_adv(domain)%vn_traj)
        CALL ser_component(o, "topflx_tra",  prep_adv(domain)%topflx_tra)
      ENDIF
    ENDIF

  END SUBROUTINE ser_manually

! supporting subroutines:

  SUBROUTINE ser_test_scalars()
    ! test openACC update mechanism for scalars
    !
    ! The main objective of this test is to make sure, that references of the
    ! scalar components in p_test are passed to ser_scalar_* such that the
    ! ACC UPDATE inside ser_scalar_* affects the components of p_test.
    TYPE(t_ser_options) :: o

    IF(lscalars_test_was_run) return

    o%abs_threshold = 8
    o%rel_threshold = 8
    o%lupdate_cpu = .TRUE.
    o%lopenacc = .TRUE.

    p_test%int = 23
    p_test%float = 3.14
    p_test%double = 3.14159
    p_test%log_a = .TRUE.
    p_test%log_b = .FALSE.

    !$ACC DATA COPYIN(p_test)

#ifdef SERIALIZE_READ_REFERENCE
    ! manipulate GPU
    !$ACC SERIAL DEFAULT(NONE) PRESENT(p_test)
    p_test%int = 42
    p_test%float = 2.71
    p_test%double = 2.71828
    p_test%log_a = .FALSE.
    p_test%log_b = .TRUE.
    !$ACC END SERIAL

    o%ser_mode = 3 ! compare

    ! These should report 5 differences (which will be captured by the __TEST__ mode)
    CALL ser_component(o, "__TEST__error__p_test%int",    p_test%int)
    CALL ser_component(o, "__TEST__error__p_test%float",  p_test%float)
    CALL ser_component(o, "__TEST__error__p_test%double", p_test%double)
    CALL ser_component(o, "__TEST__error__p_test%log_a",  p_test%log_a)
    CALL ser_component(o, "__TEST__error__p_test%log_b",  p_test%log_b)

    o%ser_mode = 1 ! read serialized values
#else
    o%ser_mode = 0 ! write
#endif
    CALL ser_component(o, "__TEST__error__p_test%int",    p_test%int)
    CALL ser_component(o, "__TEST__error__p_test%float",  p_test%float)
    CALL ser_component(o, "__TEST__error__p_test%double", p_test%double)
    CALL ser_component(o, "__TEST__error__p_test%log_a",  p_test%log_a)
    CALL ser_component(o, "__TEST__error__p_test%log_b",  p_test%log_b)

#ifdef SERIALIZE_READ_REFERENCE
    o%ser_mode = 3 ! compare values, that were read in, to reference
#endif

    ! These should report no differences
    CALL ser_component(o, "__TEST__noerror__p_test%int",    p_test%int)
    CALL ser_component(o, "__TEST__noerror__p_test%float",  p_test%float)
    CALL ser_component(o, "__TEST__noerror__p_test%double", p_test%double)
    CALL ser_component(o, "__TEST__noerror__p_test%log_a",  p_test%log_a)
    CALL ser_component(o, "__TEST__noerror__p_test%log_b",  p_test%log_b)

    !$ACC END DATA

#ifdef SERIALIZE_READ_REFERENCE
    ! Turn this test off after it was run in read-and-compare mode.
    ! This test is always run in write mode to make sure that the comparison
    ! find compare values for any save point.
    lscalars_test_was_run = .TRUE.
#endif

  END SUBROUTINE ser_test_scalars

#endif

END MODULE mo_ser_manually
