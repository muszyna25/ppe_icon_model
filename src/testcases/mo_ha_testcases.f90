!>
!! 
!! @par Revision History
!!  Modification by Hui Wan (MPI-M, 2007-09-17):
!!  - variable testtype moved to mo_hydro_contro to avoid
!!    recursive dependency.
!!  Modification by Hui Wan (MPI-M, 2008-04-04):
!!  - control varialbe testtype renamed ctest_name.
!!  - default height of the mountain in the gravity wave test
!!    reduced from 2000 m to 50 m.
!!  Modification by Almut Gassmann, MPI-M, (2008-04-24)
!!  - add mountain induced Rossby wave test for NCAR workshop
!!  Modification by Marco Giorgetta, MPI-M, (2009-03-28)
!!  - for hydrostatic testcases, the switches ldynamics,
!!    ltransport and lforcing are set as required, thus
!!    overriding the defaults or run_ctl settings
!!  Modification by Marco Giorgetta, MPI-M, (2009-03-31)
!!  - overriding of run_ctl parameters is removed, to keep a clean
!!    separation of controling the initial states of the test cases
!!    and of controling the processes in the tendency calculation.
!!    The switch lpure_advection is removed, as it can be replaced
!!    by ldynamics=.false., ltransport=.true. and itime_scheme=3
!!    (to avoid effects of the SI-correction on the circulation).
!!  Modification by Daniel Reinert, DWD, (2009-10-07)
!!  - added stationary vortex test case for tracer transport
!!  Modification by Daniel Reinert, DWD, (2010-03-26)
!!  - added new class of deformational flow testcases for tracer 
!!    transport
!! 
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!! 
!! 
MODULE mo_ha_testcases
!-------------------------------------------------------------------------  
!  
!    ProTeX FORTRAN source: Style 2  
!    modified for ICON project, DWD/MPI-M 2004                            
!  
!-------------------------------------------------------------------------  
!  
!   
! 

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message_text, message, finish
  USE mo_impl_constants,  ONLY: SUCCESS, MAX_CHAR_LENGTH, TRACER_ONLY, MAX_NTRACER
  USE mo_io_units,        ONLY: nnml, nnml_output
  USE mo_namelist,        ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_config,   ONLY: isRestart
  USE mo_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist, open_and_restore_namelist, close_tmpfile
  USE mo_model_domain,    ONLY: t_patch
  USE mo_ext_data_types,  ONLY: t_external_data
  USE mo_grid_config,     ONLY: n_dom
  USE mo_intp_data_strc,  ONLY: t_int_state
  USE mo_parallel_config, ONLY: nproma
  USE mo_run_config,      ONLY: ltransport, iqv
  USE mo_dynamics_config, ONLY: ltwotime, lshallow_water, lcoriolis, &
                                nnow,nold,nnew
  USE mo_ha_dyn_config,   ONLY: ha_dyn_config
       
  USE mo_icoham_dyn_types, ONLY: t_hydro_atm
  USE mo_ha_prog_util,     ONLY: copy_prog_state,                    &
                               & init_hydro_state_prog_isoRest,      &
                               & hydro_state_prog_add_random

  USE mo_rh_test,         ONLY: init_hydro_state_prog_rhtest,        &
                                rh_wavenum, rh_init_shift_deg
  USE mo_jw_test,         ONLY: init_hydro_state_prog_jwtest,        &
                                jw_uptb
  USE mo_mrw_test,        ONLY: init_hydro_state_prog_mrossby,       &
                                init_hydro_state_prog_mrossby2,      &
                                mountctr_lon_deg, mountctr_lat_deg,  &
                                mountctr_height, mount_half_width,   &
                                mount_u0
  USE mo_gw_test,         ONLY: init_hydro_state_prog_gwtest,        &
                                gw_brunt_vais, gw_u0,                &
                                gw_lon_deg, gw_lat_deg
  USE mo_pa_test,         ONLY: init_hydro_state_prog_patest
  USE mo_sv_test,         ONLY: init_hydro_state_prog_svtest
  USE mo_df_test,         ONLY: init_hydro_state_prog_dftest,        &
                              & df_distv_barycenter, df_cell_indices 
  USE mo_sw_test,         ONLY: init_will2_test, init_will3_test,    &
                                init_will5_test, init_will6_test,    &
                                init_usbr_test, init_swgw_test
  USE mo_mpi,             ONLY: my_process_is_stdio
  USE mo_ldf_init,        ONLY: ldf_init_prog_state
  
  IMPLICIT NONE
  PRIVATE 
  PUBLIC :: read_ha_testcase_namelist, init_testcase, rotate_axis_deg  !subroutines

  ! Namelist variables

  CHARACTER(len=MAX_CHAR_LENGTH),PUBLIC :: ctest_name   ! Test case specifier
  CHARACTER(len=MAX_CHAR_LENGTH),PUBLIC :: ape_sst_case ! SST for APE experiments
  
  REAL(wp) :: rotate_axis_deg
  INTEGER, PUBLIC :: ihs_init_type ! 0=isothermal, 1=JWw
  LOGICAL, PUBLIC :: lhs_vn_ptb
  REAL(wp),PUBLIC :: hs_vn_ptb_scale

  LOGICAL, PUBLIC :: lrh_linear_pres  !< if .TRUE., initialize relative humidity
                                      !< with a linear function of pressure
  REAL(wp),PUBLIC :: rh_at_1000hpa    !< relative humidity [0,1] at 1000 hPa.

  LOGICAL, PUBLIC :: linit_tracer_fv  !< finite volume initialization for tracer 
                                      !< fields if .TRUE.

  INTEGER, PUBLIC :: ildf_init_type   ! isothermal atmosphere at rest (=0, default),
                                      ! JWs zonal state (=1) 
  LOGICAL, PUBLIC :: ldf_symm         ! if .TRUE. local diabatic forcing symmetric
                                      ! about the equator. if .FALSE. the forcing 
                                      ! is placed at 30 N. 

  INTEGER, PUBLIC :: tracer_inidist_list(MAX_NTRACER) ! Initial distribution of nth tracer
                                                      ! Applicable to test cases
                                                      ! df_test, pa_test, jw_test

  NAMELIST/ha_testcase_nml/ ctest_name, rotate_axis_deg, ape_sst_case,  &
    &                    gw_brunt_vais, gw_u0, gw_lon_deg, gw_lat_deg,  &
    &                    rh_wavenum, rh_init_shift_deg,                 &
    &                    mountctr_lon_deg, mountctr_lat_deg,            &
    &                    mountctr_height, mount_half_width,             &
    &                    jw_uptb, mount_u0, ihs_init_type, lhs_vn_ptb,  &
    &                    hs_vn_ptb_scale, lrh_linear_pres,              &
    &                    rh_at_1000hpa, linit_tracer_fv, ldf_symm,      &
    &                    ildf_init_type, tracer_inidist_list



!-------------------------------------------------------------------------   

  CONTAINS

!-------------------------------------------------------------------------
!
!

  !>
  !!               Initialization of variables that determine.
  !! 
  !!               Initialization of variables that determine
  !!               which test case to run
  !! 
  !! @par Revision History
  !!  Initial version by Hui Wan (2007-07-20)
  !! 
  SUBROUTINE read_ha_testcase_namelist( filename )

    CHARACTER(LEN=*),INTENT(IN) :: filename
    INTEGER :: istat, funit

    !0!CHARACTER(len=*), PARAMETER :: &
    !0!  &      routine = 'mo_ha_testcases:read_ha_testcase_namelist'

    !------------------------------
    ! Default values
    !------------------------------
    rotate_axis_deg   = 0.0_wp

    ctest_name        = 'JWw'
  
    rh_wavenum        = 4
    rh_init_shift_deg = 0._wp

    gw_brunt_vais     = 0.01_wp
    gw_u0             = 0.0_wp
    gw_lon_deg        = 180._wp
    gw_lat_deg        = 0._wp

    mountctr_lon_deg  = 90._wp
    mountctr_lat_deg  = 30._wp
    mountctr_height   = 2000._wp
    mount_half_width  = 1500000.0_wp
    mount_u0          = 20.0_wp

    jw_uptb           = 1._wp
    lrh_linear_pres   = .FALSE.
    rh_at_1000hpa     = 0.75_wp

    ihs_init_type     = 1      ! JWs zonal state
    lhs_vn_ptb        = .TRUE. ! add random noise to normal wind
    hs_vn_ptb_scale   = 1._wp  ! magnitude of the random noise

    ape_sst_case      = 'sst1'

    linit_tracer_fv   = .TRUE. ! finite volume initialization for tracer

    ildf_init_type    = 0      ! isothermal atmosphere at rest
    ldf_symm          = .TRUE. ! forcing symmetric about the equator

    ! initial tracer distributions for test cases
    ! df_test, pa_test, jw_test
    tracer_inidist_list(:) = 1

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values in the previous integration.
    !----------------------------------------------------------------
    IF (isRestart()) THEN
      funit = open_and_restore_namelist('ha_testcase_nml')
      READ(funit,NML=ha_testcase_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('ha_testcase_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, ha_testcase_nml)
    END SELECT
    CALL close_nml

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=ha_testcase_nml)
      CALL store_and_close_namelist(funit, 'ha_testcase_nml')
    ENDIF
    ! Write the contents of the namelist to an ASCII file.
    ! Probably will be removed later.
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=ha_testcase_nml)

  END SUBROUTINE read_ha_testcase_namelist
  !-------------------------------------------------------------------------
  !>
  !!               Initialization routine for the testcases.
  !! 
  !!               Initialization routine for the testcases
  !!               of the 3D hydrostatic test cases.
  !! 
  !! @par Revision History
  !!  Hui Wan, MPI-M (2007-07-26)
  !! 
  SUBROUTINE init_testcase(pt_patch, pt_hydro_state, pt_int_state, ext_data)


  TYPE(t_patch), TARGET, INTENT(INOUT) :: pt_patch(n_dom)
  TYPE(t_int_state)                    :: pt_int_state(n_dom)
  TYPE(t_hydro_atm)                    :: pt_hydro_state(n_dom)
  TYPE(t_external_data), INTENT(INOUT) :: ext_data(n_dom) 

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = &
                     & '(mo_hydro_testcases) init_testcase:'

  INTEGER :: jg
  INTEGER :: nlev              !< number of full levels
  INTEGER :: ist     !< status variable
!-------------------------------------------------

  SELECT CASE (ctest_name)
  CASE('PA','SV','DF1','DF2','DF3','DF4')

    IF ((.NOT.ltwotime).OR.(ha_dyn_config%itime_scheme/=TRACER_ONLY)) THEN 

      WRITE(message_text,'(A,I3)') 'running the hydrostatic atm model in'//&
           & TRIM(ctest_name)//' mode requires ltwotime = .TRUE. and '   //&
           & 'ha_dyn_nml:itime_scheme =', TRACER_ONLY

      CALL finish(TRIM(routine),TRIM(message_text))

    END IF
  END SELECT

DO jg = 1,n_dom

  nlev = pt_patch(jg)%nlev

  IF (lshallow_water) THEN

    IF (TRIM(ctest_name)=='PA') THEN

      CALL init_hydro_state_prog_patest(pt_patch(jg),       &
           &           pt_hydro_state(jg)%prog(nnow(jg)),   &
           &           pt_hydro_state(jg)%diag,             &
           &           pt_int_state(jg), ext_data(jg),      &
           &           rotate_axis_deg, linit_tracer_fv,    &
           &           tracer_inidist_list)

      CALL init_hydro_state_prog_patest(pt_patch(jg),       &
           &           pt_hydro_state(jg)%prog(nnew(jg)),   &
           &           pt_hydro_state(jg)%diag,             &
           &           pt_int_state(jg), ext_data(jg),      &
           &           rotate_axis_deg, linit_tracer_fv,    &
           &           tracer_inidist_list)
    ENDIF

  ENDIF

  IF (.NOT. lshallow_water) THEN
 
     SELECT CASE (TRIM(ctest_name))

     CASE ('GW')
        !
        IF (.NOT.ltwotime) &
        CALL init_hydro_state_prog_gwtest( lcoriolis, pt_patch(jg),   &
             & ext_data(jg), pt_hydro_state(jg)%prog(nold(jg)))

        CALL init_hydro_state_prog_gwtest( lcoriolis, pt_patch(jg),   &
             & ext_data(jg), pt_hydro_state(jg)%prog(nnow(jg)))

     CASE ('MRW')
        !
        IF (.NOT.ltwotime) &
             & CALL init_hydro_state_prog_mrossby(pt_patch(jg),&
             & pt_hydro_state(jg)%prog(nold(jg)), ext_data(jg))
        CALL init_hydro_state_prog_mrossby(pt_patch(jg),       &
             & pt_hydro_state(jg)%prog(nnow(jg)), ext_data(jg))

     CASE ('MRW2')
        !
        IF (.NOT.ltwotime) &
             & CALL init_hydro_state_prog_mrossby2(pt_patch(jg),&
             & pt_hydro_state(jg)%prog(nold(jg)), ext_data(jg))
        CALL init_hydro_state_prog_mrossby2(pt_patch(jg),       &
             & pt_hydro_state(jg)%prog(nnow(jg)), ext_data(jg))

     CASE ('RH')
        !
        IF (.NOT.ltwotime) &
             & CALL init_hydro_state_prog_rhtest(pt_patch(jg), &
             & pt_hydro_state(jg)%prog(nold(jg)), ext_data(jg))
        CALL init_hydro_state_prog_rhtest(pt_patch(jg),        &
             & pt_hydro_state(jg)%prog(nnow(jg)), ext_data(jg))

     CASE ('JWs')
        !
        jw_uptb = 0.0_wp
        IF (.NOT.ltwotime) &
             & CALL init_hydro_state_prog_jwtest(pt_patch(jg),&
             & pt_hydro_state(jg)%prog(nold(jg)),             &
             & pt_hydro_state(jg)%diag, ext_data(jg),         &
             & rotate_axis_deg, tracer_inidist_list)
        CALL init_hydro_state_prog_jwtest(pt_patch(jg), &
             & pt_hydro_state(jg)%prog(nnow(jg)),       &
             & pt_hydro_state(jg)%diag, ext_data(jg),   &
             & rotate_axis_deg, tracer_inidist_list)

     CASE ('JWw')

        IF (.NOT.ltwotime) &
             & CALL init_hydro_state_prog_jwtest(pt_patch(jg),&
             & pt_hydro_state(jg)%prog(nold(jg)),             &
             & pt_hydro_state(jg)%diag, ext_data(jg),         &
             & rotate_axis_deg, tracer_inidist_list)
        CALL init_hydro_state_prog_jwtest(pt_patch(jg), &
             & pt_hydro_state(jg)%prog(nnow(jg)),       &
             & pt_hydro_state(jg)%diag,  ext_data(jg),  &
             & rotate_axis_deg, tracer_inidist_list)

     CASE ('PA')
        CALL init_hydro_state_prog_patest(pt_patch(jg),       &
             &           pt_hydro_state(jg)%prog(nnow(jg)),   &
             &           pt_hydro_state(jg)%diag,             &
             &           pt_int_state(jg), ext_data(jg),      &
             &           rotate_axis_deg, linit_tracer_fv,    &
             &           tracer_inidist_list)

        CALL init_hydro_state_prog_patest(pt_patch(jg),       &
             &           pt_hydro_state(jg)%prog(nnew(jg)),   &
             &           pt_hydro_state(jg)%diag,             &
             &           pt_int_state(jg), ext_data(jg),      &
             &           rotate_axis_deg, linit_tracer_fv,    &
             &           tracer_inidist_list)


     CASE ('SV')
        !
        CALL init_hydro_state_prog_svtest(pt_patch(jg),       &
             &           pt_hydro_state(jg)%prog(nnow(jg)),   &
             &           pt_hydro_state(jg)%diag,             &
             &           pt_int_state(jg), ext_data(jg))

        CALL init_hydro_state_prog_svtest(pt_patch(jg),       &
             &           pt_hydro_state(jg)%prog(nnew(jg)),   &
             &           pt_hydro_state(jg)%diag,             &
             &           pt_int_state(jg), ext_data(jg))
        !

     CASE ('DF1', 'DF2', 'DF3', 'DF4')
        !
        !
        ! allocate temporary arrays for distance vectors and upwind cells
        ALLOCATE( df_distv_barycenter(nproma,nlev,pt_patch(jg)%nblks_e,2),    &
          &       df_cell_indices(nproma,nlev,pt_patch(jg)%nblks_e,2),        &
          &       STAT=ist )
        IF (ist /= SUCCESS) THEN
          CALL finish ( TRIM(routine),                                     &
            &  'allocation for df_distv_barycenter, df_cell_indices,'  //  &
            &  'failed' )
         ENDIF

        CALL init_hydro_state_prog_dftest(pt_patch(jg),       &
             &           pt_hydro_state(jg)%prog(nnow(jg)),   &
             &           pt_hydro_state(jg)%diag,             &
             &           pt_int_state(jg), ext_data(jg),      &
             &           rotate_axis_deg, ctest_name,         &
             &           linit_tracer_fv, tracer_inidist_list )

        CALL init_hydro_state_prog_dftest(pt_patch(jg),       &
             &           pt_hydro_state(jg)%prog(nnew(jg)),   &
             &           pt_hydro_state(jg)%diag,             &
             &           pt_int_state(jg), ext_data(jg),      &
             &           rotate_axis_deg, ctest_name,         &
             &           linit_tracer_fv, tracer_inidist_list )
        !

     CASE ('HS')
        !
        ! First initialize the background state
        !
        SELECT CASE (ihs_init_type)
        CASE (1) !zonal state as defined in the JW steady state test
           !
           jw_uptb = 0.0_wp
           CALL init_hydro_state_prog_jwtest( pt_patch(jg), &
                & pt_hydro_state(jg)%prog(nnow(jg)),        &
                & pt_hydro_state(jg)%diag, ext_data(jg),    &
                & rotate_axis_deg, tracer_inidist_list)
           !
           CALL message(TRIM(routine),'Initial state used in &
                & the Held-Suarez test: JW steady')
           !
        CASE DEFAULT ! isothermal state at rest
           !
           CALL init_hydro_state_prog_isoRest( 300._wp, 100000._wp, &
                & pt_hydro_state(jg)%prog(nnow(jg)) )
           !
           CALL message(TRIM(routine),'Initial state used in &
                & the Held-Suarez test: isothermal state at rest')
           !
        END SELECT ! initial conditions of the HS test
        !
        ! Set the surface geopotential to zero
        !
        ext_data(jg)%atm%topography_c = 0._wp
        !
        ! (optional) Add random noise to the wind field
        !
        IF (lhs_vn_ptb) THEN
           CALL hydro_state_prog_add_random( pt_patch(jg), & ! input
                & pt_hydro_state(jg)%prog(nnow(jg)),     & ! in and out
                & hs_vn_ptb_scale, nproma, nlev )            ! input
           !
           CALL message(TRIM(routine),'Initial state used in the &
                & Held-Suarez test: random noised added to the normal wind')
        END IF
        !
        IF (.NOT.ltwotime) CALL copy_prog_state(    &
             & pt_hydro_state(jg)%prog(nnow(jg)), & !in
             & pt_hydro_state(jg)%prog(nold(jg)), & !out
             & .FALSE.,     &! copy temp rather than theta
             & ltransport )  ! copy tracer field if transport is on.

     CASE ('LDF')
        !
        ! First initialize the background state
        !
        SELECT CASE (ildf_init_type)
        CASE (1) !zonal state similar to JWs
           !
           CALL ldf_init_prog_state( pt_patch(jg),          &
                & pt_hydro_state(jg)%prog(nnow(jg)),        &
                & pt_hydro_state(jg)%diag, ext_data(jg),    &
                & rotate_axis_deg)
           !
           CALL message(TRIM(routine),'Initial state used in &
                & the local diabatic forcing test: balanced initial&
                & conditions similar to the JWs testcase')
           !
        CASE DEFAULT ! isothermal state at rest
           !
           CALL init_hydro_state_prog_isoRest( 300._wp, 100000._wp, &
                & pt_hydro_state(jg)%prog(nnow(jg)) )
           !
           CALL message(TRIM(routine),'Initial state used in &
                & the local diabatic forcing test: isothermal state at rest')
           !
           ! Set the surface geopotential to zero
           !
           ext_data(jg)%atm%topography_c = 0._wp
           !
        END SELECT ! initial conditions of the LDF test
        !
        !
        IF (.NOT.ltwotime) CALL copy_prog_state(    &
             & pt_hydro_state(jg)%prog(nnow(jg)), & !in
             & pt_hydro_state(jg)%prog(nold(jg)), & !out
             & .FALSE.,     &! copy temp rather than theta
             & ltransport )  ! copy tracer field if transport is on.

     CASE ('LDF-Moist')
        !
        ! First initialize the background state
        !
        SELECT CASE (ildf_init_type)
        CASE (1)
        ! Initialize the background state
        ! similar to the JWw-Moist background state
        !
         CALL ldf_init_prog_state( pt_patch(jg),          &
              & pt_hydro_state(jg)%prog(nnow(jg)),        &
              & pt_hydro_state(jg)%diag, ext_data(jg),    &
              & rotate_axis_deg,                          &
              & lrh_linear_pres, rh_at_1000hpa            )
        !
        CALL message(TRIM(routine),'Initial state used in &
                & the LDF-Moist test: balanced initial    &
                & conditions similar to the JWw-Moist testcase')
        !
        CASE DEFAULT ! isothermal state at rest
           !
           CALL init_hydro_state_prog_isoRest( 300._wp, 100000._wp, &
                & pt_hydro_state(jg)%prog(nnow(jg)) )
           !
           CALL message(TRIM(routine),'Initial state used in &
                & the local diabatic forcing test: isothermal state at rest')
           !
           ! Set the surface geopotential to zero
           !
           ext_data(jg)%atm%topography_c = 0._wp
           !
        END SELECT ! initial conditions of the LDF test
        !
        !
        IF (.NOT.ltwotime) CALL copy_prog_state(    &
             & pt_hydro_state(jg)%prog(nnow(jg)), & !in
             & pt_hydro_state(jg)%prog(nold(jg)), & !out
             & .FALSE.,     &! copy temp rather than theta
             & ltransport )  ! copy tracer field if transport is on.

     CASE ('JWw-Moist')

        CALL init_hydro_state_prog_jwtest(pt_patch(jg), &
             & pt_hydro_state(jg)%prog(nnow(jg)),       &
             & pt_hydro_state(jg)%diag,  ext_data(jg),  &
             & rotate_axis_deg, tracer_inidist_list,    &
             & lrh_linear_pres, rh_at_1000hpa )

        IF (.NOT.ltwotime) CALL copy_prog_state(  &
             & pt_hydro_state(jg)%prog(nnow(jg)), & !in
             & pt_hydro_state(jg)%prog(nold(jg)), & !out
             & .FALSE.,     &! copy temp rather than theta
             & ltransport )  ! copy tracer field if transport is on.

     CASE ('APE','APEi','APEc','RCEhydro')
        ! Initial conditions are the same as for the 'JWw-Moist' case

        SELECT CASE (ihs_init_type)
        CASE (0) ! isothermal state at rest
           !
           CALL init_hydro_state_prog_isoRest( 300._wp, 100000._wp, &
                & pt_hydro_state(jg)%prog(nnow(jg)) )
           !
           CALL message(TRIM(routine),'Initial state used in &
                & APE test: isothermal state at rest')
           !
        CASE DEFAULT 
           !
          CALL init_hydro_state_prog_jwtest(pt_patch(jg), &
              & pt_hydro_state(jg)%prog(nnow(jg)),       &
              & pt_hydro_state(jg)%diag,  ext_data(jg),  &
              & rotate_axis_deg, tracer_inidist_list,    &
              & lrh_linear_pres, rh_at_1000hpa )
           !
           CALL message(TRIM(routine),'Initial state used in &
                & the APE test: JW steady')
        
        END SELECT ! initial conditions of the HS test

        pt_hydro_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqv+1:) = 0._wp

        IF (.NOT.ltwotime) CALL copy_prog_state(  &
             & pt_hydro_state(jg)%prog(nnow(jg)), & !in
             & pt_hydro_state(jg)%prog(nold(jg)), & !out
             & .FALSE.,     &! copy temp rather than theta
             & ltransport )  ! copy tracer field if transport is on.

        ! Surface geopotential height is constantly zero.

        ext_data(jg)%atm%topography_c = 0._wp

     CASE ('AMIP')
        ! Initial conditions are the same as for the 'JWw-Moist' case

        !
          CALL init_hydro_state_prog_jwtest(pt_patch(jg), &
              & pt_hydro_state(jg)%prog(nnow(jg)),       &
              & pt_hydro_state(jg)%diag,  ext_data(jg),  &
              & rotate_axis_deg, tracer_inidist_list,    &
              & lrh_linear_pres, rh_at_1000hpa )
        !
           CALL message(TRIM(routine),'Initial state used in &
                & the AMIP test: JW steady')
        
        pt_hydro_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqv+1:) = 0._wp

        ! Surface geopotential height is constantly zero.
        
        ! For the time being, we start with a flat topography which is gradually "grown"
        ! to the external field ext_data%atm%elevation_c, during about 2 years. Set initial
        ! surface pressure to 98264 Pa so that a realistic atmospheric mass is reached after this time.
        pt_hydro_state(jg)%prog(nnow(jg))%pres_sfc(:,:) = 98264._wp
        CALL message(TRIM(routine), 'Initial state for surface pressure set to 98264 Pa')

        ext_data(jg)%atm%topography_c(:,:) = 0._wp

        IF (.NOT.ltwotime) CALL copy_prog_state(  &
             & pt_hydro_state(jg)%prog(nnow(jg)), & !in
             & pt_hydro_state(jg)%prog(nold(jg)), & !out
             & .FALSE.,     &! copy temp rather than theta
             & ltransport )  ! copy tracer field if transport is on.

     CASE default
       CALL finish(TRIM(routine),'unknown choice of TESTTYPE')
     END SELECT

  ELSE ! shallow_water

     SELECT CASE (TRIM(ctest_name))

     CASE ('Will_2')
     IF (.NOT. ltwotime) &
     CALL init_will2_test( pt_patch(jg),pt_hydro_state(jg)%prog(nold(jg)),&
                           ext_data(jg), rotate_axis_deg)
     CALL init_will2_test( pt_patch(jg),pt_hydro_state(jg)%prog(nnow(jg)),&
                           ext_data(jg), rotate_axis_deg)

     CASE ('Will_3')
     IF (.NOT. ltwotime) &
     CALL init_will3_test( pt_patch(jg),pt_hydro_state(jg)%prog(nold(jg)),&
                           ext_data(jg), rotate_axis_deg)
     CALL init_will3_test( pt_patch(jg),pt_hydro_state(jg)%prog(nnow(jg)),&
                           ext_data(jg), rotate_axis_deg)

     CASE ('Will_5')
     IF (.NOT. ltwotime) &
     CALL init_will5_test( pt_patch(jg),pt_hydro_state(jg)%prog(nold(jg)), &
       &                   ext_data(jg))
     CALL init_will5_test( pt_patch(jg),pt_hydro_state(jg)%prog(nnow(jg)), &
       &                   ext_data(jg))

     CASE ('Will_6')
     IF (.NOT. ltwotime) &
     CALL init_will6_test( pt_patch(jg),pt_hydro_state(jg)%prog(nold(jg)), &
       &                   ext_data(jg))
     CALL init_will6_test( pt_patch(jg),pt_hydro_state(jg)%prog(nnow(jg)), &
       &                   ext_data(jg))

     CASE ('USBR')
     IF (.NOT. ltwotime) &
     CALL init_usbr_test( pt_patch(jg),pt_hydro_state(jg)%prog(nold(jg)), &
       &                  ext_data(jg))
     CALL init_usbr_test( pt_patch(jg),pt_hydro_state(jg)%prog(nnow(jg)), &
       &                  ext_data(jg))

     CASE ('SW_GW')
     IF (.NOT. ltwotime) &
     CALL init_swgw_test( pt_patch(jg),pt_hydro_state(jg)%prog(nold(jg)), &
       &                  ext_data(jg))
     CALL init_swgw_test( pt_patch(jg),pt_hydro_state(jg)%prog(nnow(jg)), &
       &                  ext_data(jg))

     END SELECT

  ENDIF
ENDDO

  END SUBROUTINE init_testcase

!-------------------------------------------------------------------------
END MODULE mo_ha_testcases
