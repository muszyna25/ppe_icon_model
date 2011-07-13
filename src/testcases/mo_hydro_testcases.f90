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
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!! 
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!! 
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!! 
!! 
MODULE mo_hydro_testcases
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
  USE mo_exception,       ONLY: message, finish
  USE mo_impl_constants,  ONLY: SUCCESS, MAX_CHAR_LENGTH, TRACER_ONLY
  USE mo_io_units,        ONLY: nnml, nnml_output
  USE mo_namelist,        ONLY: position_nml, POSITIONED
  USE mo_master_nml,      ONLY: lrestart
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist, &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_model_domain,    ONLY: t_patch
  USE mo_ext_data,        ONLY: ext_data
  USE mo_model_domain_import,ONLY: n_dom
  USE mo_interpolation,   ONLY: t_int_state
  USE mo_parallel_configuration,  ONLY: nproma
  USE mo_run_config,      ONLY: num_lev, ntracer, ltransport, iqv
  USE mo_dynamics_config, ONLY: ltwotime,itime_scheme,lshallow_water,&
                                lcoriolis,nnow,nold,nnew
       
  USE mo_grid_configuration, ONLY :  global_cell_type
  
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
  USE mo_mpi,             ONLY: p_pe, p_io
  USE mo_ldf_init,        ONLY: ldf_init_prog_state
  
  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$' 
 
  PUBLIC   ! mo_hydro_testcase variables referenced in mo_io_vlist
           ! --> change from PRIVATE to PUBLIC

  CHARACTER(len=MAX_CHAR_LENGTH),PUBLIC :: ctest_name  ! Test case specifier
  CHARACTER(len=MAX_CHAR_LENGTH),PUBLIC :: ape_sst_case !SST for APE experiments
  

                            
  REAL(wp) :: rotate_axis_deg
  INTEGER  :: ihs_init_type
  LOGICAL :: lhs_vn_ptb
  REAL(wp) :: hs_vn_ptb_scale

  LOGICAL :: lrh_linear_pres  !< if .TRUE., initialize the relative humidity using 
                               !< a linear function of pressure
  LOGICAL :: linit_tracer_fv  !< finite volume initialization for tracer fields
                               !< if .TRUE.
  REAL(wp) :: rh_at_1000hpa    !< relative humidity [0,1] at 1000 hPa.

  INTEGER  :: ildf_init_type   ! isothermal atmosphere at rest (=0, default),
                               ! JWs zonal state (=1) 
  LOGICAL  :: ldf_symm         ! if .TRUE. local diabatic forcing symmetric
                               ! about the equator. if .FALSE. the forcing is placed
                               ! at 30 N. 

! Control variables setting up the configuration of the test run

  NAMELIST/testcase_ctl/ ctest_name, rotate_axis_deg, ape_sst_case,     &
    &                    gw_brunt_vais, gw_u0, gw_lon_deg, gw_lat_deg,  &
    &                    rh_wavenum, rh_init_shift_deg,                 &
    &                    mountctr_lon_deg, mountctr_lat_deg,            &
    &                    mountctr_height, mount_half_width,             &
    &                    jw_uptb, mount_u0, ihs_init_type, lhs_vn_ptb,  &
    &                    hs_vn_ptb_scale, lrh_linear_pres,              &
    &                    rh_at_1000hpa, linit_tracer_fv, ldf_symm,      &
    &                    ildf_init_type

  PUBLIC  :: setup_testcase, init_testcase, rotate_axis_deg


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
 SUBROUTINE setup_testcase
!
! !local variable
   INTEGER :: istat, funit
   INTEGER :: nlev            !< number of full levels

   CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = &
                                 & '(mo_hydro_testcases) setup_testcase'

!-----------------------------------------------------------------------
   nlev = num_lev(1)

!
! set up the default values
!
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

   ape_sst_case           = 'sst1'

   IF (global_cell_type == 3) THEN
     linit_tracer_fv   = .TRUE. ! finite volume initialization for tracer
   ELSE
     linit_tracer_fv   = .FALSE.
   ENDIF

   ildf_init_type     = 0      ! isothermal atmosphere at rest
   ldf_symm           = .TRUE. ! forcing symmetric about the equator

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values in the previous integration.
    !----------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('testcase_ctl')
      READ(funit,NML=testcase_ctl)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL position_nml ('testcase_ctl', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, testcase_ctl)
    END SELECT

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=testcase_ctl)
    CALL store_and_close_namelist(funit, 'testcase_ctl')

    ! Write the contents of the namelist to an ASCII file.
    ! Probably will be removed later.

    IF(p_pe == p_io) WRITE(nnml_output,nml=testcase_ctl)


   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   IF (.NOT. lshallow_water) THEN

      SELECT CASE (TRIM(ctest_name))

     CASE ('GW')
        CALL message(TRIM(routine),'running the gravity wave test.')
      CASE ('MRW')
        CALL message(TRIM(routine),'running the mountain induced Rossby &
                                   &wave train.')
      CASE ('MRW2')
        CALL message(TRIM(routine),'running the modified mountain induced Rossby &
                                   &wave train (smaller-scale mountain).')
      CASE ('RH')
        CALL message(TRIM(routine),'running the 3D Rossby-Haurwitz wave test.')
      CASE ('JWs')
        CALL message(TRIM(routine),'running the Jablonowski-Williamson &
                                   &steady solution test.')
      CASE ('JWw')
        CALL message(TRIM(routine),'running the Jablonowski-Williamson &
                                   &baroclinic wave test.')
      CASE ('PA')
        CALL message(TRIM(routine),'running the Pure 3D-Advection test.')

      CASE ('SV')
        CALL message(TRIM(routine),'running the stationary vortex 2D-Advection test.')

      CASE ('DF1')
        CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 1.')

      CASE ('DF2')
        CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 2.')

      CASE ('DF3')
        CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 3.')

      CASE ('DF4')
        CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 4.')

      CASE ('HS')
        CALL message(TRIM(routine),'running the Held-Suarez test.')

      CASE ('LDF')
        CALL message(TRIM(routine),'running the local diabatic forcing test.')

      CASE ('LDF-Moist')
        CALL message(TRIM(routine),'running the local diabatic forcing test &
                                   & with moist processes.')
      CASE ('APE')
        CALL message(TRIM(routine),'running the aqua-planet simulation.')

      CASE ('JWw-Moist')
        CALL message(TRIM(routine),'running the Jablonowski-Williamson baroclinic &
                                   &wave test with moist processes.')
      CASE default
        CALL finish(TRIM(routine),'unknown choice of CTEST_NAME')

      END SELECT

   ELSE ! shallow water case

      SELECT CASE (TRIM(ctest_name))

      CASE ('Will_2')
        CALL message(TRIM(routine),'running Williamson test no. 2')
      CASE ('Will_3')
        CALL message(TRIM(routine),'running Williamson test no. 3')
      CASE ('Will_5')
        CALL message(TRIM(routine),'running Williamson test no. 5')
      CASE ('Will_6')
        CALL message(TRIM(routine),'running Williamson test no. 6')
      CASE ('USBR')
        CALL message(TRIM(routine),'running unsteady solid body rotation')
      CASE ('SW_GW')
        CALL message(TRIM(routine),'running shallow water gravity wave test')
      CASE ('PA')
        CALL message(TRIM(routine),'running the Pure 2D-Advection test.')
      CASE default
        CALL finish(TRIM(routine),'unknown choice of CTEST_NAME')
      END SELECT

   ENDIF

   !---------------------------------------------------------------------------
   ! Cross-check (should be moved somewhere else)
   !---------------------------------------------------------------------------

   IF ((TRIM(ctest_name)=='GW') .AND. (nlev /= 20)) THEN
     CALL finish(TRIM(routine),'nlev MUST be 20 for the gravity-wave test case')
   ENDIF
   IF ((TRIM(ctest_name)/='GW') .AND. (nlev == 20)) THEN
     CALL finish(TRIM(routine),'nlev=20 is allowed ONLY for the gravity-wave test case')
   ENDIF
   IF ((TRIM(ctest_name)=='SV') .AND. ntracer /= 2 ) THEN
     CALL finish(TRIM(routine), &
       & 'ntracer MUST be 2 for the stationary vortex test case')
   ENDIF 
   IF ((TRIM(ctest_name)=='DF1') .AND. ntracer == 1 ) THEN
     CALL finish(TRIM(routine), &
       & 'ntracer MUST be >=2 for the deformational flow test case 1')
   ENDIF 
   IF ((TRIM(ctest_name)=='DF2') .AND. ntracer == 1 ) THEN
     CALL finish(TRIM(routine), &
       & 'ntracer MUST be >=2 for the deformational flow test case 2')
   ENDIF 
   IF ((TRIM(ctest_name)=='DF3') .AND. ntracer == 1 ) THEN
     CALL finish(TRIM(routine), &
       & 'ntracer MUST be >=2 for the deformational flow test case 3')
   ENDIF 
   IF ((TRIM(ctest_name)=='DF4') .AND. ntracer == 1 ) THEN
     CALL finish(TRIM(routine), &
       & 'ntracer MUST be >=2 for the deformational flow test case 4')
   ENDIF     


END SUBROUTINE setup_testcase
!-------------------------------------------------------------------------
!

  !>
  !!               Initialization routine for the testcases.
  !! 
  !!               Initialization routine for the testcases
  !!               of the 3D hydrostatic test cases.
  !! 
  !! @par Revision History
  !!  Hui Wan, MPI-M (2007-07-26)
  !! 
  SUBROUTINE init_testcase(pt_patch, pt_hydro_state, pt_int_state)


  TYPE(t_patch),INTENT(inout),TARGET :: pt_patch(n_dom)
  TYPE(t_int_state)   :: pt_int_state(n_dom)
  TYPE(t_hydro_atm) :: pt_hydro_state(n_dom)

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = &
                     & '(mo_hydro_testcases) init_testcase:'

  INTEGER :: jg
  INTEGER :: nlev              !< number of full levels
  INTEGER :: ist     !< status variable
!-------------------------------------------------

DO jg = 1,n_dom

  nlev = pt_patch(jg)%nlev

  IF (lshallow_water) THEN

    IF (TRIM(ctest_name)=='PA') THEN
      !
      IF ( .NOT.ltwotime .OR. (itime_scheme /= TRACER_ONLY) ) THEN 
        CALL finish(TRIM(routine),'running the model in  &
          & Pure Advection mode requires ltwotime=.TRUE. and &
          & itime_scheme=1.')
      END IF
      !
      CALL init_hydro_state_prog_patest(pt_patch(jg),       &
           &           pt_hydro_state(jg)%prog(nnow(jg)),   &
           &           pt_hydro_state(jg)%diag,             &
           &           pt_int_state(jg), ext_data(jg),      &
           &           rotate_axis_deg, linit_tracer_fv)

      CALL init_hydro_state_prog_patest(pt_patch(jg),       &
           &           pt_hydro_state(jg)%prog(nnew(jg)),   &
           &           pt_hydro_state(jg)%diag,             &
           &           pt_int_state(jg), ext_data(jg),      &
           &           rotate_axis_deg, linit_tracer_fv)
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
             & rotate_axis_deg)
        CALL init_hydro_state_prog_jwtest(pt_patch(jg), &
             & pt_hydro_state(jg)%prog(nnow(jg)),       &
             & pt_hydro_state(jg)%diag, ext_data(jg),   &
             & rotate_axis_deg)

     CASE ('JWw')

        IF (.NOT.ltwotime) &
             & CALL init_hydro_state_prog_jwtest(pt_patch(jg),&
             & pt_hydro_state(jg)%prog(nold(jg)),             &
             & pt_hydro_state(jg)%diag, ext_data(jg),         &
             & rotate_axis_deg)
        CALL init_hydro_state_prog_jwtest(pt_patch(jg), &
             & pt_hydro_state(jg)%prog(nnow(jg)),       &
             & pt_hydro_state(jg)%diag,  ext_data(jg),  &
             & rotate_axis_deg)

     CASE ('PA')
        !
        IF ( .NOT.ltwotime .OR. (itime_scheme /= TRACER_ONLY) )           &
             & CALL finish(TRIM(routine),'running the model in  &
             & Pure Advection mode requires ltwotime=.TRUE. and &
             & itime_scheme=1.')
        !
        CALL init_hydro_state_prog_patest(pt_patch(jg),       &
             &           pt_hydro_state(jg)%prog(nnow(jg)),   &
             &           pt_hydro_state(jg)%diag,             &
             &           pt_int_state(jg), ext_data(jg),      &
             &           rotate_axis_deg, linit_tracer_fv)

        CALL init_hydro_state_prog_patest(pt_patch(jg),       &
             &           pt_hydro_state(jg)%prog(nnew(jg)),   &
             &           pt_hydro_state(jg)%diag,             &
             &           pt_int_state(jg), ext_data(jg),      &
             &           rotate_axis_deg, linit_tracer_fv)


     CASE ('SV')
        !
        IF ( .NOT.ltwotime .OR. (itime_scheme /= TRACER_ONLY) )           &
             & CALL finish(TRIM(routine),'running the model in  &
             & Stationary Vortex mode requires ltwotime=.TRUE. and &
             & itime_scheme=1.')
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
        IF ( .NOT.ltwotime .OR. (itime_scheme /= TRACER_ONLY) ) &
             & CALL finish(TRIM(routine),'running the model in  &
             & deformational flow mode requires ltwotime=.TRUE. and &
             & itime_scheme=1.')
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
             &           linit_tracer_fv )

        CALL init_hydro_state_prog_dftest(pt_patch(jg),       &
             &           pt_hydro_state(jg)%prog(nnew(jg)),   &
             &           pt_hydro_state(jg)%diag,             &
             &           pt_int_state(jg), ext_data(jg),      &
             &           rotate_axis_deg, ctest_name,         &
             &           linit_tracer_fv )
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
                & rotate_axis_deg)
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
             & rotate_axis_deg,                         &
             & lrh_linear_pres, rh_at_1000hpa )

        IF (.NOT.ltwotime) CALL copy_prog_state(  &
             & pt_hydro_state(jg)%prog(nnow(jg)), & !in
             & pt_hydro_state(jg)%prog(nold(jg)), & !out
             & .FALSE.,     &! copy temp rather than theta
             & ltransport )  ! copy tracer field if transport is on.

     CASE ('APE')
        ! Initial conditions are the same as for the 'JWw-Moist' case

        CALL init_hydro_state_prog_jwtest(pt_patch(jg), &
             & pt_hydro_state(jg)%prog(nnow(jg)),       &
             & pt_hydro_state(jg)%diag,  ext_data(jg),  &
             & rotate_axis_deg,                         &
             & lrh_linear_pres, rh_at_1000hpa )

        pt_hydro_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqv+1:) = 0._wp

        IF (.NOT.ltwotime) CALL copy_prog_state(  &
             & pt_hydro_state(jg)%prog(nnow(jg)), & !in
             & pt_hydro_state(jg)%prog(nold(jg)), & !out
             & .FALSE.,     &! copy temp rather than theta
             & ltransport )  ! copy tracer field if transport is on.

        ! Surface geopotential height is constantly zero.

        ext_data(jg)%atm%topography_c = 0._wp

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
END MODULE mo_hydro_testcases
