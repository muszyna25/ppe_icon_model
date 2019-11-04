!>
!! This module controls the namelist input for:
!! - Upper-atmosphere physics
!! - Deep-atmosphere dynamics
!! - Upper-atmosphere extrapolation
!!
!! @author Guidi Zhou, MPI-M, 2016-03-03
!!         Sebastian Borchert, DWD, 2016-03-03
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_upatmo_nml

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: finish, message
  USE mo_upatmo_config,           ONLY: upatmo_dyn_config, upatmo_exp_config, &
    &                                   t_upatmo_config, imsg_thr, itmr_thr,  &
    &                                   istatus
  USE mo_namelist,                ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_io_units,                ONLY: nnml, nnml_output
  USE mo_mpi,                     ONLY: my_process_is_stdio
  USE mo_master_control,          ONLY: use_restart_namelists
  USE mo_restart_namelist,        ONLY: open_tmpfile, store_and_close_namelist, &
    &                                   open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,            ONLY: temp_defaults, temp_settings
  USE mo_impl_constants,          ONLY: SUCCESS, MAX_CHAR_LENGTH, max_dom, &
    &                                   inh_atmosphere, inoforcing, inwp, iecham, ivexpol
  USE mo_name_list_output_types,  ONLY: t_output_name_list
  USE mo_name_list_output_config, ONLY: is_variable_in_output
  USE mo_util_string,             ONLY: int2string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_upatmo_namelist
  PUBLIC :: check_upatmo

  CHARACTER(LEN = *), PARAMETER :: modname = 'mo_upatmo_nml'

  !------------------------------------------------------------
  !                    Namelist parameters
  !------------------------------------------------------------

  !---------------
  !   Dynamics
  !---------------  

  ! Specifieres for metric and approximations used in model equations
  !
  LOGICAL :: lnontrad        ! .TRUE. -> non-traditional deep-atmosphere terms 
                             ! in components of momentum equation are switched on
  LOGICAL :: lconstgrav      ! .TRUE. -> gravitational acceleration is const. 
                             ! like in case of the shallow atmosphere
  LOGICAL :: lcentrifugal    ! .TRUE. -> explicit centrifugal acceleration is switched on
  LOGICAL :: ldeepatmo2phys  ! .TRUE. -> the input fields to the ECHAM physics parameterizations 
                             ! are modified for the deep atmosphere, if required 
                             ! .FALSE. -> the input fields are computed in accordance 
                             ! with the shallow-atmosphere approximation (standard) in any case

  !---------------
  ! Extrapolation
  !---------------

  ! Variables for extrapolation methods
  !
  REAL(wp) :: expol_start_height     ! [m] Height above which extrapolation (blending) of initial data starts
  REAL(wp) :: expol_blending_scale   ! [m] Blending scale height
  REAL(wp) :: expol_vn_decay_scale   ! [m] Scale height for (exponential) decay of extrapolated 
                                     ! horizontal wind component (for stability reasons)
  REAL(wp) :: expol_temp_infty       ! [K] Climatological temperature of exosphere (z -> infinity)  
  LOGICAL  :: lexpol_sanitycheck     ! .TRUE. -> Apply sanity check to extrapolated fields

  !------------------------------------------------------------

  NAMELIST /upatmo_nml/ lnontrad,              &   
    &                   lconstgrav,            &   
    &                   lcentrifugal,          &
    &                   ldeepatmo2phys,        &
    &                   expol_start_height,    & 
    &                   expol_blending_scale,  &
    &                   expol_vn_decay_scale,  &
    &                   expol_temp_infty,      &
    &                   lexpol_sanitycheck

  !------------------------------------------------------------

CONTAINS !..................................................................................

  !------------------------------------------------------------
  !                       Subroutines
  !------------------------------------------------------------

  !>
  !! Read Namelist for upper atmosphere. 
  !!
  !! This subroutine: 
  !! - Reads the Namelist for upper atmosphere
  !! - Sets default values
  !! - Potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - Reads the user's (new) specifications
  !! - Does a consistency check
  !! - Stores the Namelist for restart
  !! - Fills the configuration state (partly) 
  !!
  SUBROUTINE read_upatmo_namelist( filename )

    ! In/out variables
    CHARACTER(LEN=*), INTENT(IN) :: filename

    ! Local variables
    INTEGER :: istat, funit
    INTEGER :: iunit
    INTEGER :: jg
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':read_upatmo_namelist'

    !------------------------------------------------------------

    !------------------------------------------------------------
    !                    Set default values
    !------------------------------------------------------------

    ! (Where not otherwise stated, the settings apply to all domains)

    !---------------
    !   Dynamics
    !---------------  

    ! (Note: to switch on the deep-atmosphere dynamics set 'dynamics_nml: ldeepatmo = .TRUE.'.)

    lnontrad       = .TRUE.   ! Non-traditional deep-atmosphere terms in components 
                              ! of momentum equation (such as f_t*w or vn/r) are active 
    lconstgrav     = .FALSE.  ! Gravitational acceleration shall be a function of radius 
                              ! grav -> grav*(a/r)^2, where a is Earth's radius and r = a + z
    lcentrifugal   = .FALSE.  ! Explicit centrifugal acceleration switched off
    ldeepatmo2phys = .FALSE.  ! The input fields to the ECHAM physics parameterizations 
                              ! are computed in accordance with the shallow-atmosphere approximation

    !---------------
    ! Extrapolation
    !---------------

    ! (Note: to switch on the upper-atmosphere extrapolation set 'initicon_nml: itype_vert_expol = 2'.)

    expol_start_height   = 70000.0_wp  ! [m] Height above which extrapolation (blending) of initial data starts
    expol_blending_scale = 10000.0_wp  ! [m] Blending scale height
    expol_vn_decay_scale = 10000.0_wp  ! [m] Horizontal wind decay scale height
    expol_temp_infty     = 400.0_wp    ! [K] Exospheric mean reference temperature 
                                       ! (e.g. climatological value from Hedin (1983) would be 1035 K)
    lexpol_sanitycheck   = .FALSE.     ! No sanity check of extrapolated fields 
                                       ! (bacause it is computationally expensive)

    !------------------------------------------------------------
    !  If this is a resumed integration, overwrite the defaults 
    !     above by values used in the previous integration
    !------------------------------------------------------------
    
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('upatmo_nml')
      READ(funit,NML=upatmo_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------
    !             Read user's (new) specifications 
    !            (Done so far by all MPI processes)
    !------------------------------------------------------------

    CALL open_nml(TRIM(filename))
    CALL position_nml ('upatmo_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, upatmo_nml)    ! Write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, upatmo_nml)     ! Overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, upatmo_nml)  ! Write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !------------------------------------------------------------
    !                     Consistency checks
    !------------------------------------------------------------

    !---------------
    ! Extrapolation
    !---------------

    IF( expol_start_height < 0._wp ) THEN 
      ! The height above which the extrapolation starts
      ! has to be equal to or greater than zero
      CALL finish( TRIM(routine), & 
        & 'Invalid value for expol_start_height, it has to be >= 0.' )
    ELSEIF( expol_blending_scale < 0._wp ) THEN 
      ! The blending scale height has to be equal to 
      ! or greater than zero
      CALL finish( TRIM(routine), & 
        & 'Invalid value for expol_blending_scale, it has to be >= 0.' )
    ELSEIF( expol_vn_decay_scale < 0._wp ) THEN 
      ! Likewise the decay scale height for the horizontal wind
      CALL finish( TRIM(routine), & 
        & 'Invalid value for expol_vn_decay_scale, it has to be >= 0.' )
    ELSEIF( expol_temp_infty < 0._wp ) THEN 
      ! The exospheric mean temperature has to be equal to 
      ! or greater than zero
      CALL finish( TRIM(routine), & 
        & 'Invalid value for expol_temp_infty, it has to be >= 0.' )
    ENDIF

    !------------------------------------------------------------
    !              Fill the configuration state
    !------------------------------------------------------------

    DO jg = 0, max_dom
      
      !---------------
      !   Dynamics
      !---------------  

      upatmo_dyn_config(jg)%lnontrad       = lnontrad 
      upatmo_dyn_config(jg)%lconstgrav     = lconstgrav
      upatmo_dyn_config(jg)%lcentrifugal   = lcentrifugal
      upatmo_dyn_config(jg)%ldeepatmo2phys = ldeepatmo2phys
      
      !---------------
      ! Extrapolation
      !---------------

      upatmo_exp_config(jg)%expol_start_height   = expol_start_height  
      upatmo_exp_config(jg)%expol_blending_scale = expol_blending_scale
      upatmo_exp_config(jg)%expol_vn_decay_scale = expol_vn_decay_scale
      upatmo_exp_config(jg)%expol_temp_infty     = expol_temp_infty  
      upatmo_exp_config(jg)%lexpol_sanitycheck   = lexpol_sanitycheck

    ENDDO  !jg

    !------------------------------------------------------------
    !              Store the namelist for restart
    !------------------------------------------------------------

    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=upatmo_nml)
      CALL store_and_close_namelist(funit, 'upatmo_nml') 
    ENDIF

    !------------------------------------------------------------
    !     Write the contents of the namelist to an ASCII file
    !------------------------------------------------------------

    IF(my_process_is_stdio()) WRITE(nnml_output,nml=upatmo_nml)        

  END SUBROUTINE read_upatmo_namelist

  !------------------------------------------------------------

  !>
  !! Check for conflicts with other namelist settings.
  !! (Note:
  !!  - This subroutine was placed in this module, to avoid circular dependencies,
  !!    it is called in 'src/configure_model/mo_nml_crosscheck'
  !!  - Please, be careful when you change something in the sequence of events below)
  !!
  SUBROUTINE check_upatmo( n_dom_start,            & !in
    &                      n_dom,                  & !in
    &                      iequations,             & !in
    &                      iforcing,               & !in
    &                      ldeepatmo,              & !in
    &                      is_plane_torus,         & !in
    &                      l_limited_area,         & !in
    &                      lart,                   & !in
    &                      ivctype,                & !in
    &                      flat_height,            & !in
    &                      itype_vert_expol,       & !in
    &                      ltestcase,              & !in
    &                      nh_test_name,           & !in
    &                      first_output_name_list, & !in
    &                      upatmo_config           ) !inout

    ! In/out variables
    INTEGER,                               INTENT(IN)    :: n_dom_start            ! Start index of domains
    INTEGER,                               INTENT(IN)    :: n_dom                  ! End index of domains
    INTEGER,                               INTENT(IN)    :: iequations             ! Switch for model equations 
                                                                                   ! (non-hydrostatic etc.)
    INTEGER,                               INTENT(IN)    :: iforcing               ! Switch for physics package 
                                                                                   ! (nwp, echam etc.)  
    LOGICAL,                               INTENT(IN)    :: ldeepatmo              ! Switch for deep-atmosphere dynamics
    LOGICAL,                               INTENT(IN)    :: is_plane_torus         ! Switch for torus mode
    LOGICAL,                               INTENT(IN)    :: l_limited_area         ! Switch for limited-area mode
    LOGICAL,                               INTENT(IN)    :: lart                   ! Switch for ART interface
    INTEGER,                               INTENT(IN)    :: ivctype                ! Type of vertical grid (SLEVE etc.)
    REAL(wp),                              INTENT(IN)    :: flat_height            ! Below 'flat_height' grid layer 
                                                                                   ! interfaces follow topography
    INTEGER,                               INTENT(IN)    :: itype_vert_expol       ! Type of vertical extrapolation 
                                                                                   ! of initial atmosphere state
    LOGICAL,                               INTENT(IN)    :: ltestcase              ! Switch for test case mode
    CHARACTER(LEN=*),                      INTENT(IN)    :: nh_test_name           ! Test case name
    TYPE(t_output_name_list),              POINTER       :: first_output_name_list ! Pointer to a linked list 
                                                                                   ! of output name lists
    TYPE(t_upatmo_config),    ALLOCATABLE, INTENT(INOUT) :: upatmo_config(:)       ! Upper-atmosphere configuration type

    ! Local variables
    INTEGER :: istat
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':check_upatmo'

    !------------------------------------------------------------

    !---------------
    !   Dynamics
    !---------------  

    IF (ldeepatmo) THEN

      ! If deep-atmosphere dynamics have been switched on ...
      IF (iequations /= inh_atmosphere) THEN
        ! ... only the non-hydrostatic set of equations is allowed
        CALL finish(TRIM(routine), &
          & "Deep-atmosphere configuration is not available for other than the non-hydrostatic equations.")
      ELSEIF(.NOT. ANY((/inoforcing, inwp, iecham/) == iforcing)) THEN
        ! ... only no physics forcing, ECHAM forcing or NWP forcing are allowed
        CALL finish(TRIM(routine), &
          & "Deep-atmosphere configuration is not available for all forcings but iforcing = "// &
          & TRIM(int2string(inoforcing))//', or '//TRIM(int2string(iecham))//', or '//TRIM(int2string(inwp))//'.')
      ELSEIF (ltestcase .AND. (.NOT. (TRIM(nh_test_name) == 'dcmip_bw_11' .OR. TRIM(nh_test_name) == 'lahade'))) THEN
        ! ... most test cases are not available for the time being
        ! (only the baroclinic wave test case of Ullrich et al. (2014),
        ! and the lahade-testcase are currently intended to test the deep-atmosphere equations)
        CALL finish(TRIM(routine), &
          & "Deep-atmosphere configuration is not available for all test cases but dcmip_bw_11 and lahade.")
      ELSEIF (is_plane_torus) THEN 
        ! ... the torus configuration is not available for the time being (-> no spherical geometry)
        CALL finish(TRIM(routine), &
          & "Deep-atmosphere configuration is not available in combination with the torus mode.")
      ENDIF

      ! ... the limited-area mode requires a warning for the time being
      IF (l_limited_area) THEN 
        CALL message(TRIM(routine), "WARNING, are the deep-atmosphere dynamics really necessary,"// &
          & " and consistent with the driving model?")
      ENDIF

      ! ... a run in combination with ART necessitates at least a warning for the time being
      IF (lart) THEN 
        CALL message(TRIM(routine), "WARNING, (dynamical) cross-consistency/compatibility of ART"// &
          & " in combination with the deep atmosphere has not been checked!")
      ENDIF

      ! ... the computation of the output variable 'potential vorticity' is not modified 
      ! for the deep atmosphere for the time being
      IF (is_variable_in_output(first_output_name_list, var_name="pv")) THEN
        CALL message(TRIM(routine),'WARNING, PV-computation is not modified for deep atmosphere!')
      ENDIF

      ! ... OpenACC-parallelization is not available in combination with the deep-atmosphere configuration
#ifdef _OPENACC
      CALL finish(TRIM(routine), "Deep-atmosphere configuration is not available in combination with Open-ACC.")
#endif

    ENDIF  !IF (ldeepatmo)
    
    !---------------
    ! Extrapolation
    !---------------

    IF (itype_vert_expol == ivexpol%upatmo) THEN

      IF ((ivctype == 2 .OR. ivctype == 12) .AND. &
        & ANY(upatmo_exp_config(:)%expol_start_height < flat_height)) THEN
        ! Upper-atmosphere extrapolation: start height above which extrapolation takes place 
        ! should not lie below 'flat_height'
        CALL finish(TRIM(routine), &
          & "Upper-atmosphere extrapolation: start height has to be above flat_height.")
      ELSEIF (l_limited_area) THEN
        ! This type of extrapolation is not intended for the limited-area mode
        CALL finish(TRIM(routine), &
          & "The limited-area mode requires: itype_vert_expol = "//int2string(ivexpol%lin))
      ENDIF

    ENDIF  !IF (itype_vert_expol == ivexpol%upatmo)

    !---------------
    ! Miscellaneous
    !---------------

    ! To save some memory, 'upatmo_config' was made allocatable. Unfortunately this entails some complications:
    ! - It cannot be allocated in 'read_upatmo_namelist' above, since 'n_dom_start' and 'n_dom' might not yet be available
    ! - An allocation in 'src/configure_model/mo_upatmo_config: configure_upatmo' would be too late, 
    !   since 'upatmo_config' is already required in 'src/configure_model/mo_initicon_config: configure_initicon' 
    !   (and below as well)
    ! However, this subroutine is called early enough, so an allocation of 'upatmo_config' here 
    ! is practically without any alternative.
    IF (ALLOCATED(upatmo_config)) THEN
      CALL finish(TRIM(routine), "Error in calling sequence: upatmo_config is already allocated.")
    ELSE
      ! Some of the variables in 'upatmo_config' might be necessary for a coarser radiation grid as well, 
      ! so the index range starts with 'n_dom_start' (which is zero, if a coarser radiation grid is used).
      ALLOCATE(upatmo_config(n_dom_start:n_dom), STAT=istat)
      IF (istat /= SUCCESS) CALL finish(TRIM(routine), "Allocation of upatmo_config failed.")
    ENDIF

    !---------------
    ! Status update
    !---------------
    
    ! Indicate that crosscheck has taken place
    upatmo_config(n_dom_start:n_dom)%dyn%l_status(istatus%checked) = .TRUE.
    upatmo_config(n_dom_start:n_dom)%phy%l_status(istatus%checked) = .TRUE.
    upatmo_config(n_dom_start:n_dom)%exp%l_status(istatus%checked) = .TRUE.
    upatmo_config(n_dom_start:n_dom)%l_status(istatus%checked)     = .TRUE.

  END SUBROUTINE check_upatmo

END MODULE mo_upatmo_nml
