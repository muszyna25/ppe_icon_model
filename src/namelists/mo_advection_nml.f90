!>
!! Namelist for scalar transport
!!
!! These subroutines are called by  read_atmo_namelists and do the transport 
!! setup.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-04-20)
!! - moved here from mo_advection_utils
!! Modification by Daniel Reinert, DWD (2011-04-20)
!! - some updates on the structure
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_advection_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_run_config,          ONLY: ntracer
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_ntracer, max_dom,      &
    &                               MIURA, FFSL_HYB_MCYCL, ippm_vcfl, ippm_v,   &
    &                               inol, ifluxl_sm, inol_v,                    &
    &                               islopel_vsm, ifluxl_vpd
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_advection_config,    ONLY: advection_config 
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_transport_namelist




  !----------------------------------!
  ! transport_nml namelist variables !
  !----------------------------------!

  CHARACTER(len=MAX_CHAR_LENGTH) :: &!< list of tracers to initialize
    &  ctracer_list


  INTEGER :: &                     !< selects horizontal transport scheme
    &  ihadv_tracer(max_ntracer)   !< 0:  no horizontal advection
                                   !< 1:  1st order upwind
                                   !< 2:  2nd order miura
                                   !< 3:  3rd order miura with quadr./cubic reconstr.
                                   !< 4:  Flux form semi lagrange (FFSL)
                                   !< 5:  hybrid FFSL Miura3
                                   !< 20: subcycling version of miura
                                   !< 22: 2nd order miura and miura_cycl
                                   !< 32: 3rd order miura with miura_cycl
                                   !< 42: FFSL with miura_cycl
                                   !< 52: FFSL_HYB with miura_cycl

  INTEGER :: &                     !< selects vertical transport scheme
    &  ivadv_tracer(max_ntracer)   !< 0 : no vertical advection
                                   !< 1 : 1st order upwind
                                   !< 3 : 3rd order PPM for CFL>
                                   !< 30: 3rd order PPM

  INTEGER :: &                     !< advection of TKE
    &  iadv_tke                    !< 0 : none
                                   !< 1 : vertical advection only
                                   !< 2 : vertical and horizontal advection


  LOGICAL :: lvadv_tracer          !< if .TRUE., calculate vertical tracer advection
  LOGICAL :: lclip_tracer          !< if .TRUE., clip negative tracer values
  LOGICAL :: lstrang               !< if .TRUE., use complete Strang splitting
                                   !< (\Delta t/2 vert)+(\Delta t hor)+(\Delta t/2 vert)

  LOGICAL :: llsq_svd              !< least squares reconstruction with 
                                   !< singular value decomposition (TRUE) or 
                                   !< QR decomposition (FALSE) of design matrix A

  INTEGER :: &                     !< parameter used to select the limiter
    &  itype_vlimit(max_ntracer)   !< for vertical transport
                               

  INTEGER :: &                     !< parameter used to select the limiter
    &  itype_hlimit(max_ntracer)   !< for horizontal transport
                                   !< 0: no limiter
                                   !< 3: monotonous flux limiter
                                   !< 4: positive definite flux limiter

  INTEGER :: niter_fct             !< number of iterations for monotone
                                   !< flux correction procedure

  REAL(wp):: beta_fct              !< factor of allowed over-/undershooting in monotonous limiter


  INTEGER :: iord_backtraj         !< parameter to select the spacial order
                                   !< of accuracy for the backward trajectory

  INTEGER :: igrad_c_miura         !< parameter used to select the gradient
                                   !< reconstruction method at cell center
                                   !< for second order miura scheme

  INTEGER :: ivcfl_max             !< determines stability range of vertical 
                                   !< ppm-scheme (approximate allowable maximum 
                                   !< CFL-number)

  INTEGER :: npassive_tracer       !< number of additional passive tracers, in addition to
                                   !< microphysical- and ART tracers. 

  CHARACTER(len=MAX_CHAR_LENGTH) :: &!< Comma separated list of initialization formulae 
    &  init_formula                  !< for passive tracers.

  NAMELIST/transport_nml/ ihadv_tracer, ivadv_tracer, lvadv_tracer,       &
    &                     itype_vlimit, ivcfl_max, itype_hlimit,          &
    &                     iadv_tke, niter_fct, beta_fct, iord_backtraj,   &
    &                     lclip_tracer, ctracer_list, igrad_c_miura,      &
    &                     lstrang, llsq_svd, npassive_tracer,             &
    &                     init_formula


CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for tracer transport. 
  !!
  !! This subroutine 
  !! - reads the Namelist for tracer transport
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  !! @par Revision History
  !! Initial Revision by Daniel Reinert, DWD (2011-05-07)
  !!
  SUBROUTINE read_transport_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_advection_nml: read_transport_nml'

    !-----------------------
    ! 1. default settings   
    !-----------------------
    ctracer_list    = ''
    ihadv_tracer(:) = MIURA     ! miura horizontal advection scheme
    itype_hlimit(:) = ifluxl_sm ! positive definite flux limiter
    ivadv_tracer(:) = ippm_vcfl ! PPM vertical advection scheme
    itype_vlimit(:) = islopel_vsm ! semi-monotonous slope limiter
    iadv_tke        = 0         ! no TKE advection
    niter_fct       = 1         ! number of FCT-iterations
    beta_fct        = 1.005_wp  ! factor of allowed over-/undershooting in monotonous limiter
    ivcfl_max       = 5         ! CFL-stability range for vertical advection
    iord_backtraj   = 1         ! 1st order backward trajectory
    lvadv_tracer    = .TRUE.    ! vertical advection yes/no
    lclip_tracer    = .FALSE.   ! clipping of negative values yes/no
    lstrang         = .FALSE.   ! Strang splitting yes/no

    igrad_c_miura   = 1         ! MIURA linear least squares reconstruction

    llsq_svd        = .TRUE.    ! apply singular-value-decomposition (FALSE: use QR-decomposition)
    npassive_tracer = 0         ! no additional passive tracers
    init_formula    = ''        ! no explizit initialization of passive tracers


    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('transport_nml')
      READ(funit,NML=transport_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('transport_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, transport_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, transport_nml)                                        ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, transport_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------

    ! flux computation methods - sanity check
    !
    IF ( ANY(ihadv_tracer(1:ntracer) > FFSL_HYB_MCYCL) .OR.            &
      &  ANY(ihadv_tracer(1:ntracer) < 0) )    THEN
      CALL finish( TRIM(routine),                                     &
        &  'incorrect settings for ihadv_tracer. Must be 0,1,2,3,4,5,'//&
        &  '20,22,32,42 or 52 ')
    ENDIF
    IF ( ANY(ivadv_tracer(1:ntracer) > ippm_v) .OR.                   &
      &  ANY(ivadv_tracer(1:ntracer) < 0)) THEN
      CALL finish( TRIM(routine),                                     &
        &  'incorrect settings for ivadv_tracer. Must be 0,1,2,3,20, or 30 ')
    ENDIF



    ! limiter - sanity check
    !
    IF ( ANY(itype_vlimit(1:ntracer) < inol_v ) .OR.                  &
      &  ANY(itype_vlimit(1:ntracer) > ifluxl_vpd)) THEN
      CALL finish( TRIM(routine),                                     &
        &  'incorrect settings for itype_vlimit. Must be 0,1,2 or 4 ')
    ENDIF
    IF ( ANY(itype_hlimit(1:ntracer) < inol ) .OR.                    &
      &  ANY(itype_hlimit(1:ntracer) > ifluxl_sm)) THEN
      CALL finish( TRIM(routine),                                     &
        &  'incorrect settings for itype_hlimit. Must be 0,1,2,3 or 4 ')
    ENDIF


    ! FCT-iterations - sanity check
    !
    IF ( niter_fct < 1 ) THEN
      CALL finish( TRIM(routine), 'niter_fct must be greater than 0 ')
    ENDIF


    ! FCT multiplicative spreading - sanity check
    ! 
    IF ( (beta_fct < 1.0_wp) .OR. (beta_fct >= 2.0_wp) ) THEN
      CALL finish( TRIM(routine),                                     & 
        & 'permissible range of values for beta_fct:  [1,2[ ')
    ENDIF


    ! TKE advection - sanity check
    !
    IF ( ALL((/0,1,2/) /= iadv_tke) ) THEN
      CALL finish( TRIM(routine),                                     &
        &  'incorrect settings for iadv_tke. Must be 0, 1 or 2')
    ENDIF

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 0,max_dom
      advection_config(jg)%ctracer_list   = ctracer_list
      advection_config(jg)%ihadv_tracer(:)= ihadv_tracer(:)
      advection_config(jg)%ivadv_tracer(:)= ivadv_tracer(:)
      advection_config(jg)%lvadv_tracer   = lvadv_tracer
      advection_config(jg)%lclip_tracer   = lclip_tracer
      advection_config(jg)%lstrang        = lstrang
      advection_config(jg)%llsq_svd       = llsq_svd
      advection_config(jg)%itype_vlimit(:)= itype_vlimit(:)
      advection_config(jg)%itype_hlimit(:)= itype_hlimit(:)
      advection_config(jg)%niter_fct      = niter_fct
      advection_config(jg)%beta_fct       = beta_fct
      advection_config(jg)%iord_backtraj  = iord_backtraj
      advection_config(jg)%igrad_c_miura  = igrad_c_miura
      advection_config(jg)%iadv_tke       = iadv_tke
      advection_config(jg)%ivcfl_max      = ivcfl_max
      advection_config(jg)%npassive_tracer= npassive_tracer 
      advection_config(jg)%init_formula   = init_formula 
    ENDDO



    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=transport_nml)                    
      CALL store_and_close_namelist(funit, 'transport_nml')             
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=transport_nml)


  END SUBROUTINE read_transport_namelist

END MODULE mo_advection_nml
