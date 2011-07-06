!>
!! Contains the setup of configuration of the dynamical core
!!
!!        
!! @par Revision History
!!   Revision History in mo_global_variables.f90 (r3611)
!!   Modification by Constantin Junk (2011-02-24)
!!     - added new module mo_dynamics_nml
!!     - separated declaration of namelist dynamics_ctl from 
!!       mo_global_variables and moved it mo_dynamics_nml
!!     - moved subroutine deallocate_timelevs from 
!!       mo_global_variables to mo_dynamics_nml
!!   Modification by Constantin Junk (2011-03-23)
!!     - moved hydrostatic_ctl variables to mo_dynamics_nml
!!     - added the consistency checks for these variables
!!       to setup_dynamics_nml
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_dynamics_nml

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, UNKNOWN, TRACER_ONLY, &
                                   LEAPFROG_EXPL, LEAPFROG_SI, AB2
  USE mo_physical_constants, ONLY: grav
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned
  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_run_nml,            ONLY: ltransport,dtime,ntracer,global_cell_type,      &
    &                              iforcing,lshallow_water,INWP, IHS_ATM_TEMP,&
    &                              IHS_ATM_THETA,iequations
  USE mo_grid_configuration,  ONLY: global_cell_type
  USE mo_master_nml,         ONLY: lrestart
  USE mo_io_restart_attributes, ONLY: get_restart_attribute
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist,   &       
                                 & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  PUBLIC

  !---------------------------------------------------------------
  ! Namelist variables and auxiliary parameters setting up the
  ! configuration of the dynamical core
  !---------------------------------------------------------------
  ! time stepping scheme 

  INTEGER  :: itime_scheme       ! parameter used to select the time stepping scheme
                                 ! = 1, explicit 2 time level scheme
                                 ! = 2, semi implicit 2 time level scheme
                                 ! = 3, explicit leapfrog
                                 ! = 4, leapfrog with semi implicit correction
                                 ! = 5, 4-stage Runge-Kutta method
                                 ! = 6, SSPRK(5,4) (Runge-Kutta) method

  INTEGER  :: ileapfrog_startup  ! choice of the first time step in
                                 ! a leapfrog time stepping scheme
                                 ! 1 = Euler forward
                                 ! 2 = several sub-steps

  REAL(wp) :: asselin_coeff      ! parameter used in Asselin filter

  REAL(wp) :: si_2tls            ! time averaging parameter
                                 ! in two time level semi implicit
                                 ! trapeziodal discretization
                                 ! (aka theta-method!!)
                                 ! WARNING: must be between 0.5 and 1.0
                                 ! for unconditional (linear!!) stability
                                 ! values around 0.6 are safe for all tests
                                 ! but result in larger damping
                                 ! smaller values can be used but

  INTEGER  :: si_expl_scheme     ! scheme for the explicit part of the
                                 ! 2-time-level semi-implicit time integration.
                                 ! See mo_impl_constants for the options.

  REAL(wp) :: si_coeff           !  = 0 : explicit scheme(for *d*,*t*,*alps*).
                                 !  = 1 : semi implicit scheme.
                                 !  in (0,1): a weighted scheme

  REAL(wp) :: si_offctr          ! weighting parameter used in calculating the
                                 ! second temporal derivatives in the semi-implicit
                                 ! correction scheme. The value read from namelist are
                                 ! assumed to be the offcentering (i.e. between 0 and 1).

  LOGICAL  :: lsi_3d             ! if .true., solve the 3D equation

  REAL(wp) :: si_rtol            ! relative tolerance

  REAL(wp) :: si_cmin            ! min. phase speed of the decomposed modes to be
                                 ! solved by the semi-implicit correction scheme

  REAL(wp) :: sw_ref_height      ! reference height to linearize around if using
                                 ! lshallow_water and semi-implicit correction

  ! way of computing the divergence operator in the triangular model -------------

  INTEGER  :: idiv_method  ! 1: Hydrostatic atmospheric model: 
                           !    Gauss integral with original normal velocity components
                           ! 1: Non-Hydrostatic atmospheric model: 
                           !    Gauss integral with averged normal velocity components
                           ! Thus, in a linear equilateral grid, methods 1 and 2 for
                           ! the non-hydrostatic model are the same.
                           ! 2: divergence averaging with bilinear averaging

  REAL(wp) :: divavg_cntrwgt  ! weight of central cell for divergence averaging


  
  LOGICAL :: ldry_dycore ! if .TRUE., ignore the effact of water vapor,
                         ! cloud liquid and cloud ice on virtual temperature.
  LOGICAL :: lref_temp   ! if .TRUE., involve the reference temperature profile
                         ! in the calculation of pressure gradient force.

  NAMELIST/dynamics_ctl/ itime_scheme, ileapfrog_startup,              &
                       & asselin_coeff, si_2tls,                       &
                       & idiv_method, divavg_cntrwgt,                  &
                       & si_coeff, si_offctr, si_rtol,                 &
                       & si_cmin, lsi_3d, si_expl_scheme,              &
                       & ldry_dycore, lref_temp

  !------------------------------------------------------------------------
  ! Dependent variables 
  !------------------------------------------------------------------------

  LOGICAL  :: ltwotime       ! if .TRUE., two time level discretizations are used,
                             ! if .FALSE., three time level discretizations are used
  REAL(wp) :: gdt            ! coefficients in
  REAL(wp) :: gsi_2tlsdt     ! semi-implicit schemes. They are not explicitly
  REAL(wp) :: gsi_2tls2dt    ! specified by the user, but derived from some
  REAL(wp) :: g1msi_2tlsdt   ! of the parameters above.

  ! time level variables

  INTEGER,ALLOCATABLE :: nold(:), nnow(:), nnew(:) ! variables denoting time levels
  INTEGER,ALLOCATABLE :: nsav1(:), nsav2(:)        ! extra 'time levels' of prognostic variables
                                                   ! needed to compute boundary tendencies and
                                                   ! feedback increments
  INTEGER,ALLOCATABLE :: nnow_rcf(:), nnew_rcf(:)  ! extra time levels for reduced
                                                   ! calling frequency (rcf)

CONTAINS

  SUBROUTINE read_dynamics_namelist()

    INTEGER :: istat, funit

    !------------------------------------------------------------
    ! 1. Set up the default values for dynamics_ctl
    !------------------------------------------------------------
 
    itime_scheme      = LEAPFROG_SI
    ileapfrog_startup = 1
    asselin_coeff     = 0.1_wp
 
    si_2tls           = 0.6_wp
    si_expl_scheme    = AB2 
 
    si_coeff          = 1.0_wp
    si_offctr         = 0.7_wp
    lsi_3d            = .FALSE.
    si_rtol           = 1.e-3_wp
    si_cmin           = 30._wp
 
    idiv_method       = 1
    divavg_cntrwgt    = 0.5_wp
 
    SELECT CASE (iequations)
    CASE (IHS_ATM_TEMP, IHS_ATM_THETA)
      ldry_dycore       = .TRUE.
      lref_temp         = .FALSE.
    END SELECT

    !------------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above by 
    !    values in the restart file
    !------------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('dynamics_ctl')
      READ(funit,NML=dynamics_ctl)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL position_nml ('dynamics_ctl', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, dynamics_ctl)
    END SELECT

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=dynamics_ctl)
    CALL store_and_close_namelist(funit, 'dynamics_ctl')

  ! ! write the contents of the namelist to an ASCII file
  ! IF(p_pe == p_io) WRITE(nnml_output,nml=dynamics_ctl)

  END SUBROUTINE read_dynamics_namelist

  !>
  !! @brief Initialization of variables that set up the dynamica core.
  !!
  !! Initialization of variables that set up the configuration
  !! of the dynamical core using values read from
  !! namelist 'dynamics_ctl'.
  !!
  !! @par Revision History
  !!  by Hui Wan, MPI-M (2007-02-23)
  !! Modification by A. Gassmann, MPI-M (2007-06-12)
  !!   changed violating criteria for diffusion coefficients according to stability
  !!   analysis (cf. Durran:"Numerical methods for wave equations in geophysical
  !!   fluid dynamics", pages 136 ff.)
  !! Modification by Jochen Foerstner, DWD (2008-05-06)
  !!   new namelist variables: lsimpson_dyn, lsimpson_tradv to switch on the
  !!   Simpson's rule for the quadrature in the formulation of the divergence
  !!   operator.
  !! Modification by Jochen Foerstner, DWD (2008-07-16)
  !!   new namelist variable: lrbf_vec_int_edge to switch on the vector RBF
  !!   reconstruction at edge midpoints (instead of averaging cell centered
  !!   values to the edges).
  !! Modification by Almut Gassmann, MPI-M (2008-09-23)
  !!  -remove dxmin: estimate the dual edge length by nroot and start_lev
  !!  -clean up dyn_ctl and remove unnecessary variables
  !!
  SUBROUTINE dynamics_nml_setup(i_ndom)

    INTEGER, INTENT(IN) :: i_ndom
 
    INTEGER :: istat, funit, jdom

    CHARACTER(len=MAX_CHAR_LENGTH),PARAMETER ::             &
             & routine = 'mo_dynamics_nml/dynamics_nml_setup'
 
    CHARACTER(len=MAX_CHAR_LENGTH) :: string
 

    !------------------------------------------------------------
    ! 4. Check the consistency of the parameters
    !------------------------------------------------------------
    ! for the time stepping scheme
  
    IF (si_offctr>1._wp.OR.si_offctr<0._wp) THEN
      CALL finish( TRIM(routine), 'Invalid offcentering parameter.'//&
                 & 'Valid range for si_offctr is [0,1].')
    ENDIF
  
    IF (lshallow_water .AND. lsi_3d) THEN
      CALL message( TRIM(routine), 'In shallow water mode lsi_3d is set to .FALSE.')
      lsi_3d=.FALSE.
    ENDIF
  
    IF((itime_scheme<=0).OR.(itime_scheme>=unknown)) THEN
      WRITE(string,'(A,i2)') &
      'wrong value of itime_scheme, must be 1 ...', unknown -1
      CALL finish( TRIM(routine),TRIM(string))
    ENDIF
  
    IF((itime_scheme==tracer_only).AND.(.NOT.ltransport)) THEN
      WRITE(string,'(A,i2,A)') &
      'itime_scheme set to ', tracer_only, 'but ltransport to .FALSE.'
      CALL finish( TRIM(routine),TRIM(string))
    ENDIF
  
    IF(ltransport .AND. ntracer <= 0 .AND. iforcing/=INWP) THEN
      CALL finish( TRIM(routine),'tracer advection not possible for ntracer=0')
      ! [for nwp forcing ntracer setting is treated in setup_transport]
    ENDIF
  
    IF(asselin_coeff<0._wp) THEN
      CALL finish( TRIM(routine),'wrong (negative) coefficient of Asselin filter')
    ENDIF
  
    IF(si_2tls<0.5_wp) THEN
      CALL finish( TRIM(routine),&
                   'off centering parameter si_2tls outside stability range')
    ENDIF
  
    IF(si_2tls>1._wp) THEN
      CALL finish( TRIM(routine),'off centering parameter si_2tls larger than 1')
    ENDIF
  
    IF (global_cell_type/=3) THEN
      idiv_method = 1
    ENDIF
    IF (idiv_method > 2 .OR. idiv_method < 1 )THEN
       CALL finish(TRIM(routine),'Error: idiv_method must be 1 or 2 !')
    ENDIF
  
    SELECT CASE (iequations)
    CASE (IHS_ATM_TEMP, IHS_ATM_THETA)
  
      IF (ldry_dycore) THEN
        CALL message(TRIM(routine),'running the DRY dynamical core')
      ELSE
        CALL message(TRIM(routine),'running the dynamical core with tracer')
      ENDIF
  
      IF (lref_temp) THEN
        CALL message(TRIM(routine), &
             'use of reference temperature switched ON in ' // &
             'calculation of pressure gradient force.')
      ELSE
        CALL message(TRIM(routine), &
             'use of reference temperature switched OFF in ' // &
             'calculation of pressure gradient force.')
      ENDIF
  
    END SELECT


    !-----------------------------------------------------------------------
    ! 6. Assign value to dependent variables 
    !-----------------------------------------------------------------------
    ALLOCATE( nnow(i_ndom),nnew(i_ndom),nold(i_ndom),nsav1(i_ndom), &
            & nsav2(i_ndom),nnow_rcf(i_ndom), nnew_rcf(i_ndom))

    IF (lrestart) THEN
      ! Read time level indices from restart file.
      ! NOTE: this part will be modified later for a proper handling 
      ! of multiple domains!!!

      IF (i_ndom>1) &
      CALL finish(TRIM(routine),'Restart functionality can not handle multiple domains (yet)')

      jdom = 1  ! only consider one domain at the moment
      !DO jdom = 1,i_ndom
        CALL get_restart_attribute( 'nold'    ,nold    (jdom) )
        CALL get_restart_attribute( 'nnow'    ,nnow    (jdom) )
        CALL get_restart_attribute( 'nnew'    ,nnew    (jdom) )
        CALL get_restart_attribute( 'nnow_rcf',nnow_rcf(jdom) )
        CALL get_restart_attribute( 'nnew_rcf',nnew_rcf(jdom) )
      !END DO

    ELSE
      nnow(:) = 1
      nnew(:) = 2
      nold(:) = 3
      nnow_rcf(:) = 1
      nnew_rcf(:) = 2
    END IF

    IF (itime_scheme/=LEAPFROG_EXPL .AND. itime_scheme/=LEAPFROG_SI ) THEN
      ltwotime = .TRUE.
    ELSE
      ltwotime = .FALSE.
    END IF

    ! For two time-level schemes, we also allocate 3 prognostic state vectors:
    ! new, now AND old. The "old" vector is used as an temporary vector,
    ! e.g., for the stages in the Runge-Kutta method.

    nsav1 = 0
    nsav2 = 4

    !-----------------------------------------------------------------------
    ! assign value to variables that are frequently used by
    ! the time stepping scheme
    !-----------------------------------------------------------------------

    gdt         = grav*dtime
    gsi_2tlsdt  = gdt*si_2tls
    gsi_2tls2dt = gsi_2tlsdt*si_2tls*dtime
    g1msi_2tlsdt= (1._wp -si_2tls)*gdt

  END SUBROUTINE dynamics_nml_setup
  !-------------
  !>
  !!
  SUBROUTINE cleanup_dyn_params 

    DEALLOCATE(nnow,nnew,nold,nsav1,nsav2,nnow_rcf,nnew_rcf)

  END SUBROUTINE cleanup_dyn_params 
  !-------------

END MODULE mo_dynamics_nml
