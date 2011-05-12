!>
!! Contains the setup of configuration of the
!! nonhydrostatic dynamical core
!!
!!        
!! @par Revision History
!!   Revision History in mo_global_variables.f90 (r3914)
!!   Modification by Constantin Junk (2011-03-28)
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
MODULE mo_nonhydrostatic_nml
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
!
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: max_char_length, max_dom
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned
  USE mo_mpi,                ONLY: p_pe, p_io

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC

  CHARACTER(len=*), PARAMETER :: modelname    = 'icon'
  CHARACTER(len=*), PARAMETER :: modelversion = 'dev'

  ! ----------------------------------------------------------------------------
  ! 1.0 Namelist variables for the nonhydrostatic core
  ! ----------------------------------------------------------------------------
  !
  INTEGER :: iadv_rcf       !if 1: no reduced calling frequency for adv. and phy.
                            !if 2: adv. and phys. are called only every 2nd
                            !      time step.
                            !if 4: called every 4th time step ...
  INTEGER :: ivctype        ! Type of vertical coordinate (Gal-Chen / SLEVE)
  REAL(wp):: htop_moist_proc ! Top height (in m) of the part of the model domain
                            ! where processes related to moist physics are computed
  INTEGER :: kstart_moist(max_dom) ! related flow control variable (NOT a namelist variable)
  REAL(wp):: htop_qvadv     ! Top height (in m) up to which water vapor is advected
                            ! workaround to circumvent CFL instability in the stratopause
                            ! region for aquaplanet experiments
  INTEGER :: kstart_qv(max_dom) ! related flow control variable (NOT a namelist variable)

  ! Parameters active with i_cell_type=3 only
  REAL(wp):: damp_height(max_dom)    ! height at which damping starts
  REAL(wp):: rayleigh_coeff(max_dom) ! Rayleigh damping coefficient in w-equation
  REAL(wp):: vwind_offctr   ! Off-centering in vertical wind solver
  INTEGER :: iadv_rhotheta  ! Advection scheme used for density and pot. temperature
  INTEGER :: igradp_method  ! Method for computing the horizontal presure gradient
  REAL(wp):: exner_expol    ! Temporal extrapolation of Exner for computation of
                            ! horizontal pressure gradient
  LOGICAL :: l_open_ubc     ! .true.: open upper boundary condition (w=0 otherwise)
  LOGICAL :: l_nest_rcf     ! .true.: call nests only with rcf frequency
  LOGICAL :: l_masscorr_nest! Apply mass conservation correction also to nested domain
  LOGICAL :: l_zdiffu_t     ! .true.: apply truly horizontal temperature diffusion over steep slopes
  REAL(wp):: thslp_zdiffu   ! threshold slope above which temperature diffusion is applied
  REAL(wp):: thhgtd_zdiffu  ! threshold height difference between adjacent model grid points
                            ! above which temperature diffusion is applied

  ! Parameters active with i_cell_type=6 only
  REAL(wp) :: gmres_rtol_nh ! relative tolerance for gmres convergence
  LOGICAL  :: ltheta_up_hori! horizontal 3rd order advection of theta_v
  REAL(wp) :: upstr_beta    ! =1 for 3rd order upstream, =0 for 4th order centered
                            ! theta advection
  LOGICAL :: l_impl_vert_adv! implicit vertical advection for horizontal velocity
                            ! Note: Implicit vertical advections for rho and w are
                            ! boiled into the 5band matrix vertical solver 
  REAL(wp) :: k2_updamp_coeff ! 2nd order additional horizontal diffusion
                            ! coefficient in the upper damping zone


  NAMELIST/nonhydrostatic_ctl/ rayleigh_coeff, damp_height, iadv_rhotheta, &
                               vwind_offctr, igradp_method, exner_expol,   &
                               ltheta_up_hori, l_impl_vert_adv, gmres_rtol_nh, &
                               iadv_rcf, ivctype, upstr_beta, l_open_ubc,  &
                               l_nest_rcf, l_zdiffu_t, thslp_zdiffu, thhgtd_zdiffu, &
                               k2_updamp_coeff, l_masscorr_nest, htop_moist_proc, &
                               htop_qvadv
  !
  CONTAINS



!-------------------------------------------------------------------------
!
!
 !!
 !>
 !!               Initialization of variables that determine 
 !!
 !!               Initialization of variables that determine
 !!               some settings of the non-hydrostatic dynamical core.
 !!
 !! @par Revision History
 !!  Initial version by Almut Gassmann (2009-03-04)
 !!
 SUBROUTINE nonhydrostatic_nml_setup

  CHARACTER(len=max_char_length), PARAMETER :: &
            routine = '(mo_nonhydrostatic_nml/nonhydrostatic_nml_setup:'
  CHARACTER(len=max_char_length) :: string


  !local variable
  INTEGER :: i_status, jg

  !------------------------------------------------------------
  ! 2.0 set up the default values for dynamics_ctl
  !------------------------------------------------------------
  !
  ! reduced calling frequency for transport
  iadv_rcf = 1  ! no reduced calling frequency

  ! Type of vertical coordinate (1: Gal-Chen, 2: SLEVE)
  ivctype        = 1

  ! Top height of partial domain where moist physics is computed
  ! (set to 200 km, which in practice means that moist physics is
  ! computed everywhere by default)
  htop_moist_proc = 200000._wp
  htop_qvadv      = 250000._wp
  kstart_moist(:) = 1
  kstart_qv(:)    = 1

  ! Settings for icell_type=3
  damp_height(1)    = 17500.0_wp
  rayleigh_coeff(1) = 0.05_wp
  vwind_offctr      = 0.05_wp
  iadv_rhotheta     = 2
  igradp_method     = 1
  l_open_ubc        = .FALSE.
  l_nest_rcf        = .TRUE.
  l_masscorr_nest   = .FALSE.
  exner_expol       = 0.5_wp

  ! dummy values for nested domains; will be reset to value of domain 1 
  ! if not specified explicitly in the namelist
  damp_height(2:max_dom)    = -1.0_wp
  rayleigh_coeff(2:max_dom) = -1.0_wp

  ! truly horizontal temperature diffusion
  l_zdiffu_t     = .FALSE. ! not used by default
  thslp_zdiffu   = 0.025_wp ! slope threshold 0.025
  thhgtd_zdiffu  = 200._wp ! threshold for height difference between adjacent grid points 200 m

  ! Settings for icell_type=6
  ltheta_up_hori =.FALSE.
  gmres_rtol_nh  = 1.0e-6_wp
  upstr_beta     = 1.0_wp
  l_impl_vert_adv=.TRUE.
  k2_updamp_coeff= 2.0e6_wp
  !
  !
  !------------------------------------------------------------
  ! 3.0 Read the nonhydrostatic namelist.
  !------------------------------------------------------------
  ! (done so far by all MPI processes)

  CALL position_nml ('nonhydrostatic_ctl', status=i_status)
  SELECT CASE (i_status)
  CASE (positioned)
     READ (nnml, nonhydrostatic_ctl)
  END SELECT

  ! Reset values of rayleigh_coeff and damp_height to that of the global domain if not specified
  DO jg = 2, max_dom
    IF (damp_height(jg) < 0.0_wp) THEN
      damp_height(jg) = damp_height(1)
    ENDIF
    IF (rayleigh_coeff(jg) < 0.0_wp) THEN
      rayleigh_coeff(jg) = rayleigh_coeff(1)
    ENDIF
  ENDDO

  !------------------------------------------------------------
  ! 4.0 check the consistency of the parameters
  !------------------------------------------------------------
  !

  ! reset l_nest_rcf to false if iadv_rcf = 1
  IF (iadv_rcf == 1) l_nest_rcf = .FALSE.

  IF (upstr_beta > 1.0_wp .OR. upstr_beta < 0.0_wp) THEN
    CALL finish(TRIM(routine), 'upstr_beta out of range 0..1')
  ENDIF

  ! for reduced calling frequency of tracer advection / fast physics:
  ! odd values of iadv_rcf are allowed only if nest calls are synchronized with advection
  IF ( .NOT. l_nest_rcf .AND. MOD(iadv_rcf,2) /= 0 .AND. iadv_rcf /= 1 .OR. iadv_rcf == 0) THEN
      CALL finish( TRIM(routine), 'Invalid reduced-calling-frequency parameter& '//&
        &'Value must be even or 1 if l_nest_rcf=.FALSE.')
  ENDIF

  ! write the contents of the namelist to an ASCII file

  IF(p_pe == p_io) WRITE(nnml_output,nml=nonhydrostatic_ctl)


END SUBROUTINE nonhydrostatic_nml_setup

END MODULE mo_nonhydrostatic_nml
