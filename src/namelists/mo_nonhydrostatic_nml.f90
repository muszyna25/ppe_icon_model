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

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: max_char_length, max_dom
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_nh_dyn_config       !ONLY: nh_dyn_config
  USE mo_master_nml,          ONLY: lrestart
  USE mo_mpi,                 ONLY: p_pe, p_io
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE
  PUBLIC

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !-----------------------------------------------------------------------------
  ! Namelist variables
  !-----------------------------------------------------------------------------

  INTEGER :: nml_iadv_rcf              !if 1: no reduced calling frequency for adv. and phy.
                                       !if 2: adv. and phys. are called only every 2nd
                                       !      time step.
                                       !if 4: called every 4th time step ...
  INTEGER :: nml_ivctype               ! Type of vertical coordinate (Gal-Chen / SLEVE)
  REAL(wp):: nml_htop_moist_proc       ! Top height (in m) of the part of the model domain
                                       ! where processes related to moist physics are computed
  INTEGER :: nml_kstart_moist(max_dom) ! related flow control variable (NOT a namelist variable)
  INTEGER :: nml_kstart_qv(max_dom)    ! related flow control variable (NOT a namelist variable)
  REAL(wp):: nml_htop_qvadv            ! Top height (in m) up to which water vapor is advected
                                       ! workaround to circumvent CFL instability in the 
                                       ! stratopause region for aquaplanet experiments
 

  ! Parameters active with cell_type=3 only
  REAL(wp):: nml_damp_height(max_dom)    ! height at which w-damping and sponge layer start
  REAL(wp):: nml_damp_height_u           ! height at which Rayleigh damping of u starts
  REAL(wp):: nml_rayleigh_coeff(max_dom) ! Rayleigh damping coefficient in w-equation
  REAL(wp):: nml_damp_timescale_u ! damping time scale for u in uppermost layer (in seconds)
  REAL(wp):: nml_vwind_offctr   ! Off-centering in vertical wind solver
  INTEGER :: nml_iadv_rhotheta  ! Advection scheme used for density and pot. temperature
  INTEGER :: nml_igradp_method  ! Method for computing the horizontal presure gradient
  REAL(wp):: nml_exner_expol    ! Temporal extrapolation of Exner for computation of
                            ! horizontal pressure gradient
  LOGICAL :: nml_l_open_ubc     ! .true.: open upper boundary condition (w=0 otherwise)
  LOGICAL :: nml_l_nest_rcf     ! .true.: call nests only with rcf frequency
  LOGICAL :: nml_l_masscorr_nest! Apply mass conservation correction also to nested domain
  LOGICAL :: nml_l_zdiffu_t     ! .true.: apply truly horizontal temperature diffusion over steep slopes
  REAL(wp):: nml_thslp_zdiffu   ! threshold slope above which temperature diffusion is applied
  REAL(wp):: nml_thhgtd_zdiffu  ! threshold height difference between adjacent model grid points
                            ! above which temperature diffusion is applied

  ! Parameters active with cell_type=6 only
  REAL(wp) :: nml_gmres_rtol_nh ! relative tolerance for gmres convergence
  LOGICAL  :: nml_ltheta_up_hori! horizontal 3rd order advection of theta_v
  REAL(wp) :: nml_upstr_beta    ! =1 for 3rd order upstream, =0 for 4th order centered
                            ! theta advection
  LOGICAL  :: nml_ltheta_up_vert ! upwind vertical advection of theta
  REAL(wp) :: nml_k2_updamp_coeff ! 2nd order additional horizontal diffusion
                            ! coefficient in the upper damping zone


  NAMELIST/nonhydrostatic_nml/ nml_rayleigh_coeff, nml_damp_height, nml_iadv_rhotheta, &
                               nml_vwind_offctr, nml_igradp_method, nml_exner_expol,   &
                               nml_ltheta_up_hori, nml_ltheta_up_vert, nml_gmres_rtol_nh, &
                               nml_iadv_rcf, nml_ivctype, nml_upstr_beta, nml_l_open_ubc,  &
                               nml_l_nest_rcf, nml_l_zdiffu_t, nml_thslp_zdiffu, &
                               & nml_thhgtd_zdiffu, &
                               nml_k2_updamp_coeff, nml_l_masscorr_nest, nml_htop_moist_proc, &
                               nml_htop_qvadv, nml_damp_timescale_u, nml_damp_height_u

CONTAINS
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


  !local variable
  INTEGER :: i_status, jg

  !------------------------------------------------------------
  ! 2.0 set up the default values for dynamics_nml
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
  damp_height(1)    = 30000.0_wp
  rayleigh_coeff(1) = 0.05_wp
  damp_timescale_u  = 3._wp*86400._wp ! 3 days
  damp_height_u     = 100000._wp
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
  ltheta_up_vert = .FALSE.
  k2_updamp_coeff= 2.0e6_wp
  !
  !
  !------------------------------------------------------------
  ! 3.0 Read the nonhydrostatic namelist.
  !------------------------------------------------------------
  ! (done so far by all MPI processes)

  CALL position_nml ('nonhydrostatic_nml', status=i_status)
  SELECT CASE (i_status)
  CASE (positioned)
     READ (nnml, nonhydrostatic_nml)
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

  IF(p_pe == p_io) WRITE(nnml_output,nml=nonhydrostatic_nml)


END SUBROUTINE nonhydrostatic_nml_setup


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for nonhydrostatic core. 
  !!
  !! This subroutine 
  !! - reads the Namelist for nonhydrostatic core
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-06-07)
  !!
  SUBROUTINE read_nonhydrostatic_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg           ! loop index

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nonhydrostatic_nml: read_nonhydrostatic_namelist'

    !-----------------------
    ! 1. default settings
    !-----------------------

    ! reduced calling frequency for transport
    nml_iadv_rcf = 1  ! no reduced calling frequency

    ! Type of vertical coordinate (1: Gal-Chen, 2: SLEVE)
    nml_ivctype        = 1

    ! Top height of partial domain where moist physics is computed
    ! (set to 200 km, which in practice means that moist physics is
    ! computed everywhere by default)
    nml_htop_moist_proc = 200000._wp
    nml_htop_qvadv      = 250000._wp
    nml_kstart_moist(:) = 1
    nml_kstart_qv(:)    = 1

    ! Settings for icell_type=3
    nml_damp_height(1)    = 30000.0_wp
    nml_rayleigh_coeff(1) = 0.05_wp
    nml_damp_timescale_u  = 3._wp*86400._wp ! 3 days
    nml_damp_height_u     = 100000._wp
    nml_vwind_offctr      = 0.05_wp
    nml_iadv_rhotheta     = 2
    nml_igradp_method     = 1
    nml_l_open_ubc        = .FALSE.
    nml_l_nest_rcf        = .TRUE.
    nml_l_masscorr_nest   = .FALSE.
    nml_exner_expol       = 0.5_wp

    ! dummy values for nested domains; will be reset to value of domain 1 
    ! if not specified explicitly in the namelist
    nml_damp_height(2:max_dom)    = -1.0_wp
    nml_rayleigh_coeff(2:max_dom) = -1.0_wp

    ! truly horizontal temperature diffusion
    nml_l_zdiffu_t     = .FALSE. ! not used by default
    nml_thslp_zdiffu   = 0.025_wp ! slope threshold 0.025
    nml_thhgtd_zdiffu  = 200._wp ! threshold for height difference between adjacent grid points 200 m

    ! Settings for icell_type=6
    nml_ltheta_up_hori =.FALSE.
    nml_gmres_rtol_nh  = 1.0e-6_wp
    nml_upstr_beta     = 1.0_wp
    nml_ltheta_up_vert = .FALSE.
    nml_k2_updamp_coeff= 2.0e6_wp

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('nonhydrostatic_nml')
      READ(funit,NML=nonhydrostatic_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('nonhydrostatic_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, nonhydrostatic_nml)
    END SELECT
    CALL close_nml

    !--------------------------------------------------------------------
    ! Sanity and internal checks
    !--------------------------------------------------------------------
    ! Reset values of rayleigh_coeff and damp_height to that of the 
    ! global domain if not specified

    DO jg = 2, max_dom
      IF (nml_damp_height(jg) < 0.0_wp) THEN
        nml_damp_height(jg) = nml_damp_height(1)
      ENDIF
      IF (nml_rayleigh_coeff(jg) < 0.0_wp) THEN
        nml_rayleigh_coeff(jg) = nml_rayleigh_coeff(1)
      ENDIF
    ENDDO

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
! KF postponed work

!    DO jg = 1,max_dom
!      nh_dyn_config(jg)%rayleigh_coeff  = rayleigh_coeff(jg)
!      nh_dyn_config(jg)%damp_height     = damp_height(jg)
!      nh_dyn_config(jg)%iadv_rhotheta   = iadv_rhotheta
!      nh_dyn_config(jg)%vwind_offctr    = vwind_offctr
!      nh_dyn_config(jg)%igradp_method   = igradp_method
!      nh_dyn_config(jg)%exner_expol     = exner_expol
!      nh_dyn_config(jg)%ltheta_up_hori  = ltheta_up_hori
!      nh_dyn_config(jg)%ltheta_up_vert  = ltheta_up_vert
!      nh_dyn_config(jg)%gmres_rtol_nh   = gmres_rtol_nh
!      nh_dyn_config(jg)%iadv_rcf        = iadv_rcf
!      nh_dyn_config(jg)%ivctype         = ivctype
!      nh_dyn_config(jg)%upstr_beta      = upstr_beta
!      nh_dyn_config(jg)%l_open_ubc      = l_open_ubc
!      nh_dyn_config(jg)%l_nest_rcf      = l_nest_rcf
!      nh_dyn_config(jg)%l_zdiffu_t      = l_zdiffu_t
!      nh_dyn_config(jg)%thslp_zdiffu    = thslp_zdiffu
!      nh_dyn_config(jg)%thhgtd_zdiffu   = thhgtd_zdiffu
!      nh_dyn_config(jg)%k2_updamp_coeff = k2_updamp_coeff
!      nh_dyn_config(jg)%l_masscorr_nest = l_masscorr_nest
!      nh_dyn_config(jg)%htop_moist_proc = htop_moist_proc
!      nh_dyn_config(jg)%htop_qvadv      = htop_qvadv
!      nh_dyn_config(jg)%damp_timescale_u= damp_timescale_u
!      nh_dyn_config(jg)%damp_height_u   = damp_height_u
!    ENDDO
!
       rayleigh_coeff(:)  = nml_rayleigh_coeff(:)
       damp_height   (:)  = nml_damp_height   (:)
       kstart_moist  (:)  = nml_kstart_moist  (:)
       kstart_qv     (:)  = nml_kstart_qv     (:)
       iadv_rhotheta   = nml_iadv_rhotheta
       vwind_offctr    = nml_vwind_offctr
       igradp_method   = nml_igradp_method
       exner_expol     = nml_exner_expol
       ltheta_up_hori  = nml_ltheta_up_hori
       ltheta_up_vert  = nml_ltheta_up_vert
       gmres_rtol_nh   = nml_gmres_rtol_nh
       iadv_rcf        = nml_iadv_rcf
       ivctype         = nml_ivctype
       upstr_beta      = nml_upstr_beta
       l_open_ubc      = nml_l_open_ubc
       l_nest_rcf      = nml_l_nest_rcf
       l_zdiffu_t      = nml_l_zdiffu_t
       thslp_zdiffu    = nml_thslp_zdiffu
       thhgtd_zdiffu   = nml_thhgtd_zdiffu
       k2_updamp_coeff = nml_k2_updamp_coeff
       l_masscorr_nest = nml_l_masscorr_nest
       htop_moist_proc = nml_htop_moist_proc
       htop_qvadv      = nml_htop_qvadv
       damp_timescale_u= nml_damp_timescale_u
       damp_height_u   = nml_damp_height_u
  
    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=nonhydrostatic_nml)                    
    CALL store_and_close_namelist(funit, 'nonhydrostatic_nml') 

    ! 6. write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=nonhydrostatic_nml)

  END SUBROUTINE read_nonhydrostatic_namelist


END MODULE mo_nonhydrostatic_nml
