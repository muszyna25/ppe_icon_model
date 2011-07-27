!>
!! Contains the setup of configuration of the
!! nonhydrostatic dynamical core
!!        
!! @par Revision History
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

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish
  USE mo_impl_constants,        ONLY: max_dom
  USE mo_io_units,              ONLY: nnml, nnml_output
  USE mo_namelist,              ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,        ONLY: is_restart_run
  USE mo_mpi,                   ONLY: my_process_is_stdio
  USE mo_io_restart_namelist,   ONLY: open_tmpfile, store_and_close_namelist,  &
                                    & open_and_restore_namelist, close_tmpfile

  USE mo_nonhydrostatic_config, ONLY: &
                                    ! from namelist
                                    & config_itime_scheme     => itime_scheme     , &
                                    & config_iadv_rcf         => iadv_rcf         , &
                                    & config_ivctype          => ivctype          , &
                                    & config_htop_moist_proc  => htop_moist_proc  , &
                                    & config_htop_qvadv       => htop_qvadv       , &
                                    & config_damp_height      => damp_height      , &
                                    & config_damp_height_u    => damp_height_u    , &
                                    & config_rayleigh_coeff   => rayleigh_coeff   , &
                                    & config_damp_timescale_u => damp_timescale_u , &
                                    & config_vwind_offctr     => vwind_offctr     , &
                                    & config_iadv_rhotheta    => iadv_rhotheta    , &
                                    & config_igradp_method    => igradp_method    , &
                                    & config_exner_expol      => exner_expol      , &
                                    & config_l_open_ubc       => l_open_ubc       , &
                                    & config_l_nest_rcf       => l_nest_rcf       , &
                                    & config_l_masscorr_nest  => l_masscorr_nest  , &
                                    & config_l_zdiffu_t       => l_zdiffu_t       , &
                                    & config_thslp_zdiffu     => thslp_zdiffu     , &
                                    & config_thhgtd_zdiffu    => thhgtd_zdiffu    , &
                                    & config_gmres_rtol_nh    => gmres_rtol_nh    , &
                                    & config_ltheta_up_hori   => ltheta_up_hori   , &
                                    & config_upstr_beta       => upstr_beta       , &
                                    & config_ltheta_up_vert   => ltheta_up_vert   , &
                                    & config_k2_updamp_coeff  => k2_updamp_coeff  , &
                                    ! not from namelist
                                    & config_kstart_moist     => kstart_moist     , &
                                    & config_kstart_qv        => kstart_qv

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: read_nonhydrostatic_namelist

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !-----------------------------------------------------------------------------
  ! Namelist variables
  !-----------------------------------------------------------------------------


  INTEGER  :: itime_scheme   ! parameter used to select the time stepping scheme
                             ! = 1, explicit 2 time level scheme for tracer
                             ! = 3, Matsuno, comp of velocity tendencies on corretor step only
                             ! = 4, Matsuno scheme
                             ! = 5, under development
                             ! = 6, Matsuno with velocitiy tendendcies averaged over 2 time steps

  INTEGER :: iadv_rcf                ! if 1: no reduced calling frequency for adv. and phy.
                                     ! if 2: adv. and phys. are called every 2nd time step.
                                     ! if 4: ... every 4th time step.
  INTEGER :: ivctype                 ! Type of vertical coordinate (Gal-Chen / SLEVE)
  REAL(wp):: htop_moist_proc         ! Top height (in m) of the part of the model domain
                                     ! where processes related to moist physics are computed
  REAL(wp):: htop_qvadv              ! Top height (in m) up to which water vapor is advected
                                     ! workaround to circumvent CFL instability in the 
                                     ! stratopause region for aquaplanet experiments

  ! Parameters active with cell_type=3 only
  REAL(wp):: damp_height(max_dom)    ! height at which w-damping and sponge layer start
  REAL(wp):: damp_height_u           ! height at which Rayleigh damping of u starts
  REAL(wp):: rayleigh_coeff(max_dom) ! Rayleigh damping coefficient in w-equation
  REAL(wp):: damp_timescale_u        ! damping time scale for u in uppermost layer (in seconds)
  REAL(wp):: vwind_offctr            ! Off-centering in vertical wind solver
  INTEGER :: iadv_rhotheta           ! Advection scheme used for density and pot. temperature
  INTEGER :: igradp_method           ! Method for computing the horizontal presure gradient
  REAL(wp):: exner_expol             ! Temporal extrapolation of Exner for computation of
                                     ! horizontal pressure gradient
  LOGICAL :: l_open_ubc              ! .true.: open upper boundary condition (w=0 otherwise)
  LOGICAL :: l_nest_rcf              ! .true.: call nests only with rcf frequency
  LOGICAL :: l_masscorr_nest         ! Apply mass conservation correction also to nested domain
  LOGICAL :: l_zdiffu_t              ! .true.: apply truly horizontal temperature diffusion
                                     !         over steep slopes
  REAL(wp):: thslp_zdiffu            ! threshold slope above which temperature diffusion is applied
  REAL(wp):: thhgtd_zdiffu           ! threshold height diff. between adjacent model grid points
                                     ! above which temperature diffusion is applied

  ! Parameters active with cell_type=6 only
  REAL(wp) :: gmres_rtol_nh          ! relative tolerance for gmres convergence
  LOGICAL  :: ltheta_up_hori         ! horizontal 3rd order advection of theta_v
  REAL(wp) :: upstr_beta             ! =1 for 3rd order upstream, =0 for 4th order centered
                                     ! theta advection
  LOGICAL  :: ltheta_up_vert         ! upwind vertical advection of theta
  REAL(wp) :: k2_updamp_coeff        ! 2nd order additional horizontal diffusion
                                     ! coefficient in the upper damping zone

  ! Reated parameters not part of the namelist
  INTEGER :: kstart_moist(max_dom)   ! related flow control variable (NOT a namelist variable)
  INTEGER :: kstart_qv(max_dom)      ! related flow control variable (NOT a namelist variable)


  NAMELIST /nonhydrostatic_nml/ itime_scheme, iadv_rcf, ivctype, htop_moist_proc,          &
                              & htop_qvadv, damp_height, damp_height_u, rayleigh_coeff,    &
                              & damp_timescale_u, vwind_offctr, iadv_rhotheta,             &
                              & igradp_method, exner_expol, l_open_ubc, l_nest_rcf,        &
                              & l_masscorr_nest, l_zdiffu_t, thslp_zdiffu, thhgtd_zdiffu,  &
                              & gmres_rtol_nh, ltheta_up_hori, upstr_beta, ltheta_up_vert, &
                              & k2_updamp_coeff

CONTAINS
  !-------------------------------------------------------------------------
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
  !!  by Daniel Reinert, DWD (2011-07-06)
  !!
  SUBROUTINE read_nonhydrostatic_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg           ! loop index

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nonhydrostatic_nml:read_nonhydrostatic_namelist'

    !-----------------------
    ! 1. default settings
    !-----------------------

    ! Time scheme chose for the nh-model
    itime_scheme = 4 ! Matsuno scheme

    ! reduced calling frequency for transport
    iadv_rcf = 1  ! no reduced calling frequency

    ! Type of vertical coordinate (1: Gal-Chen, 2: SLEVE)
    ivctype  = 1

    ! Top height of partial domain where moist physics is computed
    ! (set to 200 km, which in practice means that moist physics is
    ! computed everywhere by default)
    htop_moist_proc = 200000._wp
    htop_qvadv      = 250000._wp
    kstart_moist(:) = 1
    kstart_qv(:)    = 1

    ! Settings for icell_type=3
    damp_height(1)    = 30000.0_wp
    damp_height_u     = 100000._wp
    rayleigh_coeff(1) = 0.05_wp
    damp_timescale_u  = 3._wp*86400._wp ! 3 days
    vwind_offctr      = 0.05_wp
    iadv_rhotheta     = 2
    igradp_method     = 1
    exner_expol       = 0.5_wp
    l_open_ubc        = .FALSE.
    l_nest_rcf        = .TRUE.
    l_masscorr_nest   = .FALSE.

    ! dummy values for nested domains; will be reset to value of domain 1 
    ! if not specified explicitly in the namelist
    damp_height(2:max_dom)    = -1.0_wp
    rayleigh_coeff(2:max_dom) = -1.0_wp

    ! truly horizontal temperature diffusion
    l_zdiffu_t     = .FALSE.  ! not used by default
    thslp_zdiffu   = 0.025_wp ! slope threshold 0.025
    thhgtd_zdiffu  = 200._wp  ! threshold for height difference between adjacent grid points 200 m

    ! Settings for icell_type=6
    gmres_rtol_nh  = 1.0e-6_wp
    ltheta_up_hori =.FALSE.
    upstr_beta     = 1.0_wp
    ltheta_up_vert = .FALSE.
    k2_updamp_coeff= 2.0e6_wp

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
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
      IF (damp_height(jg) < 0.0_wp) THEN
        damp_height(jg) = damp_height(1)
      ENDIF
      IF (rayleigh_coeff(jg) < 0.0_wp) THEN
        rayleigh_coeff(jg) = rayleigh_coeff(1)
      ENDIF
    ENDDO

    ! reset l_nest_rcf to false if iadv_rcf = 1
    IF (iadv_rcf == 1) l_nest_rcf = .FALSE.

    IF (upstr_beta > 1.0_wp .OR. upstr_beta < 0.0_wp) THEN
      CALL finish(TRIM(routine), 'upstr_beta out of range 0..1')
    ENDIF

    ! for reduced calling frequency of tracer advection / fast physics:
    ! odd values of iadv_rcf are allowed only if nest calls are 
    ! synchronized with advection
    IF ( .NOT. l_nest_rcf .AND. MOD(iadv_rcf,2) /= 0 .AND. iadv_rcf /= 1 &
         .OR. iadv_rcf == 0) THEN
        CALL finish( TRIM(routine), 'Invalid reduced-calling-frequency parameter& '//&
          &'Value must be even or 1 if l_nest_rcf=.FALSE.')
    ENDIF

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
! KF postponed work

!    DO jg = 1,max_dom
!      nonhydrostatic_config(jg)%rayleigh_coeff  = rayleigh_coeff(jg)
!      nonhydrostatic_config(jg)%damp_height     = damp_height(jg)
!      nonhydrostatic_config(jg)%iadv_rhotheta   = iadv_rhotheta
!      nonhydrostatic_config(jg)%vwind_offctr    = vwind_offctr
!      nonhydrostatic_config(jg)%igradp_method   = igradp_method
!      nonhydrostatic_config(jg)%exner_expol     = exner_expol
!      nonhydrostatic_config(jg)%ltheta_up_hori  = ltheta_up_hori
!      nonhydrostatic_config(jg)%ltheta_up_vert  = ltheta_up_vert
!      nonhydrostatic_config(jg)%gmres_rtol_nh   = gmres_rtol_nh
!      nonhydrostatic_config(jg)%iadv_rcf        = iadv_rcf
!      nonhydrostatic_config(jg)%ivctype         = ivctype
!      nonhydrostatic_config(jg)%upstr_beta      = upstr_beta
!      nonhydrostatic_config(jg)%l_open_ubc      = l_open_ubc
!      nonhydrostatic_config(jg)%l_nest_rcf      = l_nest_rcf
!      nonhydrostatic_config(jg)%l_zdiffu_t      = l_zdiffu_t
!      nonhydrostatic_config(jg)%thslp_zdiffu    = thslp_zdiffu
!      nonhydrostatic_config(jg)%thhgtd_zdiffu   = thhgtd_zdiffu
!      nonhydrostatic_config(jg)%k2_updamp_coeff = k2_updamp_coeff
!      nonhydrostatic_config(jg)%l_masscorr_nest = l_masscorr_nest
!      nonhydrostatic_config(jg)%htop_moist_proc = htop_moist_proc
!      nonhydrostatic_config(jg)%htop_qvadv      = htop_qvadv
!      nonhydrostatic_config(jg)%damp_timescale_u= damp_timescale_u
!      nonhydrostatic_config(jg)%damp_height_u   = damp_height_u
!    ENDDO
!
       config_rayleigh_coeff(:) = rayleigh_coeff(:)
       config_damp_height   (:) = damp_height   (:)
       config_kstart_moist  (:) = kstart_moist  (:)
       config_kstart_qv     (:) = kstart_qv     (:)
       config_iadv_rhotheta     = iadv_rhotheta
       config_vwind_offctr      = vwind_offctr
       config_igradp_method     = igradp_method
       config_exner_expol       = exner_expol
       config_ltheta_up_hori    = ltheta_up_hori
       config_ltheta_up_vert    = ltheta_up_vert
       config_gmres_rtol_nh     = gmres_rtol_nh
       config_iadv_rcf          = iadv_rcf
       config_itime_scheme      = itime_scheme
       config_ivctype           = ivctype
       config_upstr_beta        = upstr_beta
       config_l_open_ubc        = l_open_ubc
       config_l_nest_rcf        = l_nest_rcf
       config_l_zdiffu_t        = l_zdiffu_t
       config_thslp_zdiffu      = thslp_zdiffu
       config_thhgtd_zdiffu     = thhgtd_zdiffu
       config_k2_updamp_coeff   = k2_updamp_coeff
       config_l_masscorr_nest   = l_masscorr_nest
       config_htop_moist_proc   = htop_moist_proc
       config_htop_qvadv        = htop_qvadv
       config_damp_timescale_u  = damp_timescale_u
       config_damp_height_u     = damp_height_u
  
    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=nonhydrostatic_nml)                    
    CALL store_and_close_namelist(funit, 'nonhydrostatic_nml') 

    ! 6. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=nonhydrostatic_nml)

  END SUBROUTINE read_nonhydrostatic_namelist


END MODULE mo_nonhydrostatic_nml
