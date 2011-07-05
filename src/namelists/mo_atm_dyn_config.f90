!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_atm_dyn_config

  USE mo_kind,     ONLY: wp
  USE mo_impl_constants, ONLY: MAX_DOM

  IMPLICIT NONE

  ! All types declarations are private

  PRIVATE :: atm_dyn_config
  PRIVATE :: t_atm_dyn_config
  PRIVATE :: t_ha_dyn_config
  PRIVATE :: t_nh_dyn_config
  PRIVATE :: t_sleve_config
  PRIVATE :: t_hdiff_config

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  ! All other entities, i.e., the "get_xyz" functions, are public.

  PUBLIC  

  !>
  !!--------------------------------------------------------------------------
  !!        Control variables specific to the hydrostatic atm model
  !!--------------------------------------------------------------------------
  TYPE :: t_ha_dyn_config

    INTEGER  :: ileapfrog_startup  !<
    REAL(wp) :: asselin_coeff      !<

    INTEGER  :: si_expl_scheme     !< Scheme for the explicit part of the
                                   !< 2-time-level semi-implicit time integration.
                                   !< See mo_atm_constants for the options.
    REAL(wp) :: si_coeff           !< = 0 : explicit scheme
                                   !< = 1 : semi implicit scheme.
                                   !< in (0,1): a weighted scheme
    REAL(wp) :: si_offctr          !< Weighting parameter used in calculating the
                                   !< second temporal derivatives in the semi-implicit
                                   !< correction scheme. The value read from namelist are
                                   !< assumed to be the offcentering (i.e. between 0 and 1).
    REAL(wp) :: si_cmin            !< Min. phase speed of the decomposed modes to be
                                   !< solved by the semi-implicit correction scheme
    REAL(wp) :: si_rtol            !< Relative tolerance
    LOGICAL  :: lsi_3d             !< If .true., solve the 3D equation

    LOGICAL :: ldry_dycore !< If .TRUE., ignore the effact of water vapor,
                           !< cloud liquid and cloud ice on virtual temperature.
    LOGICAL :: lref_temp   !< If .TRUE., involve the reference temperature profile
                           !< in the calculation of pressure gradient force.

  END TYPE t_ha_dyn_config 

  !>
  !!--------------------------------------------------------------------------
  !!        Control variables specific to the nonhydrostatic atm model
  !!--------------------------------------------------------------------------
  TYPE :: t_nh_dyn_config

    INTEGER :: iadv_rcf        !if 1: no reduced calling frequency for adv. and phy.
                               !if 2: adv. and phys. are called only every 2nd
                               !      time step.
                               !if 4: called every 4th time step ...
    INTEGER :: ivctype         ! Type of vertical coordinate (Gal-Chen / SLEVE)
    REAL(wp):: htop_moist_proc ! Top height (in m) of the part of the model domain
                               ! where processes related to moist physics are computed
    INTEGER :: kstart_moist(max_dom) ! related flow control variable (NOT a namelist variable)
    REAL(wp):: htop_qvadv            ! Top height (in m) up to which water vapor is advected
                                     ! workaround to circumvent CFL instability in the stratopause
                                     ! region for aquaplanet experiments
    INTEGER :: kstart_qv(max_dom)    ! related flow control variable (NOT a namelist variable)
  
    ! Parameters active with i_cell_type=3 only

    REAL(wp):: damp_height(max_dom)    ! height at which w-damping and sponge layer start
    REAL(wp):: damp_height_u           ! height at which Rayleigh damping of u starts
    REAL(wp):: rayleigh_coeff(max_dom) ! Rayleigh damping coefficient in w-equation
    REAL(wp):: damp_timescale_u ! damping time scale for u in uppermost layer (in seconds)
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

    REAL(wp) :: gmres_rtol_nh   ! relative tolerance for gmres convergence
    LOGICAL  :: ltheta_up_hori  ! horizontal 3rd order advection of theta_v
    REAL(wp) :: upstr_beta      ! =1 for 3rd order upstream, =0 for 4th order centered
                                ! theta advection
    LOGICAL :: ltheta_up_vert   ! upwind vertical advection of theta
    REAL(wp) :: k2_updamp_coeff ! 2nd order additional horizontal diffusion
                                ! coefficient in the upper damping zone

  END TYPE t_nh_dyn_config 

  !>
  !!--------------------------------------------------------------------------
  !!                       SLEVE vertical coordinate
  !!--------------------------------------------------------------------------
  TYPE :: t_sleve_config

    ! a) Parameters specifying the distrubution of the coordinate surfaces
    !     (the initializations are a workaround for a NEC compiler bug)

    REAL(wp):: min_lay_thckn = 1._wp  ! Layer thickness of lowermost level
    REAL(wp):: stretch_fac   = 1._wp  ! Factor for stretching/squeezing the model layer distribution
    REAL(wp):: top_height    = 1._wp  ! Height of model top
  
    ! b) Parameters for SLEVE definition

    REAL(wp):: decay_scale_1 = 1._wp  ! Decay scale for large-scale topography component
    REAL(wp):: decay_scale_2 = 1._wp  ! Decay scale for small-scale topography component
    REAL(wp):: decay_exp     = 1._wp  ! Exponent for decay function
    REAL(wp):: flat_height   = 1._wp  ! Height above which the coordinate surfaces are exactly flat
                                      ! additional feature not available in the standard
                                      ! SLEVE definition
  END TYPE t_sleve_config
  !>
  !!--------------------------------------------------------------------------
  !!                           Horizontal Diffusion
  !!--------------------------------------------------------------------------
  TYPE :: t_hdiff_config

    LOGICAL :: lhdiff_temp  ! if .TRUE., apply horizontal diffusion to thermodynamic variable.
    LOGICAL :: lhdiff_vn    ! if .TRUE., apply horizontal diffusion to momentum.

    INTEGER :: hdiff_order  ! order of horizontal diffusion
                            ! 2: 2nd order linear diffusion on all vertical levels 
                            ! 3: Smagorinsky diffusion for hexagonal model
                            ! 4: 4th order linear diffusion on all vertical levels 
                            ! 5: Smagorinsky diffusion for triangular model
                            ! 24 or 42: 2nd order linear diffusion for upper levels,
                            !           4th order for lower levels
  
    REAL(wp) :: k2_pres_max  ! (relevant only when hdiff_order = 24 or 42)
                             ! pressure (in Pa) specified by the user
                             ! to determine the lowest vertical level 
                             ! to which 2nd order linear diffusion is applied.
                             ! For the levels with pressure > k2_pres_max, 
                             ! 4th order linear diffusion is applied. 
  
    INTEGER :: k2_klev_max  ! (relevant only when hdiff_order = 24 or 42)
                            ! vertical level index specified by the user
                            ! to determine the lowest vertical level 
                            ! to which 2nd order linear diffusion is applied.
                            ! For the levels with k > k2_klev_max, 
                            ! 4th order linear diffusion is applied. 

    REAL(wp) ::               &
      & hdiff_efdt_ratio,     &! ratio of e-folding time to (2*)time step
      & hdiff_min_efdt_ratio, &! minimum value of hdiff_efdt_ratio (for upper sponge layer)
      & hdiff_tv_ratio,       &! the ratio of diffusion coefficient: temp:mom
      & hdiff_smag_fac,       &! scaling factor for Smagorinsky diffusion
      & hdiff_multfac          ! multiplication factor of normalized diffusion coefficient
                               ! for nested domains
    REAL(wp), ALLOCATABLE :: &
      & k6(:), k4(:), k2(:)  ! numerical diffusion coefficients
                             ! Values for these parameters are not directly
                             ! specified by the user, but derived from the ratio 
                             ! between the e-folding time and the model time step
                             ! (hdiff_efdt_ratio above), and the horizontal 
                             ! resolution of the model
  
    INTEGER k2s, k2e, k4s, k4e  ! indices defining to which vertical levels
                                ! 2nd and 4th linear diffusion are applied.
                                ! The values are not specified by the user via namelist,
                                ! but determined from k2_klev_max, k2_pres_max
                                ! and the configuration of the vertical coordinate
  
  END TYPE t_hdiff_config
  !>
  !!--------------------------------------------------------------------------
  !!          Top-level configuration state for atm dynamics
  !!--------------------------------------------------------------------------
  TYPE :: t_atm_dyn_config

    ! namelist variables

    INTEGER :: ieqautions      !< Choice of governing equation set
    INTEGER :: itime_scheme    !< Choice of time stepping scheme
    INTEGER :: i_cell_type     !< Shape of control volume. 3 = triangle, 6 = hexagon/pentagon
    INTEGER :: idiv_method     !< Divergence operator
    INTEGER :: divavg_cntrwgt  !< Weight of central cell for divergence averaging

    ! derived variables

    LOGICAL :: ltwotime
    INTEGER,ALLOCATABLE :: nold(:), nnow(:), nnew(:) ! variables denoting time levels
    INTEGER,ALLOCATABLE :: nsav1(:), nsav2(:)        ! extra 'time levels' of prognostic variables
                                                     ! needed to compute boundary tendencies and
                                                     ! feedback increments
    INTEGER,ALLOCATABLE :: nnow_rcf(:), nnew_rcf(:)  ! extra time levels for reduced
                                                     ! calling frequency (rcf)

    TYPE(t_ha_dyn_config) :: ha
    TYPE(t_nh_dyn_config) :: nh
    TYPE(t_sleve_config)  :: sleve
    TYPE(t_hdiff_config)  :: hdiff

  END TYPE t_atm_dyn_config

  TYPE(t_atm_dyn_config),ALLOCATABLE :: atm_dyn_config(:)  !< the variable

END MODULE mo_atm_dyn_config

