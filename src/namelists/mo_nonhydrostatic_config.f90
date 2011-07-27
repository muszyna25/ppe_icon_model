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
MODULE mo_nonhydrostatic_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_dom

  IMPLICIT NONE

  PUBLIC


  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !>
  !!----------------------------------------------------------------------------
  !! Derived type containing control variables specific to the nonhydrostatic 
  !! atm model
  !!----------------------------------------------------------------------------
!  TYPE :: t_nonhydrostatic_config

    INTEGER  :: itime_scheme            !< Choice of time stepping scheme

    INTEGER :: iadv_rcf                 !if 1: no reduced calling frequency for adv. and phy.
                                        !if 2: adv. and phys. are called only every 2nd
                                        !      time step.
                                        !if 4: called every 4th time step ...
    INTEGER :: ivctype                  ! Type of vertical coordinate (Gal-Chen / SLEVE)
    REAL(wp):: htop_moist_proc          ! Top height (in m) of the part of the model domain
                                        ! where processes related to moist physics are computed
    INTEGER :: kstart_moist(max_dom)    ! related flow control variable (NOT a namelist variable)
    INTEGER :: kstart_qv(max_dom)       ! related flow control variable (NOT a namelist variable)
    REAL(wp):: htop_qvadv               ! Top height (in m) up to which water vapor is advected
                                        ! workaround to circumvent CFL instability in the 
                                        ! stratopause region for aquaplanet experiments
  
    ! Parameters active with cell_type=3 only

    REAL(wp):: damp_height(max_dom)    ! height at which w-damping and sponge layer start
    REAL(wp):: damp_height_u           ! height at which Rayleigh damping of u starts
    REAL(wp):: rayleigh_coeff(max_dom) ! Rayleigh damping coefficient in w-equation
    REAL(wp):: damp_timescale_u ! damping time scale for u in uppermost layer (in seconds)
    REAL(wp):: vwind_offctr     ! Off-centering in vertical wind solver
    INTEGER :: iadv_rhotheta    ! Advection scheme used for density and pot. temperature
    INTEGER :: igradp_method    ! Method for computing the horizontal presure gradient
    REAL(wp):: exner_expol      ! Temporal extrapolation of Exner for computation of
                                ! horizontal pressure gradient
    LOGICAL :: l_open_ubc       ! .true.: open upper boundary condition (w=0 otherwise)
    LOGICAL :: l_nest_rcf       ! .true.: call nests only with rcf frequency
    LOGICAL :: l_masscorr_nest  ! Apply mass conservation correction also to nested domain
    LOGICAL :: l_zdiffu_t       ! .true.: apply truly horizontal temperature diffusion 
                                ! over steep slopes
    REAL(wp):: thslp_zdiffu     ! threshold slope above which temperature diffusion is applied
    REAL(wp):: thhgtd_zdiffu    ! threshold height difference between adjacent model grid points
                                ! above which temperature diffusion is applied
  
    ! Parameters active with cell_type=6 only

    REAL(wp) :: gmres_rtol_nh   ! relative tolerance for gmres convergence
    LOGICAL  :: ltheta_up_hori  ! horizontal 3rd order advection of theta_v
    REAL(wp) :: upstr_beta      ! =1 for 3rd order upstream, =0 for 4th order centered
                                ! theta advection
    LOGICAL :: ltheta_up_vert   ! upwind vertical advection of theta
    REAL(wp) :: k2_updamp_coeff ! 2nd order additional horizontal diffusion
                                ! coefficient in the upp‚per damping zone

!  END TYPE t_nonhydrostatic_config 
  !>
  !!
!  TYPE(t_nonhydrostatic_config) :: nonhydrostatic_config(max_dom) ! config state 


END MODULE mo_nonhydrostatic_config
