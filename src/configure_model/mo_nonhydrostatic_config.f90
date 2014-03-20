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

  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: max_dom, MAX_CHAR_LENGTH !, & 
!    &                                MATSUNO_DEF,MATSUNO_COR,MATSUNO_AVE  
  USE mo_exception,            ONLY: message, message_text
  USE mo_vertical_coord_table, ONLY: vct_a, vct_b

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
    LOGICAL :: lvadv_tke                ! if true. perform vertical   advection of TKE

    LOGICAL :: lhdiff_rcf               ! if true: compute horizontal diffusion also at the large time step
    LOGICAL :: lextra_diffu             ! if true: apply additional diffusion at grid points close 
                                        ! to the CFL stability limit for vertical advection
    LOGICAL :: lbackward_integr         ! if true: integrate backward in time (needed for testing DFI)
    REAL(wp):: divdamp_fac              ! Scaling factor for divergence damping (used if lhdiff_rcf = true)
    REAL(wp):: divdamp_fac_o2           ! Scaling factor for second-order divergence damping 
                                        ! (derived variable; used if lhdiff_rcf = true and divdamp_order = 2 or 24)
    INTEGER :: divdamp_order            ! Order of divergence damping
    INTEGER :: divdamp_type             ! Type of divergence damping (2D or 3D divergence)
    INTEGER :: ivctype                  ! Type of vertical coordinate (Gal-Chen / SLEVE)
    REAL(wp):: htop_moist_proc          ! Top height (in m) of the part of the model domain
                                        ! where processes related to moist physics are computed
    INTEGER :: kstart_moist(max_dom)    ! related flow control variable (NOT a namelist variable)
    REAL(wp):: hbot_qvsubstep           ! Bottom height (in m) down to which water vapor is 
                                        ! advected with internal substepping (to circumvent CFL 
                                        ! instability in the stratopause region).
    INTEGER :: kend_qvsubstep(max_dom)  ! related flow control variable (NOT a namelist variable)
    INTEGER :: ih_clch(max_dom)         ! end index for levels contributing to high-level clouds, clch
    INTEGER :: ih_clcm(max_dom)         ! end index for levels contributing to mid-level clouds, clcm



    REAL(wp):: rayleigh_type    ! type of Rayleigh damping (1: CLASSIC, 2: Klemp (2008))
    REAL(wp):: damp_height(max_dom)    ! height at which w-damping and sponge layer start
    REAL(wp):: rayleigh_coeff(max_dom) ! Rayleigh damping coefficient in w-equation
    REAL(wp):: vwind_offctr     ! Off-centering in vertical wind solver
    REAL(wp):: rhotheta_offctr  ! Off-centering for density and potential temperature at interface levels
    REAL(wp):: veladv_offctr    ! Off-centering for velocity advection
    INTEGER :: iadv_rhotheta    ! Advection scheme used for density and pot. temperature
    INTEGER :: igradp_method    ! Method for computing the horizontal presure gradient
    REAL(wp):: exner_expol      ! Temporal extrapolation of Exner for computation of
                                ! horizontal pressure gradient
    LOGICAL :: l_open_ubc       ! .true.: open upper boundary condition (w=0 otherwise)
    LOGICAL :: l_nest_rcf       ! .true.: call nests only with rcf frequency
    INTEGER :: nest_substeps           ! the number of dynamics substeps for the child patches
    LOGICAL :: l_masscorr_nest  ! Apply mass conservation correction also to nested domain
    LOGICAL :: l_zdiffu_t       ! .true.: apply truly horizontal temperature diffusion 
                                ! over steep slopes
    REAL(wp):: thslp_zdiffu     ! threshold slope above which temperature diffusion is applied
    REAL(wp):: thhgtd_zdiffu    ! threshold height difference between adjacent model grid points
                                ! above which temperature diffusion is applied


!  END TYPE t_nonhydrostatic_config 
  !>
  !!
!  TYPE(t_nonhydrostatic_config) :: nonhydrostatic_config(max_dom) ! config state 


CONTAINS

  !>
  !! Setup of additional nonhydrostatic control variables
  !!
  !! Setup of additional nonhydrostatic control variables depending on the 
  !! nonhydrostatic-NAMELIST and potentially other namelists. This routine is 
  !! called, after all namelists have been read and a synoptic consistency 
  !! check has been done.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-08-04)
  !!
  SUBROUTINE configure_nonhydrostatic(jg, nlev, nshift_total)
  !
    INTEGER,  INTENT(IN) :: jg           !< patch 
    INTEGER,  INTENT(IN) :: nlev         !< number of full vertical levels 
    INTEGER,  INTENT(IN) :: nshift_total 

    INTEGER :: jk, jk1

    REAL(wp), PARAMETER :: hbase_clch = 7185.44_wp  ! height in m of 400 hPa level in US standard atmosphere
    REAL(wp), PARAMETER :: hbase_clcm = 1948.99_wp  ! height in m of 800 hPa level in US standard atmosphere


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nonhydrostatic_config:configure_nonhydrostatic'

    !-----------------------------------------------------------------------

    ! Determine start level for moist physics processes (specified by htop_moist_proc)
    kstart_moist(jg) = 1
    DO jk = 1, nlev
      jk1 = jk + nshift_total
      IF (0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < htop_moist_proc) THEN
        kstart_moist(jg) = jk
        EXIT
      ENDIF
    ENDDO

    IF ( kstart_moist(jg) >= 1 ) THEN
      WRITE(message_text,'(2(a,i4))') 'Domain', jg, &
        '; computation of moist physics processes starts in layer ', kstart_moist(jg)
      CALL message(TRIM(routine),message_text)
    ENDIF


    ! Determine end level for qv-advection substepping (specified by hbot_qvsubstep)
    kend_qvsubstep(jg) = 0
    DO jk = nlev, 1, -1
      jk1 = jk + nshift_total
      IF (0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) > hbot_qvsubstep) THEN
        kend_qvsubstep(jg) = jk
        EXIT
      ENDIF
    ENDDO

    IF ( kend_qvsubstep(jg) >= 1 ) THEN
      WRITE(message_text,'(2(a,i4))') 'Domain', jg, &
        '; QV substepping ends in layer ', kend_qvsubstep(jg)
      CALL message(TRIM(routine),message_text)
    ELSE
      WRITE(message_text,'(a,i4,a)') 'Domain', jg, &
        '; No QV substepping'
      CALL message(TRIM(routine),message_text)
    ENDIF

    ! height indices for cloud classification
    !
    ih_clch(jg) = kstart_moist(jg)
    DO jk = 1, nlev
      jk1 = jk + nshift_total
      IF (0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < hbase_clch) THEN
        ih_clch(jg) = jk
        EXIT
      ENDIF
    ENDDO
    ih_clcm(jg) = ih_clch(jg) + 1
    DO jk = ih_clch(jg) + 1, nlev
      jk1 = jk + nshift_total
      IF (0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < hbase_clcm) THEN
        ih_clcm(jg) = jk
        EXIT
      ENDIF
    ENDDO
    WRITE(message_text,'(2(a,i4),i4)') 'Domain', jg, &
      '; high- and mid-level clouds in layers above ', ih_clch(jg), ih_clcm(jg)
    CALL message(TRIM(routine),message_text)


  END SUBROUTINE configure_nonhydrostatic


END MODULE mo_nonhydrostatic_config
