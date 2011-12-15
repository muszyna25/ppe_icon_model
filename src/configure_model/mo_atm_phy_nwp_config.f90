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
MODULE mo_atm_phy_nwp_config

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: max_dom, MAX_CHAR_LENGTH, itconv, itccov,  &
    &                               itrad, itradheat, itsso, itgscp, itsatad,  &
    &                               itupdate, itturb, itsfc, itgwd, itfastphy, &
    &                               iphysproc, iphysproc_short
  USE mo_exception,           ONLY: message, message_text !, finish

  USE mo_icoham_sfc_indices,  ONLY: init_sfc_indices

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: atm_phy_nwp_config, t_atm_phy_nwp_config, tcall_phy
  PUBLIC :: configure_atm_phy_nwp

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !!--------------------------------------------------------------------------
  !! Basic configuration setup for atm dynamics
  !!--------------------------------------------------------------------------
  TYPE :: t_atm_phy_nwp_config

    ! namelist variables

    INTEGER ::  inwp_gscp        !> microphysics
    INTEGER ::  inwp_satad       !! saturation adjustment
    INTEGER ::  inwp_convection  !! convection
    INTEGER ::  inwp_radiation   !! radiation
    INTEGER ::  inwp_sso         !! sso
    INTEGER ::  inwp_gwd         !! non-orographic gravity wave drag
    INTEGER ::  inwp_cldcover    !! cloud cover
    INTEGER ::  inwp_turb        !! turbulence
    INTEGER ::  inwp_surface     !! surface including soil, ocean, ice,lake

    REAL(wp) :: dt_conv    !> field element for convection
    REAL(wp) :: dt_ccov    !! field element for subscale cloud cover
    REAL(wp) :: dt_rad     !! "-"                     radiation
    REAL(wp) :: dt_sso     !! "-"  for subscale orographic gravity waves
    REAL(wp) :: dt_gwd     !! "-"  for subscale gravity waves
    REAL(wp) :: dt_fastphy !! field element for fast physics processes
                         !! microphysics, saturation adjustment, turbulence, 
                         !! surface (in addition: update and radheat)

    INTEGER :: imode_turb, itype_wcld, icldm_turb, itype_tran
    LOGICAL :: limpltkediff, ltkesso, lexpcor
    REAL(wp):: tur_len, pat_len, a_stab, tkhmin, tkmmin, c_diff, &
      &        rlam_heat, rlam_mom, rat_sea 

    REAL(wp) :: qi0, qc0
    REAL(wp) :: ustart_raylfric    !! velocity at which extra Rayleigh friction starts
    REAL(wp) :: efdt_min_raylfric  !! e-folding time corresponding to maximum relaxation 
                                   !! coefficient
    LOGICAL  :: latm_above_top     !! use extra layer above model top for radiation 
                                   !! (reduced grid only)


    ! Derived variables
    !
    LOGICAL :: lproc_on(iphysproc) !> contains information about status of 
                                   !! corresponding physical process
                                    !! ON: TRUE; OFF: FALSE

  END TYPE t_atm_phy_nwp_config

  !>
  !!
  TYPE(t_atm_phy_nwp_config) :: atm_phy_nwp_config(max_dom) !< shape: (n_dom)


  REAL(wp) ::  &                          !> Field of calling-time interval (seconds) for
    &  tcall_phy(max_dom,iphysproc_short) !! each domain and phys. process


CONTAINS

SUBROUTINE configure_atm_phy_nwp( n_dom, pat_level, ltestcase, dtime_adv )
 !-------------------------------------------------------------------------
  !
  !>
  !! Setup NWP physics
  !!
  !! Read namelist for physics. Choose the physical package and subsequent
  !! parameters.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-10-06)
  !! revision for restructuring by Kristina Froehlich MPI-M (2011-07-12)

  INTEGER, INTENT(IN) :: n_dom
  INTEGER, INTENT(IN) :: pat_level(n_dom)
  REAL(wp),INTENT(IN) :: dtime_adv
  LOGICAL, INTENT(IN) :: ltestcase

  INTEGER :: jg
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
    &      routine = 'mo_atm_phy_nwp_config:configure_atm_phy_nwp'


    CALL message(TRIM(routine), '')

    DO jg = 1,n_dom
 
     atm_phy_nwp_config(jg)%dt_fastphy = (dtime_adv/2._wp**(pat_level(jg) &
          &                            -  pat_level(1)))                  !seconds
    ENDDO

    tcall_phy(:,:) = 0._wp


    DO jg = 1,n_dom

      ! initialize lproc_on
      atm_phy_nwp_config(jg)%lproc_on(1:iphysproc)    = .FALSE.


      ! Fill derived variable lproc_on
      IF ( atm_phy_nwp_config(jg)%inwp_satad > 0 )               &
        &  atm_phy_nwp_config(jg)%lproc_on(itupdate)  = .TRUE.
 
      IF ( atm_phy_nwp_config(jg)%inwp_satad > 0 )               &
        &  atm_phy_nwp_config(jg)%lproc_on(itsatad)   = .TRUE. 

      IF ( atm_phy_nwp_config(jg)%inwp_convection > 0 )          &
        &  atm_phy_nwp_config(jg)%lproc_on(itconv)    = .TRUE. 

      IF ( atm_phy_nwp_config(jg)%inwp_cldcover > 0 )            &
        &  atm_phy_nwp_config(jg)%lproc_on(itccov)    = .TRUE.

      IF ( atm_phy_nwp_config(jg)%inwp_radiation > 0 )           &
        &  atm_phy_nwp_config(jg)%lproc_on(itrad)     = .TRUE.

      IF ( atm_phy_nwp_config(jg)%inwp_sso > 0 )                 &
        &  atm_phy_nwp_config(jg)%lproc_on(itsso)     = .TRUE.

      IF ( atm_phy_nwp_config(jg)%inwp_gscp > 0 )                &
        &  atm_phy_nwp_config(jg)%lproc_on(itgscp)    = .TRUE.

      IF ( atm_phy_nwp_config(jg)%inwp_turb > 0 )                &
        &  atm_phy_nwp_config(jg)%lproc_on(itturb)    = .TRUE.

      IF ( atm_phy_nwp_config(jg)%inwp_radiation > 0 )           &
        &  atm_phy_nwp_config(jg)%lproc_on(itradheat) = .TRUE.

      IF ( atm_phy_nwp_config(jg)%inwp_surface > 0 )             &
        &  atm_phy_nwp_config(jg)%lproc_on(itsfc)     = .TRUE.

      IF ( atm_phy_nwp_config(jg)%inwp_gwd > 0 )                 &
        &  atm_phy_nwp_config(jg)%lproc_on(itgwd)     = .TRUE.


      ! Slow physics:
      ! currently for each domain the same time intervals are set
      !
      ! Fast physics:
      ! for each fast physics process the time interval is set equal to the 
      ! time interval for advection.
      !

      ! Slow physics
      !
      tcall_phy(jg,itconv) = atm_phy_nwp_config(jg)% dt_conv    ! sec

      tcall_phy(jg,itsso)  = atm_phy_nwp_config(jg)% dt_sso     ! sec

      tcall_phy(jg,itgwd)  = atm_phy_nwp_config(jg)% dt_gwd     ! sec

      tcall_phy(jg,itrad)  = atm_phy_nwp_config(jg)% dt_rad     ! sec

      !> KF always call clouds after convection
      !! to ensure the proper output of Qx_tot
      IF ( atm_phy_nwp_config(jg)%lproc_on(itccov)   .AND.  &
        &  atm_phy_nwp_config(jg)%lproc_on(itconv) ) THEN
        tcall_phy(jg,itccov) = atm_phy_nwp_config(jg)% dt_conv    ! sec
      ELSE
        tcall_phy(jg,itccov) = atm_phy_nwp_config(jg)% dt_fastphy ! sec
      ENDIF


      ! Fast physics
      !
      tcall_phy(jg,itfastphy)   = atm_phy_nwp_config(jg)%dt_fastphy ! sec

    ENDDO  ! jg loop



    IF( atm_phy_nwp_config(1)%inwp_turb == 2) THEN
       CALL init_sfc_indices( ltestcase, 'APE' ) !call of a hydrostatic testcase
                                             ! to obtain the demanded parameters
    ENDIF



    ! issue a warning, if advective and convective timesteps are not synchronized
    !
    ! DR: a clean implementation would require to put the following lines into 
    ! a jg-loop.
    IF( MOD( dtime_adv,atm_phy_nwp_config(1)%dt_conv) /= 0._wp )  THEN
      WRITE(message_text,'(a,2F9.1)') &
        &'WARNING: convective and advective timesteps not synchronized: ', &
        & tcall_phy(1,itconv), dtime_adv
      CALL message(TRIM(routine), TRIM(message_text))
      WRITE(message_text,'(a,2F9.1)') &
        &'implicit synchronization in time_ctrl_physics: dt_conv !=!', &
        & REAL((FLOOR(atm_phy_nwp_config(1)%dt_conv/dtime_adv) + 1),wp) &
        & * dtime_adv
      CALL message(TRIM(routine), TRIM(message_text))
    ENDIF
 
END SUBROUTINE configure_atm_phy_nwp


END MODULE mo_atm_phy_nwp_config
