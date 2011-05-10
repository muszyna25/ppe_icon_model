!>
!!  Namelist for NWP physics
!!
!!  these Subroutines are called by control model and construct the
!!  physics composition
!!
!! @author <Kristina Froehlich, DWD>
!!
!!
!! @par Revision History
!! First implementation by Kristina Froehlich, DWD (2010-06-20>)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_atm_phy_nwp_nml

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: max_dom,MAX_CHAR_LENGTH,itconv,itccov,&
    &                               itrad,itradheat, itsso,itgscp,itsatad,itupdate,&
    &                               itturb, itsfc,  iphysproc
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_namelist,            ONLY: position_nml, POSITIONED
  USE mo_mpi,                 ONLY: p_pe, p_io
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_nonhydrostatic_nml,  ONLY: iadv_rcf
  USE mo_dynamics_nml,        ONLY: ldry_dycore
  USE mo_io_nml,              ONLY: lwrite_extra
  USE mo_run_nml,             ONLY: dtime, ltestcase
  USE mo_model_domain,        ONLY: t_patch
  USE mo_model_domain_import, ONLY: n_dom
  USE mo_radiation_nml,       ONLY: read_radiation_nml, irad_o3
  USE mo_data_turbdiff,       ONLY: imode_turb,                              &
    &                               limpltkediff, ltkesso, lexpcor,          &
    &                               tur_len, pat_len, a_stab,                &
    &                               tkhmin, tkmmin, c_diff,                  &
    &                               itype_wcld, icldm_turb,                  &
    &                               itype_tran, rlam_heat, rlam_mom, rat_sea
  USE mo_echam_vdiff_nml,     ONLY: echam_vdiff_nml_setup
  USE mo_icoham_sfc_indices,  ONLY: init_sfc_indices, nsfc_type

  IMPLICIT NONE

  PRIVATE

!> Action Variables for physical schemes
! --------------------------------------
  INTEGER, PARAMETER ::   &
&      iupdate       = 1 ,& !! update moist tracers
&      isatad        = 2 ,& !! saturation adjustment
&      islow_physics = 3 ,& !! 'slow' physical processes
       ifast_physics = 4    !> cloud microphysics

   !
   ! user defined calling intervals
   !
  REAL(wp) :: dt_conv(max_dom)   !> field element for convection
  REAL(wp) :: dt_ccov(max_dom)   !! field element for subscale cloud cover
  REAL(wp) :: dt_rad(max_dom)    !! "-"                     radiation
  REAL(wp) :: dt_radheat(max_dom)!! "-" rad. heating from radiative fluxes with updated cosmu0 
  REAL(wp) :: dt_sso(max_dom)    !! "-"  for subscale orographic gravity waves
  REAL(wp) :: dt_gscp(max_dom)   !! field element for microphysics
  REAL(wp) :: dt_turb(max_dom)   !! field element for turbulence
  REAL(wp) :: dt_sfc(max_dom)    !! field element for surface
  REAL(wp) :: dt_satad(max_dom)  !! field element for sat. adjustment
  REAL(wp) :: dt_update(max_dom) !! field element for tracer phys update


  REAL(wp) ::  &                    !> Field of calling-time interval (seconds) for
    &  tcall_phy(max_dom,iphysproc) !! each domain and phys. process

  INTEGER ::  inwp_gscp        !> microphysics
  INTEGER ::  inwp_satad       !! saturation adjustment
  INTEGER ::  inwp_convection  !! convection
  INTEGER ::  inwp_radiation   !! radiation
  INTEGER ::  inwp_sso         !! sso
  INTEGER ::  inwp_cldcover    !! cloud cover
  INTEGER ::  inwp_turb        !! turbulence
  INTEGER ::  inwp_surface     !! surface including soil, ocean, ice,lake
  INTEGER ::  nlevs, nztlev    !! number of soil layers, time integration scheme
  INTEGER ::  nlev_snow        !! number of snow layers
  INTEGER ::  nsfc_subs        !! number of TILES

  INTEGER :: jg
                               !KF should be moved to run_ctl
!  INTEGER :: inextra_2d        !> number of extra output fields for debugging
!  INTEGER :: inextra_3d        !> number of extra output fields for debugging

  LOGICAL ::       &
       lseaice,    & !> forecast with sea ice model
       llake,      & !! forecst with lake model FLake
       l3dturb,    & !! 3D-turbulence (additional horizontal diffusion)
       lprog_tke,  & !! prognostic treatment of TKE (for itype_turb= 5/7)
       lmelt     , & !! soil model with melting process
       lmelt_var , & !! freezing temperature dependent on water content
       lmulti_snow   !! run the multi-layer snow model

!> Variables for hydci_pp
! --------------------------------------

  REAL (KIND=wp), PARAMETER ::              &
    qi0 = 0.0_wp ,  & !> cloud ice threshold for autoconversion
&   qc0 = 0.0_wp      !! cloud water threshold for autoconversion



!--------------------------------------------------------------------
! nwp forcing (right hand side)
!--------------------------------------------------------------------

  NAMELIST/nwp_phy_ctl/inwp_gscp, inwp_satad, inwp_convection, &
    &                  inwp_radiation, inwp_sso, inwp_cldcover, &
    &                  inwp_turb, inwp_surface, nlevs, nztlev,  &
    &                  nlev_snow, nsfc_subs,                    &
    &                  dt_conv, dt_ccov, dt_rad,                &
    &                  dt_radheat,                              &
    &                  dt_sso, dt_gscp, dt_satad,               &
    &                  dt_turb, dt_sfc,                         & !, inextra_2d,inextra_3d
    &                  imode_turb,                              &
    &                  limpltkediff, ltkesso, lexpcor,          &
    &                  tur_len, pat_len, a_stab,                &
    &                  tkhmin, tkmmin, c_diff,                  &
    &                  itype_wcld, icldm_turb,                  &
    &                  itype_tran, rlam_heat, rlam_mom, rat_sea

   PUBLIC :: ifast_physics, isatad, islow_physics, iupdate
   PUBLIC :: setup_nwp_phy !, set_inwp_nml, read_inwp_nml
   PUBLIC :: inwp_gscp, inwp_satad, inwp_convection, inwp_radiation
   PUBLIC :: inwp_sso, inwp_cldcover, inwp_turb, inwp_surface
   PUBLIC :: nlevs, nztlev
   PUBLIC :: nlev_snow,nsfc_subs
   PUBLIC :: dt_conv, dt_ccov, dt_rad, dt_radheat, dt_sso, dt_gscp
   PUBLIC :: dt_satad, dt_update,   dt_turb, dt_sfc,tcall_phy
   PUBLIC :: qi0, qc0
   PUBLIC :: lseaice, llake, l3dturb,  lprog_tke, lmulti_snow


 CONTAINS

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
  !!
  SUBROUTINE setup_nwp_phy( p_patch )

    TYPE(t_patch), OPTIONAL, INTENT(IN) :: p_patch(:)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = 'mo_atm_phy_nwp_nml'
    !-------------------------------------------------------------------------

    ldry_dycore     = .FALSE.

    !> set default physics switches and values
    IF (PRESENT(p_patch)) THEN
      CALL set_inwp_nml( p_patch )
    ELSE
      CALL set_inwp_nml
    ENDIF

    !
    !> final settings via namelist

    CALL read_inwp_nml


    IF ( inwp_radiation > 0 )  THEN
      
      CALL read_radiation_nml
      
      SELECT CASE (irad_o3)
      CASE (0,6)
        ! ok
      CASE default
        CALL finish('setup_nwp_phy: radiation_nml','irad_o3 currently has to be 0 or 6.')
      END SELECT

    ENDIF

    IF(inwp_turb == 2) THEN
       CALL echam_vdiff_nml_setup
       CALL init_sfc_indices( ltestcase, 'APE' ) !call of a hydrostatic testcase
                                             ! to obtain the demanded parameters
    ENDIF

      CALL message(TRIM(routine), 'nwp_physics namelist read in')
  
  END SUBROUTINE setup_nwp_phy



  !-------------------------------------------------------------------------
  !
  !>
  !! Set default values for physics Namelist
  !!
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!
  SUBROUTINE set_inwp_nml( p_patch )

    TYPE(t_patch), OPTIONAL, INTENT(IN) :: p_patch(:)

  !-----------------------------------------------------------------------

    inwp_gscp       = 0           !> 0 = no microphysics
    inwp_satad      = 0           !> 1 = saturation adjustment on
    inwp_convection = 0           !> 0 = no convection
    inwp_radiation  = 0           !> 0 = no radiation
    inwp_sso        = 0           !> 0 = no sso
    inwp_cldcover   = 1           !> 1 = use grid-scale clouds for radiation
    inwp_turb       = 0           !> 0 = no turbulence,1= cosmo/turbdiff,2=echam/vdiff
    inwp_surface    = 0           !> 0 = no surface, 1 =  cosmo surface

    nlevs           = 7           !> 7 = default value for number of soil layers
    nztlev          = 2           !> 2 = default value for time integration scheme
    nlev_snow       = 1           !> 1 = default value for number of snow layers
    nsfc_subs       = 1           !> 1 = default value for number of TILES

    ! initialize the following values with zero to allow for namelist output
    dt_conv (:)  = 0.0_wp
    dt_ccov (:)  = 0.0_wp
    dt_rad  (:)  = 0.0_wp
    dt_sso  (:)  = 0.0_wp
    dt_gscp (:)  = 0.0_wp
    dt_turb (:)  = 0.0_wp
    dt_sfc  (:)  = 0.0_wp
    dt_satad(:)  = 0.0_wp
    dt_update(:) = 0.0_wp
    dt_radheat(:)= 0.0_wp

    ! default values of time interval for each physical paramterization
    IF (PRESENT(p_patch)) THEN
      DO jg=1,n_dom
        dt_conv (jg) = 600._wp      !seconds
        dt_ccov (jg) = dt_conv(jg)  !presently not used; cloud cover is synchronized with radiation
        dt_rad  (jg) = 1800._wp     !seconds
        dt_sso  (jg) = dt_conv(jg)  !seconds

        !> intervals coupled to advective timestep
        !! all three processes have to take place at the same
        !! time step and frequency

        dt_gscp (jg) = ( REAL(iadv_rcf,wp)              &
          &             * (dtime/2._wp**(p_patch(jg)%level &
          &             - p_patch(1)%level)) )          !seconds
        dt_turb (jg) = ( REAL(iadv_rcf,wp)              &
          &             * (dtime/2._wp**(p_patch(jg)%level &
          &             - p_patch(1)%level)) )          !seconds
        dt_sfc  (jg) = ( REAL(iadv_rcf,wp)              &
          &             * (dtime/2._wp**(p_patch(jg)%level &
          &             - p_patch(1)%level)) )          !seconds
        ! satad is coupled to advective timestep
        dt_satad(jg) = ( REAL(iadv_rcf,wp)              &
                        * dtime/2._wp**(p_patch(jg)%level &
           &              - p_patch(1)%level))          !seconds
        ! coupled to advective timestep
        dt_update(jg) =  dt_satad (jg)
        dt_radheat(jg)=  dt_update(jg)

      ENDDO
    ELSE ! patch information is not available for call from io_async
      DO jg=1,n_dom
        dt_conv (jg) = 600._wp      !seconds
        dt_ccov (jg) = dt_conv(jg)  !presently not used; cloud cover is synchronized with radiation
        dt_rad  (jg) = 1800._wp     !seconds
        dt_sso  (jg) = 3600._wp     !seconds
        dt_gscp (jg) = 100._wp      !seconds
        dt_turb (jg) = 100._wp      !seconds
        dt_sfc  (jg) = 100._wp      !seconds
        dt_satad(jg) = 100._wp      !seconds
        dt_update(jg) =  dt_satad (jg)
        dt_radheat(jg)=  dt_update(jg)
      ENDDO
    ENDIF

  !> KF  current settings to get NWP turbulence running
        lseaice    = .FALSE.
        llake      = .FALSE.
        l3dturb    = .FALSE.
        lprog_tke  = .FALSE.!> prognostic treatment of TKE (for itype_turb=5/7)
        lmulti_snow= .FALSE.

  END SUBROUTINE set_inwp_nml

  !-------------------------------------------------------------------------
  !
  !>
  !! Read physics Namelist
  !!
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!
 SUBROUTINE read_inwp_nml

  INTEGER :: jg
  INTEGER :: i_status

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine =  &
                              'mo_atm_phy_nwp_nml'
  !-----------------------------------------------------------------------

    CALL position_nml ('nwp_phy_ctl', status=i_status)
    !
    SELECT CASE (i_status)
    CASE (POSITIONED)
      READ (nnml, nwp_phy_ctl)
    END SELECT
  !  write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=nwp_phy_ctl)

    tcall_phy(:,:) = 0._wp

!    consistency check

    IF( (inwp_convection >0 ) .OR. (inwp_gscp > 0)) inwp_satad = 1


    DO jg = 1,n_dom
      ! Slow physics:
      ! currently for each domain the same time intervals are set
      !
      ! Fast physics:
      ! time interval are set equal to the time interval for advection
      !
      ! note that time intervalls are set to 0 if the physical process
      ! is not considered at all.
      IF ( inwp_convection == 0 ) THEN   ! 0 = no convection
        tcall_phy(jg,itconv) = 0._wp
      ELSE
        tcall_phy(jg,itconv) = dt_conv(jg)    ! seconds
      ENDIF

       !> KF always call clouds after convection
       !! to ensure the proper output of Qx_tot

      IF ( inwp_convection > 0 ) THEN   ! 0 = no convection
        tcall_phy(jg,itccov) = dt_conv(jg)    ! seconds
      ELSE 
        tcall_phy(jg,itccov) = dt_gscp(jg)    ! seconds
      ENDIF
      !really switch off the clouds
      IF ( inwp_cldcover == 0 ) THEN     ! 0 = no cloud cover
        tcall_phy(jg,itccov) = 0._wp
      ENDIF

!      ELSE
!        tcall_phy(jg,itccov) = dt_rad(jg)           ! cloud cover
!                                                    ! should be coupled
!      ENDIF                                         ! to convection
                                                    ! and radiation!

      IF ( inwp_radiation == 0 ) THEN    ! 0 = no radiation
        tcall_phy(jg,itrad)     =  0._wp
        tcall_phy(jg,itradheat) =  0._wp
      ELSE
        tcall_phy(jg,itrad)     =  dt_rad    (jg)    ! seconds
        tcall_phy(jg,itradheat) =  dt_radheat(jg)    ! seconds       
      ENDIF

      IF ( inwp_sso == 0 ) THEN          ! 0 = no sso
        tcall_phy(jg,itsso) =  0._wp
      ELSE
        tcall_phy(jg,itsso) =  dt_sso(jg)           ! seconds
      ENDIF

      IF ( inwp_gscp == 0 ) THEN         ! 0 = no microphysics
        tcall_phy(jg,itgscp) =  0._wp
      ELSE

        tcall_phy(jg,itgscp) =  dt_gscp(jg)     ! seconds
      ENDIF

      IF ( inwp_satad == 0 ) THEN         ! 0 = no satad
        tcall_phy(jg,itsatad)  =  0._wp
        tcall_phy(jg,itupdate) =  0._wp   ! no moist update needed if no satad
      ELSE
        tcall_phy(jg,itsatad)  =  dt_satad(jg)   !seconds
        tcall_phy(jg,itupdate) =  dt_update(jg)  !seconds
      ENDIF

      IF ( inwp_turb == 0 ) THEN         ! 0 = no turbulence
        tcall_phy(jg,itturb) =  0._wp 
      ELSE
        tcall_phy(jg,itturb) =  dt_turb(jg)   !seconds
      ENDIF

      IF ( inwp_surface == 0 ) THEN         ! 0 = no soil
        tcall_phy(jg,itsfc) =  0._wp 
      ELSE
        tcall_phy(jg,itsfc) =  dt_turb(jg)   !seconds
      ENDIF

    ENDDO


! consistency checks

     IF( MOD( REAL(iadv_rcf,wp)*dtime, dt_conv(1)) /= 0._wp )  THEN
       WRITE(message_text,'(a,I4,2F10.2)') &
          &'advective and convective timesteps are not- but will be synchronized ',&
      &     1, REAL(iadv_rcf,wp)*dtime,tcall_phy(1,itconv)
      CALL message(TRIM(routine), TRIM(message_text))
     ENDIF

   IF( (inwp_gscp==0) .AND. (inwp_convection==0) .AND. (inwp_radiation==0) &
     & .AND. (inwp_sso==0)  .AND. (inwp_surface == 0) .AND. (inwp_turb> 0) )   &
     CALL message(TRIM(routine),' WARNING! NWP forcing set but only turbulence selected!')


 END SUBROUTINE read_inwp_nml
!==============================================================================

END MODULE mo_atm_phy_nwp_nml

