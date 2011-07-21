!>
!! Namelist for scalar transport
!!
!! These subroutines are called by control_model and do the transport 
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
MODULE mo_advection_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_nml,          ONLY: lrestart
  USE mo_run_config,          ONLY: ntracer
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_ntracer, max_dom,      &
    &                               imiura, imiura3, iup3, ippm_vcfl, ippm_v,   &
    &                               inol, islopel_sm, islopel_m, ifluxl_m,      &
    &                               ifluxl_sm, inol_v, islopel_vsm, islopel_vm, &
    &                               ifluxl_vpd
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: p_pe, p_io
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_advection_config,    ONLY: advection_config 

  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_transport_namelist


  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !----------------------------------!
  ! transport_nml namelist variables !
  !----------------------------------!

  CHARACTER(len=MAX_CHAR_LENGTH) :: &!< list of tracers to initialize
    &  ctracer_list


  INTEGER :: &                     !< selects horizontal transport scheme
    &  ihadv_tracer(max_ntracer)   !< 0: no horizontal advection
                                   !< 1: 1st order upwind
                                   !< 2: 2nd order muscl
                                   !< 3: 2nd order miura
                                   !< 4: 3rd order miura with quadr. reconstr.


  INTEGER :: &                     !< selects vertical transport scheme
    &  ivadv_tracer(max_ntracer)   !< 0 : no vertical advection
                                   !< 1 : 1st order upwind
                                   !< 2 : 2nd order muscl
                                   !< 20: 2nd order muscl for CFL>1
                                   !< 3 : 3rd order PPM
                                   !< 30: 3rd order PPM for CFL>1


  LOGICAL :: lvadv_tracer          !< if .TRUE., calculate vertical tracer advection
  LOGICAL :: lclip_tracer          !< if .TRUE., clip negative tracer values
  LOGICAL :: lstrang               !< if .TRUE., use complete Strang splitting
                                   !< (\Delta t/2 vert)+(\Delta t hor)+(\Delta t/2 vert)

  LOGICAL :: llsq_svd              !< least squares reconstruction with 
                                   !< singular value decomposition (TRUE) or 
                                   !< QR decomposition (FALSE) of design matrix A

  INTEGER :: &                     !< parameter used to select the limiter
    &  itype_vlimit(max_ntracer)   !< for vertical transport
                               

  INTEGER, TARGET :: &             !< parameter used to select the limiter
    &  itype_hlimit(max_ntracer)   !< for horizontal transport
                                   !< 0: no limiter
                                   !< 1: semi-monotonous slope limiter
                                   !< 2: monotonous slope limiter
                                   !< 3: monotonous flux limiter

  INTEGER :: iord_backtraj         !< parameter to select the spacial order
                                   !< of accuracy for the backward trajectory

  INTEGER :: igrad_c_miura         !< parameter used to select the gradient
                                   !< reconstruction method at cell center
                                   !< for second order miura scheme

  INTEGER :: ivcfl_max             !< determines stability range of vertical 
                                   !< ppm-scheme (approximate allowable maximum 
                                   !< CFL-number)

  REAL(wp):: upstr_beta_adv        !< later, it should be combined with 
                                   !< upstr_beta in non-hydrostatic namelist


  NAMELIST/transport_nml/ ihadv_tracer, ivadv_tracer, lvadv_tracer,       &
    &                     itype_vlimit, ivcfl_max, itype_hlimit,          &
    &                     iord_backtraj, lclip_tracer, ctracer_list,      &
    &                     igrad_c_miura, lstrang, upstr_beta_adv, llsq_svd


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
  !!  by Daniel Reinert, DWD (2011-05-07)
  !!
  SUBROUTINE read_transport_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index
    INTEGER :: jt          !< tracer loop index
    INTEGER :: z_nogo(2)   !< for consistency check

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_advection_nml: read_transport_nml'

    !-----------------------
    ! 1. default settings   
    !-----------------------
    ctracer_list    = ''
    ihadv_tracer(:) = imiura    ! miura horizontal advection scheme
    itype_hlimit(:) = ifluxl_m  ! monotonous flux limiter
    ivadv_tracer(:) = ippm_vcfl ! PPM vertical advection scheme
    itype_vlimit(:) = islopel_vsm ! semi-monotonous slope limiter
    ivcfl_max       = 5         ! CFL-stability range for vertical advection
    iord_backtraj   = 1         ! 1st order backward trajectory
    lvadv_tracer    = .TRUE.    ! vertical advection yes/no
    lclip_tracer    = .FALSE.   ! clipping of negative values yes/no
    lstrang         = .FALSE.   ! Strang splitting yes/no

    igrad_c_miura   = 1         ! MIURA linear least squares reconstruction

    upstr_beta_adv  = 1.0_wp    ! =1.0 selects 3rd order advection in up3
                                ! =0.0 selects 4th order advection in up3
    llsq_svd        = .FALSE.   ! apply QR-decomposition (FALSE)

!DR    iadv_slev(:,:)      = 1         ! vertical start level !DR should go into 
!DR                                    ! the derived variables-section since it 
!DR                                    ! is not part of our namelist

!DR special settings for global_cell_type=6 will be done during the 
!DR crosscheck. At this point it is important to get rid of any 
!DR dependencies. 
!DR    SELECT CASE (global_cell_type)
!DR    CASE (3)
!DR      ihadv_tracer(:) = imiura    ! miura horizontal advection scheme
!DR      itype_hlimit(:) = ifluxl_m  ! monotonous flux limiter
!DR    CASE (6)
!DR      ihadv_tracer(:) = iup3      ! 3rd order upwind horizontal advection scheme
!DR      itype_hlimit(:) = ifluxl_sm ! semi monotonous flux limiter
!DR    END SELECT


    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('transport_nml')
      READ(funit,NML=transport_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('transport_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, transport_nml)
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------

    ! flux computation methods - sanity check
    !
    IF ( ANY(ihadv_tracer(1:ntracer) > iup3) .OR.                     &
      &  ANY(ihadv_tracer(1:ntracer) < 0) )    THEN
      CALL finish( TRIM(routine),                                     &
        &  'incorrect settings for ihadv_tracer. Must be 0,1,2,3, or 4 ')
    ENDIF
    IF ( ANY(ivadv_tracer(1:ntracer) > ippm_v) .OR.                   &
      &  ANY(ivadv_tracer(1:ntracer) < 0)) THEN
      CALL finish( TRIM(routine),                                     &
        &  'incorrect settings for ivadv_tracer. Must be 0,1,2,3,20, or 30 ')
    ENDIF
    z_nogo(1) = islopel_sm
    z_nogo(2) = islopel_m
    DO jt=1,max_ntracer
      IF ( ihadv_tracer(jt)==imiura3 .AND. ANY(z_nogo==itype_hlimit(jt)) ) THEN
        CALL finish( TRIM(routine),                                   &
          &  'incorrect settings for MIURA3. No slope limiter available ')
      ENDIF
    END DO
    IF (upstr_beta_adv > 1.0_wp .OR. upstr_beta_adv < 0.0_wp) THEN
      CALL finish( TRIM(routine),                                     &
        &  'incorrect settings for upstr_beta_adv. Must be in [0,1] ')
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
      advection_config(jg)%iord_backtraj  = iord_backtraj
      advection_config(jg)%igrad_c_miura  = igrad_c_miura
      advection_config(jg)%ivcfl_max      = ivcfl_max
      advection_config(jg)%upstr_beta_adv = upstr_beta_adv
    ENDDO



    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=transport_nml)                    
    CALL store_and_close_namelist(funit, 'transport_nml')             


    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(p_pe == p_io) WRITE(nnml_output,nml=transport_nml)


  END SUBROUTINE read_transport_namelist

END MODULE mo_advection_nml
