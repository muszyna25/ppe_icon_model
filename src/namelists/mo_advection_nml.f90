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
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_io_units,            ONLY: nnml,nnml_output
  USE mo_master_nml,          ONLY: lrestart
  USE mo_run_config,          ONLY: ntracer, ntracer_static, num_lev, nlev, &
    &                               iforcing, io3, iqcond, lvert_nest, iqv
  USE mo_grid_config,        ONLY: n_dom, global_cell_type
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_ntracer, max_dom,  &
    &                               ino_hadv, iup, imiura, imiura3, iup3,   &
    &                               ino_vadv, iup_v, imuscl_vcfl, imuscl_v, &
    &                               ippm_vcfl, ippm_v, inol, islopel_sm,    &
    &                               islopel_m, ifluxl_m, ifluxl_sm, inol_v, &
    &                               islopel_vsm, islopel_vm, ifluxl_vpd,    &
    &                               ino_flx, izero_grad, iparent_flx,       &
    &                               inoforcing, iheldsuarez, iecham, inwp,  &
    &                               ildf_dry, ildf_echam
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: p_pe, p_io
  USE mo_radiation_config,    ONLY: irad_o3
  USE mo_nonhydrostatic_config,  ONLY: l_open_ubc, kstart_moist, kstart_qv
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist, &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_advection_config 

  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: transport_nml_setup, setup_transport, read_transport_namelist

  PUBLIC :: iadv_slev, iubc_adv, iup, imiura, imiura3,                     &
    &       inol, islopel_sm, islopel_m, ifluxl_m, ifluxl_sm, iup_v,       &
    &       imuscl_v, imuscl_vcfl, ippm_v, ippm_vcfl, inol_v, islopel_vsm, & 
    &       islopel_vm, ifluxl_vpd, t_compute, t_cleanup, lcompute,        &
    &       lcleanup, iup3, ino_flx, izero_grad, iparent_flx

  PUBLIC :: shape_func, zeta, eta, wgt_zeta, wgt_eta


  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !----------------------------------!
  ! transport_nml namelist variables !
  !----------------------------------!

  CHARACTER(len=MAX_CHAR_LENGTH) :: &!< list of tracers to initialize
    &  ctracer_list


  INTEGER :: &                       !< selects horizontal transport scheme
    &  ihadv_tracer(max_ntracer)     !< 0: no horizontal advection
                                     !< 1: 1st order upwind
                                     !< 2: 2nd order muscl
                                     !< 3: 2nd order miura
                                     !< 4: 3rd order miura with quadr. reconstr.


  INTEGER :: &                       !< selects vertical transport scheme
    &  ivadv_tracer(max_ntracer)     !< 0 : no vertical advection
                                     !< 1 : 1st order upwind
                                     !< 2 : 2nd order muscl
                                     !< 20: 2nd order muscl for CFL>1
                                     !< 3 : 3rd order PPM
                                     !< 30: 3rd order PPM for CFL>1


  LOGICAL :: lvadv_tracer         !< if .TRUE., calculate vertical tracer advection
  LOGICAL :: lclip_tracer         !< if .TRUE., clip negative tracer values
  LOGICAL :: lstrang              !< if .TRUE., use complete Strang splitting
                                  !< (\Delta t/2 vert)+(\Delta t hor)+(\Delta t/2 vert)

  LOGICAL :: llsq_svd             !< least squares reconstruction with 
                                  !< singular value decomposition (TRUE) or 
                                  !< QR decomposition (FALSE) of design matrix A

  INTEGER :: &                    !< parameter used to select the limiter
    &  itype_vlimit(max_ntracer)  !< for vertical transport
                               

  INTEGER, TARGET :: &            !< parameter used to select the limiter
    &  itype_hlimit(max_ntracer)  !< for horizontal transport
                                  !< 0: no limiter
                                  !< 1: semi-monotonous slope limiter
                                  !< 2: monotonous slope limiter
                                  !< 3: monotonous flux limiter

  INTEGER :: iord_backtraj        !< parameter to select the spacial order
                                  !< of accuracy for the backward trajectory

  INTEGER :: igrad_c_miura        !< parameter used to select the gradient
                                  !< reconstruction method at cell center
                                  !< for second order miura scheme

  INTEGER :: ivcfl_max            !< determines stability range of vertical 
                                  !< ppm-scheme (approximate allowable maximum 
                                  !< CFL-number)

  REAL(wp) :: upstr_beta_adv      !< later, it should be combined with 
                                  !< upstr_beta in non-hydrostatic namelist


  NAMELIST/transport_nml/ ihadv_tracer, ivadv_tracer, lvadv_tracer,       &
    &                     itype_vlimit, ivcfl_max, itype_hlimit,          &
    &                     iord_backtraj, lclip_tracer, ctracer_list,      &
    &                     igrad_c_miura, lstrang, upstr_beta_adv, llsq_svd


  !-----------------------------!
  ! dependent control variables !
  !-----------------------------!
  REAL(wp) :: cSTR             !< if complete Strang-splitting is used,
                               !< this constant adapts the time step

  INTEGER  :: &                !< selects upper boundary condition 
    &  iubc_adv(max_dom)       !< for tracer transport
                               !< 0: no flux
                               !< 1: zero gradient 
                               !< 2: interpolated flux from parent grid

  INTEGER ::  &                !< selects vertical start level for each patch
    &  iadv_slev(max_dom,max_ntracer) !< and each tracer.

  REAL(wp) :: coeff_grid       !< parameter which is used to make the vertical
                               !< advection scheme applicable to a height 
                               !< based coordinate system (coeff_grid=-1)

  REAL(wp) :: shape_func(4,4)  !< shape functions for mapping the FFSL departure
                               !< region onto the standard rectangle (miura3 only)

  REAL(wp) ::              &   !< Gauss quadrature points in \zeta-\eta space
    &  zeta(4), eta(4)         !< (miura3 only)

  REAL(wp) :: wgt_zeta(4), &   !< Gauss quadrature weights for zeta and eta
    &         wgt_eta(4)       !< points (miura3 only)


  ! allowing for the onetime computation of tracer independent parts
  TYPE t_compute
    LOGICAL :: muscl_v (max_ntracer)
    LOGICAL :: ppm_v   (max_ntracer)
    LOGICAL :: miura_h (max_ntracer)
    LOGICAL :: miura3_h(max_ntracer)
  END TYPE t_compute

  TYPE t_cleanup
    LOGICAL :: muscl_v (max_ntracer)
    LOGICAL :: ppm_v   (max_ntracer)
    LOGICAL :: miura_h (max_ntracer)
    LOGICAL :: miura3_h(max_ntracer)
  END TYPE t_cleanup

  TYPE(t_compute) :: lcompute
  TYPE(t_cleanup) :: lcleanup


  REAL(wp), POINTER ::  &
    &  ptr_delp_mc_now(:,:,:) => NULL() !< pointer to old layer thickness
                                        !< at cell center
  REAL(wp), POINTER ::  &
    &  ptr_delp_mc_new(:,:,:) => NULL() !< pointer to new layer thickness
                                        !< at cell center


CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !!   Set up the configuration for tracer transport.
  !!
  !!
  !! @par Revision History
  !!  by Hui Wan, MPI-M (2007-02-23).
  !!  Modification by A. Gassmann, (2008-09-30)
  !!  - moved this subroutine from hydro_control here
  !!  Modification by Daniel Reinert, DWD (2009-12-16)
  !!  - included option for full Strang-splitting (lstrang=.TRUE.)
  !!  Modifications by Daniel Reinert, DWD (2010-01-29)
  !!  - removed MPDATA specific options
  !!  Modifications by Daniel Reinert, DWD (2010-03-05)
  !!  - moved this subroutine from mo_advection_stepping here
  !!  Modification by Daniel Reinert, DWD (2010-03-24)
  !!  - included option for backward trajectory of second order
  !!    accuracy.
  !!  Modification by Daniel Reinert, DWD (2010-04-23)
  !!  - included variable which is set to 1 in the case of an
  !!    eta coordinate or -1 for a Gal-Chen coordinate in
  !!    the vertical direction.
  !!  Modification by Almut Gassmann, MPI-M (2010-11-18)
  !!  - adaptions for including the hexagonal C-grid option
  !!  Modification by Daniel Reinert, DWD (2011-04-20)
  !!  - moved this subroutine from mo_advection_utils here
  !!
  SUBROUTINE transport_nml_setup
    !
    INTEGER :: istat, i_listlen, funit
    INTEGER :: jt          !< tracer loop index
    INTEGER :: z_nogo(2)   !< for consistency check

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_nml: transport_nml_setup'

    !-----------------------------------------------------------------------

    !
    ! 1. default settings
    !
    ctracer_list = ''
    SELECT CASE (global_cell_type)
    CASE (3)
      ihadv_tracer(:) = imiura    ! miura horizontal advection scheme
      itype_hlimit(:) = ifluxl_m  ! monotonous flux limiter
    CASE (6)
      ihadv_tracer(:) = iup3      ! 3rd order upwind horizontal advection scheme
      itype_hlimit(:) = ifluxl_sm ! semi monotonous flux limiter
    END SELECT
    ivadv_tracer(:) = ippm_vcfl   ! PPM vertical advection scheme
    itype_vlimit(:) = islopel_vsm ! semi-monotonous slope limiter
    ivcfl_max   = 5               ! CFL-stability range for vertical advection
    iadv_slev(:,:)  = 1           ! vertical start level
    iord_backtraj = 1             ! 1st order backward trajectory
    lvadv_tracer= .TRUE.          ! vertical advection yes/no
    lclip_tracer= .FALSE.         ! clipping of negative values yes/no
    lstrang     = .FALSE.         ! Strang splitting yes/no

    igrad_c_miura = 1             ! MIURA linear least squares reconstruction

    upstr_beta_adv = 1.0_wp       ! =1.0 selects 3rd order advection in up3
                                  ! =0.0 selects 4th order advection in up3
    llsq_svd    = .FALSE.         ! apply QR-decomposition

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !----------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('transport_nml')
      READ(funit,NML=transport_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL position_nml ('transport_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, transport_nml)
    END SELECT

    !
    ! 3. check consistency of NAMELIST-parameters
    !
    i_listlen = LEN_TRIM(ctracer_list)

    SELECT CASE ( iforcing )
    CASE ( inwp )
      
!!$      SELECT CASE (inwp_radiation)
!!$      CASE (0)
!!$        IF ( ntracer /= iqcond ) THEN
!!$          ntracer = iqcond
!!$          WRITE(message_text,'(a,i3)') 'Attention: according to physics, ntracer is set to',iqcond
!!$          CALL message(TRIM(routine),message_text)
!!$        ENDIF
!!$        IF ( i_listlen /= iqcond ) THEN
!!$          DO jt=1,ntracer
!!$            WRITE(ctracer_list(jt:jt),'(i1.1)')jt
!!$          ENDDO
!!$          WRITE(message_text,'(a)') &
!!$            & 'Attention: according to physics, ctracer_list is set to ',&
!!$            & ctracer_list(1:ntracer)
!!$          CALL message(TRIM(routine),message_text)
!!$        ENDIF
!!$      CASE (1)
!!$        ntracer_static = 1
!!$        IF ( ntracer /= iqcond ) THEN
!!$          ntracer = iqcond
!!$          WRITE(message_text,'(a,i3)') &
!!$            &  'Attention: according to physics, ntracer is set to', iqcond   
!!$          CALL message(TRIM(routine),message_text)
!!$          WRITE(message_text,'(a)') &
!!$            &  'In addition, there is one static tracer for O3'
!!$          CALL message(TRIM(routine),message_text)
!!$        ENDIF
!!$        IF ( i_listlen /= ntracer ) THEN
!!$          DO jt=1,ntracer
!!$            WRITE(ctracer_list(jt:jt),'(i1.1)')jt
!!$          ENDDO
!!$          WRITE(message_text,'(a)') &
!!$            & 'Attention: according to physics, ctracer_list is set to ',&
!!$            &   ctracer_list(1:ntracer)
!!$          CALL message(TRIM(routine),message_text)
!!$        ENDIF
!!$      CASE (2)
!!$        SELECT CASE (irad_o3)
!!$        CASE (0)
!!$          IF ( ntracer /= iqcond  ) THEN
!!$            ntracer = iqcond
!!$            WRITE(message_text,'(a,i3)') &
!!$              &  'Attention: according to physics, ntracer is set to', iqcond      
!!$            CALL message(TRIM(routine),message_text)
!!$          ENDIF
!!$          IF ( i_listlen /= ntracer ) THEN
!!$            DO jt=1,ntracer
!!$              WRITE(ctracer_list(jt:jt),'(i1.1)')jt
!!$            ENDDO
!!$            WRITE(message_text,'(a)') &
!!$              & 'Attention: according to physics, ctracer_list is set to ', &
!!$              &  ctracer_list(1:ntracer)
!!$            CALL message(TRIM(routine),message_text)
!!$          ENDIF
!!$        CASE (6)
!!$          ntracer_static = 1
!!$          IF ( ntracer /= iqcond  ) THEN
!!$            ntracer = iqcond
!!$            WRITE(message_text,'(a,i3)') &
!!$              &  'Attention: according to physics, ntracer is set to', iqcond
!!$            CALL message(TRIM(routine),message_text)           
!!$            WRITE(message_text,'(a)') &
!!$              &  'In addition, there is one static tracer for O3'
!!$            CALL message(TRIM(routine),message_text)           
!!$          ENDIF
!!$          IF ( i_listlen /= ntracer ) THEN
!!$            DO jt=1,ntracer
!!$              WRITE(ctracer_list(jt:jt),'(i1.1)')jt
!!$            ENDDO
!!$            WRITE(message_text,'(a)') &
!!$              & 'Attention: according to physics with radiation and O3 ', &
!!$              &  'ctracer_list is set to ', &
!!$              &  ctracer_list(1:ntracer)
!!$            CALL message(TRIM(routine),message_text)
!!$          ENDIF
!!$        END SELECT
!!$      END SELECT
!!$
!!$
!!$      IF ( ( inwp_radiation > 0 ) .AND. (irad_o3==0 .OR. irad_o3==6) ) THEN
!!$        IF ( ihadv_tracer(io3) /= 0 ) THEN
!!$          ihadv_tracer(io3) = 0
!!$          WRITE(message_text,'(a,i1,a)') &
!!$            & 'Attention: Since irad_o3 is set to ',irad_o3,', ihadv_tracer(io3) is set to 0.'
!!$          CALL message(TRIM(routine),message_text)
!!$        ENDIF
!!$        IF ( ivadv_tracer(io3) /= 0 ) THEN
!!$          ivadv_tracer(io3) = 0
!!$          WRITE(message_text,'(a,i1,a)') &
!!$            & 'Attention: Since irad_o3 is set to ',irad_o3,', ivadv_tracer(io3) is set to 0.'
!!$          CALL message(TRIM(routine),message_text)
!!$        ENDIF
!!$      ENDIF


    CASE (inoforcing, iheldsuarez, iecham, ildf_dry, ildf_echam)
    
      IF ( i_listlen < ntracer .AND. i_listlen /= 0 ) THEN
        ntracer = i_listlen
        CALL message(TRIM(routine),'number of tracers is adjusted according to given list')
      END IF
    
    END SELECT



    ! flux compuation methods - consistency check
    !
    IF ( ANY(ihadv_tracer(1:ntracer) > 4) .OR.                    &
      &  ANY(ihadv_tracer(1:ntracer) < 0))            THEN
      CALL finish( TRIM(routine),                                       &
           'incorrect settings for ihadv_tracer. Must be 0,1,2,3, or 4 ')
    ENDIF
    SELECT CASE (global_cell_type)
    CASE (3)
      IF ( ANY(ihadv_tracer(1:ntracer) > 3))            THEN
        CALL finish( TRIM(routine),                                       &
             'incorrect settings for TRI-C grid ihadv_tracer. Must be 0,1,2, or 3 ')
      ENDIF
    CASE (6)
      IF (ANY(ihadv_tracer(1:ntracer) == 2) .OR.                    &
       &  ANY(ihadv_tracer(1:ntracer) == 3))            THEN 
        CALL finish( TRIM(routine),                                       &
             'incorrect settings for HEX-C grid ihadv_tracer. Must be 0,1, or 4 ')
      ENDIF
    END SELECT
    IF ( ANY(ivadv_tracer(1:ntracer) > ippm_v) .OR.                 &
      &  ANY(ivadv_tracer(1:ntracer) < 0)) THEN
      CALL finish( TRIM(routine),                                       &
           'incorrect settings for ivadv_tracer. Must be 0,1,2,3,20, or 30 ')
    ENDIF
    z_nogo(1) = islopel_sm
    z_nogo(2) = islopel_m
    DO jt=1,ntracer
      IF ( ihadv_tracer(jt) == imiura3 .AND. ANY( z_nogo == itype_hlimit(jt)) ) THEN
        CALL finish( TRIM(routine),                                     &
         'incorrect settings for MIURA3. No slope limiter available ')
      ENDIF
    END DO
    IF (upstr_beta_adv > 1.0_wp .OR. upstr_beta_adv < 0.0_wp) THEN
      CALL finish( TRIM(routine),                                       &
           'incorrect settings for upstr_beta_adv. Must be in [0,1] ')
    ENDIF


    ! limiter - consistency check
    !
    IF ( ANY(itype_vlimit(1:ntracer) < inol_v ) .OR.                &
      &  ANY(itype_vlimit(1:ntracer) > ifluxl_vpd)) THEN
      CALL finish( TRIM(routine),                                       &
       'incorrect settings for itype_vlimit. Must be 0,1,2 or 4 ')
    ENDIF
    IF ( ANY(itype_hlimit(1:ntracer) < inol ) .OR.                      &
      &  ANY(itype_hlimit(1:ntracer) > ifluxl_sm)) THEN
      CALL finish( TRIM(routine),                                       &
       'incorrect settings for itype_hlimit. Must be 0,1,2,3 or 4 ')
    ENDIF
    IF (global_cell_type == 6) THEN
      IF ( ANY(itype_hlimit(1:ntracer) == islopel_sm ) .OR.             &
        &  ANY(itype_hlimit(1:ntracer) == islopel_m  ) .OR.             &
        &  ANY(itype_hlimit(1:ntracer) == ifluxl_m   )) THEN
        CALL finish( TRIM(routine),                                     &
         'incorrect settings for itype_hlimit and hexagonal grid. Must be 0 or 4 ')
      ENDIF
    ENDIF

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=transport_nml)                                                             
    CALL store_and_close_namelist(funit, 'transport_nml')                                      

    ! 4. write the contents of the namelist to an ASCII file
    !
    IF(p_pe == p_io) WRITE(nnml_output,nml=transport_nml)


  END SUBROUTINE transport_nml_setup



  !>
  !! setup components of the transport scheme depending on this namelist
  !!
  !! Setup of additional transport control variables depending on the 
  !! transport-NAMELIST and potentially other namelists. This routine is 
  !! called, after all namelists have been read and a synoptic consistency 
  !! check has been done.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-04-20)
  !!
  SUBROUTINE setup_transport( iequations )
  !
    INTEGER,INTENT(IN) :: iequations
    INTEGER :: jg          !< loop index for gauss quadrature points
    INTEGER :: jt          !< tracer loop index
    !-----------------------------------------------------------------------


    !
    ! set dependent transport variables/model components, depending on 
    ! the transport namelist and potentially other namelsists.
    !

    ! check, whether Strang-splitting has been chosen and adapt cSTR accordingly
    IF ( lstrang ) THEN
      cSTR = 0.5_wp
    ELSE
      cSTR = 1._wp
    ENDIF


    ! Set grid-coefficient according to the applied vertical grid.
    !
    ! coeff_grid=1   : pressure based vertical coordinate system
    ! coeff_grid=-1  : height based vertical coordinate system
    !
    IF (iequations == 3) THEN  ! non-hydrostatic equation-set
      coeff_grid = -1._wp
    ELSE
      coeff_grid = 1._wp
    ENDIF


    !
    ! set vertical start level for each patch and each tracer
    !
    IF (iforcing == inwp) THEN
      DO jg = 1, n_dom
        ! Set iadv_slev to kstart_moist for all tracers but QV
        iadv_slev(jg,:)   = kstart_moist(jg)
        iadv_slev(jg,iqv) = kstart_qv(jg)
      ENDDO
    ELSE
      iadv_slev(1:n_dom,:) = 1
    ENDIF


    ! set boundary condition for vertical transport
    !
    IF (iequations == 3) THEN  ! non-hydrostatic equation-set
      DO jg = 1, n_dom

        IF (.NOT. lvert_nest ) THEN ! no vertical nesting

          IF (l_open_ubc) THEN
            iubc_adv(jg) = izero_grad ! zero gradient ubc
          ELSE
            iubc_adv(jg) = ino_flx    ! no flux ubc
          ENDIF

        ELSE ! vertical nesting

          IF (num_lev(jg) < nlev) THEN
            iubc_adv(jg) = iparent_flx
          ELSE IF ( (num_lev(jg) >= nlev) .AND. l_open_ubc) THEN
            iubc_adv(jg) = izero_grad
          ELSE IF ( (num_lev(jg) >= nlev) .AND. .NOT. l_open_ubc) THEN
            iubc_adv(jg) = ino_flx
          ENDIF
        ENDIF
      ENDDO

    ELSE ! hydrostatic or shallow water equation set
      iubc_adv(:) = ino_flx    ! no flux ubc
    ENDIF


    ! MUSCL_V[CFL] specific settings (vertical transport)
    !
    lcompute%muscl_v(:) = .FALSE.
    lcleanup%muscl_v(:) = .FALSE.

    IF ( ANY(ivadv_tracer == imuscl_v) .OR. ANY(ivadv_tracer == imuscl_vcfl)  ) THEN
      ! Search for the first tracer jt for which vertical advection of
      ! type MUSCL has been selected.
      DO jt=1,ntracer
        IF ( ivadv_tracer(jt) == imuscl_v .OR. ivadv_tracer(jt) == imuscl_vcfl ) THEN
          lcompute%muscl_v(jt) = .TRUE.
          exit
        ENDIF
      ENDDO

      ! Search for the last tracer jt for which vertical advection of
      ! type MUSCL has been selected.
      DO jt=ntracer,1,-1
        IF ( ivadv_tracer(jt) == imuscl_v .OR. ivadv_tracer(jt) == imuscl_vcfl ) THEN
          lcleanup%muscl_v(jt) = .TRUE.
          exit
        ENDIF
      ENDDO
    END IF


    ! PPM_V[CFL] specific settings (vertical transport)
    !
    lcompute%ppm_v(:)   = .FALSE.
    lcleanup%ppm_v(:)   = .FALSE.

    IF ( ANY(ivadv_tracer == ippm_v) .OR. ANY(ivadv_tracer == ippm_vcfl)  ) THEN
      ! Search for the first tracer jt for which vertical advection of
      ! type PPM has been selected.
      DO jt=1,ntracer
        IF ( ivadv_tracer(jt) == ippm_v .OR. ivadv_tracer(jt) == ippm_vcfl ) THEN
          lcompute%ppm_v(jt) = .TRUE.
          exit
        ENDIF
      ENDDO

      ! Search for the last tracer jt for which vertical advection of
      ! type PPM has been selected.
      DO jt=ntracer,1,-1
        IF ( ivadv_tracer(jt) == ippm_v .OR. ivadv_tracer(jt) == ippm_vcfl ) THEN
          lcleanup%ppm_v(jt) = .TRUE.
          exit
        ENDIF
      ENDDO
    END IF


    !
    ! MIURA specific settings
    !
    lcompute%miura_h(:) = .FALSE.
    lcleanup%miura_h(:) = .FALSE.

    IF ( ANY(ihadv_tracer(1:ntracer) == imiura) ) THEN
      ! Search for the first tracer jt for which horizontal advection of
      ! type MIURA has been selected.
      DO jt=1,ntracer
        IF ( ihadv_tracer(jt) == imiura ) THEN
          lcompute%miura_h(jt) = .TRUE.
          exit
        ENDIF
      ENDDO

      ! Search for the last tracer jt for which horizontal advection of
      ! type MIURA has been selected.
      DO jt=ntracer,1,-1
        IF ( ihadv_tracer(jt) == imiura ) THEN
          lcleanup%miura_h(jt) = .TRUE.
          exit
        ENDIF
      ENDDO
    END IF


    !
    ! MIURA3 specific settings
    !
    lcompute%miura3_h(:) = .FALSE.
    lcleanup%miura3_h(:) = .FALSE.

    ! Search for the first tracer jt for which horizontal advection of
    ! type MIURA3 has been selected.
    DO jt=1,ntracer
      IF ( ihadv_tracer(jt) == imiura3 ) THEN
        lcompute%miura3_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO

    ! Search for the last tracer jt for which horizontal advection of
    ! type MIURA3 has been selected.
    DO jt=ntracer,1,-1
      IF ( ihadv_tracer(jt) == imiura3 ) THEN
        lcleanup%miura3_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO


    !
    ! Compute shape functions for mapping the departure region onto the
    ! standard rectangle. It is assumed that a second order Gauss-Legendre
    ! quadrature is applied
    !

    ! Coordinates of integration points (in \zeta,\eta-System)
    zeta(1) = -1._wp/SQRT(3._wp)
    zeta(2) =  1._wp/SQRT(3._wp)
    zeta(3) =  1._wp/SQRT(3._wp)
    zeta(4) = -1._wp/SQRT(3._wp)

    eta(1)  = -1._wp/SQRT(3._wp)
    eta(2)  = -1._wp/SQRT(3._wp)
    eta(3)  =  1._wp/SQRT(3._wp)
    eta(4)  =  1._wp/SQRT(3._wp)


    ! shape functions for mapping
    DO jg = 1,4
      shape_func(1,jg) = 0.25_wp * (1._wp-zeta(jg))*(1._wp-eta(jg))
      shape_func(2,jg) = 0.25_wp * (1._wp+zeta(jg))*(1._wp-eta(jg))
      shape_func(3,jg) = 0.25_wp * (1._wp+zeta(jg))*(1._wp+eta(jg))
      shape_func(4,jg) = 0.25_wp * (1._wp-zeta(jg))*(1._wp+eta(jg))
    END DO

    ! Gauss quadrature weights
    wgt_zeta(1) = 1._wp
    wgt_zeta(2) = 1._wp
    wgt_zeta(3) = 1._wp
    wgt_zeta(4) = 1._wp

    wgt_eta(1)  = 1._wp
    wgt_eta(2)  = 1._wp
    wgt_eta(3)  = 1._wp
    wgt_eta(4)  = 1._wp


  END SUBROUTINE setup_transport


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
    llsq_svd        = .FALSE.   ! apply QR-decomposition

    iadv_slev(:,:)      = 1         ! vertical start level !DR should go into 
                                    ! the derived variables-section since it 
                                    ! is not part of our namelist

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


    ! flux compuation methods - sanity check
    !
    IF ( ANY(ihadv_tracer(1:ntracer) > 4) .OR.                    &
      &  ANY(ihadv_tracer(1:ntracer) < 0) )    THEN
      CALL finish( TRIM(routine),                                     &
        &  'incorrect settings for ihadv_tracer. Must be 0,1,2,3, or 4 ')
    ENDIF
    IF ( ANY(ivadv_tracer(1:ntracer) > ippm_v) .OR.               &
      &  ANY(ivadv_tracer(1:ntracer) < 0)) THEN
      CALL finish( TRIM(routine),                                     &
        &  'incorrect settings for ivadv_tracer. Must be 0,1,2,3,20, or 30 ')
    ENDIF
    z_nogo(1) = islopel_sm
    z_nogo(2) = islopel_m
    DO jt=1,max_ntracer
      IF ( ihadv_tracer(jt) == imiura3 .AND.                      &
        &  ANY( z_nogo == itype_hlimit(jt)) ) THEN
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
    IF ( ANY(itype_vlimit(1:ntracer) < inol_v ) .OR.              &
      &  ANY(itype_vlimit(1:ntracer) > ifluxl_vpd)) THEN
      CALL finish( TRIM(routine),                                     &
        &  'incorrect settings for itype_vlimit. Must be 0,1,2 or 4 ')
    ENDIF
    IF ( ANY(itype_hlimit(1:ntracer) < inol ) .OR.                &
      &  ANY(itype_hlimit(1:ntracer) > ifluxl_sm)) THEN
      CALL finish( TRIM(routine),                                     &
        &  'incorrect settings for itype_hlimit. Must be 0,1,2,3 or 4 ')
    ENDIF



    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 1,max_dom
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
