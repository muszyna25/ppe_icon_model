!>
!! Some utilities which are specific to the transport algorithm.
!!
!! Module contains some functions and procedures which are specifically related
!! to the transport schemes. These subroutines or functions are needed at
!! various places within the transport scheme. Therefore outsourcing these
!! routines protects from possible circular dependencies.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2010-03-04)
!! Modification by Daniel Reinert, DWD (2010-04-23)
!! - implemented generalized Lax-Friedrich flux function
!!   laxfr_upflux_v, which allows to use the same transport
!!   code for pressure and height based vertical coordinate
!!   systems.
!! Modification by Daniel Reinert, DWD (2010-05-17)
!! - added subroutines back_traj_dreg_o1, prep_gauss_quadrature and function
!!   jac which are part of the Gauss-Legendre quadrature apllied in the
!!   Miura-scheme.
!! Modification by Daniel Reinert, DWD (2010-10-14)
!! - added subroutine prep_gauss_quadrature_c for integrating a cubic polynomial.
!!   Renamed old prep_gauss_quadrature to prep_gauss_quadrature_q
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
MODULE mo_advection_utils

  USE mo_atm_phy_nwp_nml,     ONLY: inwp_radiation
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_model_domain,        ONLY: t_patch
  USE mo_interpolation,       ONLY: t_int_state
  USE mo_io_units,            ONLY: nnml,nnml_output
  USE mo_run_nml,             ONLY: ntracer, num_lev, nlev, nproma, iequations, & 
    &                               i_cell_type, iforcing, inwp, io3, iqt,      &
    &                               iqcond, ntracer_static, lvert_nest
  USE mo_grid_nml,            ONLY: n_dom
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_ntracer, max_dom,      &
    &                               min_rlcell_int, min_rledge_int
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_namelist,            ONLY: position_nml, POSITIONED
  USE mo_mpi,                 ONLY: p_pe, p_io
  USE mo_radiation_nml,       ONLY: irad_o3
  USE mo_nonhydrostatic_nml,  ONLY: l_open_ubc

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  REAL(wp) :: cSTR          !< if complete Strang-splitting is used,
                            !< this constant adapts the time step

  !
  ! transport_ctl
  !
  CHARACTER(len=MAX_CHAR_LENGTH) :: ctracer_list
                        !< list of tracers to initialize

  INTEGER :: ihadv_tracer(max_ntracer) !< selects horizontal transport scheme
                                       !< 0: no horizontal advection
                                       !< 1: 1st order upwind
                                       !< 2: 2nd order muscl
                                       !< 3: 2nd order miura
                                       !< 4: 3rd order miura with quadr. reconstr.

  INTEGER :: ivadv_tracer(max_ntracer) !< selects vertical transport scheme
                                       !< 0 : no vertical advection
                                       !< 1 : 1st order upwind
                                       !< 2 : 2nd order muscl
                                       !< 20: 2nd order muscl for CFL>1
                                       !< 3 : 3rd order PPM
                                       !< 30: 3rd order PPM for CFL>1

  LOGICAL :: lvadv_tracer  !< if .TRUE., calculate vertical tracer advection
  LOGICAL :: lclip_tracer  !< if .TRUE., clip negative tracer values
  LOGICAL :: lstrang       !< if .TRUE., use complete Strang splitting
                           !< (\Delta t/2 vert)+(\Delta t hor)+(\Delta t/2 vert)

  INTEGER :: itype_vlimit(max_ntracer) !< parameter used to select the limiter
                                       !< for vertical transport

  INTEGER, TARGET :: itype_hlimit(max_ntracer) 
                                       !< parameter used to select the limiter
                                       !< for horizontal transport
                                       !< 0: no limiter
                                       !< 1: semi-monotonous slope limiter
                                       !< 2: monotonous slope limiter
                                       !< 3: monotonous flux limiter

  INTEGER :: iadv_slev(max_ntracer)    !< selects vertical start level

  INTEGER :: iord_backtraj !< parameter to select the spacial order
                           !< of accuracy for the backward trajectory

  INTEGER :: igrad_c_miura   !< parameter used to select the gradient
                             !< reconstruction method at cell center
                             !< for second order miura scheme

  INTEGER :: iubc_adv(max_dom) !< selects upper boundary condition 
                               !< for tracer transport
                               !< 0: no flux
                               !< 1: zero gradient 
                               !< 2: interpolated flux from parent grid

  ! auxiliary variables for selecting horizontal transport scheme
  INTEGER, PARAMETER :: ino_hadv= 0
  INTEGER, PARAMETER :: iup     = 1
  INTEGER, PARAMETER :: imiura  = 2
  INTEGER, PARAMETER :: imiura3 = 3
  INTEGER, PARAMETER :: iup3    = 4

  ! auxiliary variables for selecting vertical transport scheme
  INTEGER, PARAMETER :: ino_vadv    = 0
  INTEGER, PARAMETER :: iup_v       = 1
  INTEGER, PARAMETER :: imuscl_vcfl = 2
  INTEGER, PARAMETER :: imuscl_v    = 20
  INTEGER, PARAMETER :: ippm_vcfl   = 3
  INTEGER, PARAMETER :: ippm_v      = 30

  ! auxiliary variables for selecting horizontal limiter
  INTEGER, PARAMETER :: inol       = 0
  INTEGER, PARAMETER :: islopel_sm = 1
  INTEGER, PARAMETER :: islopel_m  = 2
  INTEGER, PARAMETER :: ifluxl_m   = 3
  INTEGER, PARAMETER :: ifluxl_sm  = 4

  ! auxiliary variables for selecting vertical limiter
  INTEGER, PARAMETER :: inol_v      = 0
  INTEGER, PARAMETER :: islopel_vsm = 1
  INTEGER, PARAMETER :: islopel_vm  = 2
  INTEGER, PARAMETER :: ifluxl_vpd  = 4

  ! auxiliary variables for selecting upper boundary condition (ubc)
  INTEGER, PARAMETER :: ino_flx     = 0
  INTEGER, PARAMETER :: izero_grad  = 1
  INTEGER, PARAMETER :: iparent_flx = 2

  ! auxiliary variables allowing for the onetime computation of tracer
  ! independent parts
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

  REAL(wp) :: upstr_beta_adv ! later, it should be combined with upstr_beta in
                             ! non-hydrostatic namelist

  NAMELIST/transport_ctl/ ihadv_tracer, ivadv_tracer, lvadv_tracer,   &
    &                     itype_vlimit, itype_hlimit, iord_backtraj,  &
    &                     lclip_tracer, ctracer_list, igrad_c_miura,  &
    &                     lstrang, upstr_beta_adv


  ! In order to avoid circular dependencies these two pointers
  ! have been moved from mo_advection_stepping to this module.
  REAL(wp), POINTER ::  &
    &  ptr_delp_mc_now(:,:,:) => NULL() !< pointer to old layer thickness
                                        !< at cell center
  REAL(wp), POINTER ::  &
    &  ptr_delp_mc_new(:,:,:) => NULL() !< pointer to new layer thickness
                                        !< at cell center

  REAL(wp) :: coeff_grid  !< parameter which is used to make the vertical
                          !< advection scheme applicable to a height based
                          !< coordinate system (coeff_grid=-1)

  REAL(wp) :: shape_func(4,4) !< shape functions for mapping the FFSL departure
                              !< region onto the standard rectangle (miura3 only)

  REAL(wp) ::          &      !< Gauss quadrature points in \zeta-\eta space
    &  zeta(4), eta(4)        !< (miura3 only)

  REAL(wp) :: wgt_zeta(4), &  !< Gauss quadrature weights for zeta and eta
    &         wgt_eta(4)      !< points (miura3 only)

  PUBLIC :: tupdate_tracer, laxfr_upflux, laxfr_upflux_v, ptr_delp_mc_now,   &
    &       ptr_delp_mc_new, transport_ctl, ctracer_list, ihadv_tracer,      &
    &       ivadv_tracer, lvadv_tracer, lclip_tracer, lstrang, itype_vlimit, &
    &       itype_hlimit, iord_backtraj, igrad_c_miura, iadv_slev, iubc_adv, &
    &       cSTR, coeff_grid, setup_transport, back_traj_o1,                 &
    &       back_traj_dreg_o1, back_traj_o2, prep_gauss_quadrature_q,        &
    &       prep_gauss_quadrature_cpoor, prep_gauss_quadrature_c,            &
    &       upstr_beta_adv

  PUBLIC :: iup, imiura, imiura3, inol, islopel_sm, islopel_m, ifluxl_m,     &
    &       ifluxl_sm, iup_v, imuscl_v, imuscl_vcfl, ippm_v, ippm_vcfl,      &
    &       inol_v, islopel_vsm, islopel_vm, ifluxl_vpd, t_compute,          &
    &       t_cleanup, lcompute, lcleanup, iup3, ino_flx, izero_grad,        &
    &       iparent_flx

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
  !!
  SUBROUTINE setup_transport
    !
    INTEGER :: i_status, i_listlen
    INTEGER :: jg          !< loop index for gass quadrature points
    INTEGER :: jt          !< tracer loop index
    INTEGER :: z_nogo(2)   !< for consistency check

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_utils: setup_transport'

    !-----------------------------------------------------------------------

    ! default settings
    ctracer_list = ''
    SELECT CASE (i_cell_type)
    CASE (3)
      ihadv_tracer(:) = imiura    ! miura horizontal advection scheme
      itype_hlimit(:) = ifluxl_m  ! monotonous flux limiter
    CASE (6)
      ihadv_tracer(:) = iup3      ! 3rd order upwind horizontal advection scheme
      itype_hlimit(:) = ifluxl_sm ! semi monotonous flux limiter
    END SELECT
    ivadv_tracer(:) = ippm_vcfl   ! PPM vertical advection scheme
    itype_vlimit(:) = islopel_vsm ! semi-monotonous slope limiter
    iadv_slev(:)    = 1           ! vertical start level
    iord_backtraj   = 1           ! backward trajectory of 1st order accuracy
    lvadv_tracer    = .TRUE.      ! vertical advection yes/no
    lclip_tracer    = .FALSE.     ! clipping of negative values yes/no
    lstrang         = .FALSE.     ! Strang splitting yes/no

    igrad_c_miura = 1             ! MIURA specific settings (linear least squares)

    upstr_beta_adv = 1.0_wp       ! =1.0 selects 3rd order advection in up3
                                  ! =0.0 selects 4th order advection in up3

    ! read namelist
    CALL position_nml ('transport_ctl', status=i_status)
    SELECT CASE (i_status)
    CASE (POSITIONED)
      READ (nnml, transport_ctl)
    END SELECT

    i_listlen = LEN_TRIM(ctracer_list)

    IF (iforcing == inwp) THEN
      
      SELECT CASE (inwp_radiation)
      CASE (0)
        IF ( ntracer /= iqcond ) THEN
          ntracer = iqcond
          WRITE(message_text,'(a,i3)') 'Attention: according to physics, ntracer is set to',iqcond
          CALL message(TRIM(routine),message_text)
        ENDIF
        IF ( i_listlen /= iqcond ) THEN
          DO jt=1,ntracer
            WRITE(ctracer_list(jt:jt),'(i1.1)')jt
          ENDDO
          WRITE(message_text,'(a)') &
            & 'Attention: according to physics, ctracer_list is set to '//ctracer_list(1:ntracer)
          CALL message(TRIM(routine),message_text)
        ENDIF
      CASE (1)
        ntracer_static = 1
        IF ( ntracer /= iqcond ) THEN
          ntracer = iqcond
          WRITE(message_text,'(a,i3)') &
            &  'Attention: according to physics, ntracer is set to', iqcond   
          CALL message(TRIM(routine),message_text)
          WRITE(message_text,'(a)') &
            &  'In addition, there is one static tracer for O3'
          CALL message(TRIM(routine),message_text)
        ENDIF
        IF ( i_listlen /= ntracer ) THEN
          DO jt=1,ntracer
            WRITE(ctracer_list(jt:jt),'(i1.1)')jt
          ENDDO
          WRITE(message_text,'(a)') &
            & 'Attention: according to physics, ctracer_list is set to '// &
            &   ctracer_list(1:ntracer)
          CALL message(TRIM(routine),message_text)
        ENDIF
      CASE (2)
        SELECT CASE (irad_o3)
        CASE (0)
          IF ( ntracer /= iqcond  ) THEN
            ntracer = iqcond
            WRITE(message_text,'(a,i3)') &
              &  'Attention: according to physics, ntracer is set to', iqcond      
            CALL message(TRIM(routine),message_text)
          ENDIF
          IF ( i_listlen /= ntracer ) THEN
            DO jt=1,ntracer
              WRITE(ctracer_list(jt:jt),'(i1.1)')jt
            ENDDO
            WRITE(message_text,'(a)') &
              & 'Attention: according to physics, ctracer_list is set to '// &
              &  ctracer_list(1:ntracer)
            CALL message(TRIM(routine),message_text)
          ENDIF
        CASE (6)
          ntracer_static = 1
          IF ( ntracer /= iqcond  ) THEN
            ntracer = iqcond
            WRITE(message_text,'(a,i3)') &
              &  'Attention: according to physics, ntracer is set to', iqcond
            CALL message(TRIM(routine),message_text)           
            WRITE(message_text,'(a)') &
              &  'In addition, there is one static tracer for O3'
            CALL message(TRIM(routine),message_text)           
          ENDIF
          IF ( i_listlen /= ntracer ) THEN
            DO jt=1,ntracer
              WRITE(ctracer_list(jt:jt),'(i1.1)')jt
            ENDDO
            WRITE(message_text,'(a)') &
              & 'Attention: according to physics with radiation and O3, ctracer_list is set to ' &
              &  //ctracer_list(1:ntracer)
            CALL message(TRIM(routine),message_text)
          ENDIF
        END SELECT
      END SELECT
      
    ELSE ! iforcing /= inwp:
    
      IF ( i_listlen < ntracer .AND. i_listlen /= 0 ) THEN
        ntracer = i_listlen
        CALL message(TRIM(routine),'number of tracers is adjusted according to given list')
      END IF
    
    ENDIF !iforcing

    SELECT CASE ( iforcing )
    CASE ( inwp )
      IF ( ( inwp_radiation > 0 ) .AND. (irad_o3==0 .OR. irad_o3==6) ) THEN
        IF ( ihadv_tracer(io3) /= 0 ) THEN
          ihadv_tracer(io3) = 0
          WRITE(message_text,'(a,i1,a)') &
            & 'Attention: Since irad_o3 is set to ',irad_o3,', ihadv_tracer(io3) is set to 0.'
          CALL message(TRIM(routine),message_text)
        ENDIF
        IF ( ivadv_tracer(io3) /= 0 ) THEN
          ivadv_tracer(io3) = 0
          WRITE(message_text,'(a,i1,a)') &
            & 'Attention: Since irad_o3 is set to ',irad_o3,', ivadv_tracer(io3) is set to 0.'
          CALL message(TRIM(routine),message_text)
        ENDIF
      ENDIF
      ! CASE ( iecham )
      !
    END SELECT
    
    ! write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=transport_ctl)


    ! check, whether Strang-splitting is used and adapt cSTR accordingly
    IF (lstrang) THEN
      cSTR = 0.5_wp
    ELSE
      cSTR = 1._wp
    ENDIF


    !
    ! Set coefficient according to the applied vertical grid.
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
    ! check consistency
    !
    ! I: flux computation
    IF ( ANY(ihadv_tracer(1:ntracer) > 4) .OR.                    &
      &  ANY(ihadv_tracer(1:ntracer) < 0))            THEN
      CALL finish( TRIM(routine),                                       &
           'incorrect settings for ihadv_tracer. Must be 0,1,2,3, or 4 ')
    ENDIF
    SELECT CASE (i_cell_type)
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
    IF ( ANY(ivadv_tracer(1:ntracer) > ippm_v) .OR.                     &
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

    !
    ! II: limiter
    IF ( ANY(itype_vlimit(1:ntracer) < inol_v ) .OR.                    &
      &  ANY(itype_vlimit(1:ntracer) > ifluxl_vpd)) THEN
      CALL finish( TRIM(routine),                                       &
       'incorrect settings for itype_vlimit. Must be 0,1,2 or 4 ')
    ENDIF
    IF ( ANY(itype_hlimit(1:ntracer) < inol ) .OR.                      &
      &  ANY(itype_hlimit(1:ntracer) > ifluxl_sm)) THEN
      CALL finish( TRIM(routine),                                       &
       'incorrect settings for itype_hlimit. Must be 0,1,2,3 or 4 ')
    ENDIF
    IF (i_cell_type == 6) THEN
      IF ( ANY(itype_hlimit(1:ntracer) == islopel_sm ) .OR.             &
        &  ANY(itype_hlimit(1:ntracer) == islopel_m  ) .OR.             &
        &  ANY(itype_hlimit(1:ntracer) == ifluxl_m   )) THEN
        CALL finish( TRIM(routine),                                     &
         'incorrect settings for itype_hlimit and hexagonal grid. Must be 0 or 4 ')
      ENDIF
    ENDIF



    !
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


    !
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

!DR    IF ( ANY(ihadv_tracer(1:ntracer) == imiura3) ) THEN
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

!DR    END IF


  END SUBROUTINE setup_transport


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Integrates tracer continuity equation from old time step to new time step
  !!
  !! This subroutine integrates the tracer continuity equation using a simple
  !! forward time step.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-02-24)
  !!
  SUBROUTINE tupdate_tracer( p_patch, p_dtime, p_tracer_now, p_density_c_now, &
    &                        p_density_c_new, p_fluxdiv_c, p_tracer_new,      &
    &                        opt_rlstart, opt_rlend )

    TYPE(t_patch), TARGET, INTENT(IN) ::   & !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) :: p_dtime      !< time step

    REAL(wp), INTENT(IN) ::     &        !< tracer field at current time
      &  p_tracer_now(:,:,:)             !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::     &        !< density (or layer thickness) at current time
      &  p_density_c_now(:,:,:)          !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::     &        !< density (or layer thickness) at new time
      &  p_density_c_new(:,:,:)          !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::     &        !< flux divergence at current time
      &  p_fluxdiv_c(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(INOUT) ::  &        !< tracer field at current time
      &  p_tracer_new(:,:,:)             !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL :: &   !< optional: refinement control start level
     &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: &   !< optional: refinement control end level
     &  opt_rlend                        !< (to avoid calculation of halo points)

    INTEGER :: nlev                      !< number of full levels
    INTEGER :: jb, jk, jc                !< loop indices
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_rlstart, i_rlend, i_nchdom !< start and end values of refined grid
    INTEGER :: i_startidx, i_endidx
   !-----------------------------------------------------------------------

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell_int
    ENDIF

    ! number of vertical levels
    nlev = p_patch%nlev

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
          DO jb = i_startblk, i_endblk

           CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

           DO jk = 1, nlev

             DO jc = i_startidx, i_endidx

               p_tracer_new(jc,jk,jb) =                                     &
                 &   ( p_tracer_now(jc,jk,jb) * p_density_c_now(jc,jk,jb)   &
                 &    - p_dtime * p_fluxdiv_c(jc,jk,jb) )                   &
                 &    / p_density_c_new(jc,jk,jb)

             ENDDO
           ENDDO
         ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE tupdate_tracer



  !-------------------------------------------------------------------------
  !>
  !! Computation of first order backward trajectories for FFSL transport scheme
  !!
  !! Computes backward trajectories in order to determine the barycenter of the
  !! departure region. Here, a simple first order method is used. Computations are
  !! performed on a plane tangent to the edge midpoint. Coordinate axes point into
  !! the local normal and tangential direction.
  !! Once the barycenter of the departure region is known, the distance vector
  !! between the circumcenter of the upstream cell and the barycenter is computed.
  !! In a final step, this vector is transformed into a different coordinate system
  !! which has its origin at the circumcenter. The coordinate axes point to the local
  !! east and local north. Note that this subroutine has specifically been designed
  !! for the MIURA scheme with second order reconstruction of the subgrid distribution.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-03-17)
  !!
  SUBROUTINE back_traj_o1( ptr_p, ptr_int, p_vn, p_vt, p_dthalf, p_cell_indices, &
    &                     p_distv_bary, opt_rlstart, opt_rlend )

    TYPE(t_patch), TARGET, INTENT(in) ::      &  !< patch on which computation is performed
      &  ptr_p

    TYPE(t_int_state), TARGET, INTENT(in) ::  &  !< pointer to data structure for interpolation
      &  ptr_int

    REAL(wp), INTENT(IN)    ::  &  !< normal component of velocity vector at edge midpoints
      &  p_vn(:,:,:)

    REAL(wp), INTENT(IN)    ::  &  !< tangential component of velocity vector at
      &  p_vt(:,:,:)               !< edge midpoints

    REAL(wp), INTENT(IN)    ::  &  !< $0.5 \Delta t$
      &  p_dthalf

    REAL(wp), INTENT(OUT)   ::  &  !< distance vectors cell center --> barycenter of
      &  p_distv_bary(:,:,:,:)     !< departure region. (geographical coordinates)
                                   !< dim: (nproma,nlev,ptr_p%nblks_e,2)

    INTEGER, INTENT(OUT)    ::  &  !< line and block indices of cell centers in which the
      &  p_cell_indices(:,:,:,:)   !< calculated barycenters are located.
                                   !< dim: (nproma,nlev,ptr_p%nblks_e,2)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    REAL(wp) :: pos_barycenter(2),  &   !< position of barycenter and distance vector
      &         z_ntdistv_bary(2)       !< cell center --> barycenter in 'normal' and
                                        !< 'tangential' coordinates.

    INTEGER :: je, jk, jb        !< index of edge, vert level, block
    INTEGER :: nlev              !< number of full levels
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart, i_rlend, i_nchdom

  !-------------------------------------------------------------------------

    ! number of vertical levels
    nlev   = ptr_p%nlev

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    ! number of child domains
    i_nchdom   = MAX(1,ptr_p%n_childdom)

    i_startblk = ptr_p%edges%start_blk(i_rlstart,1)
    i_endblk   = ptr_p%edges%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,pos_barycenter,z_ntdistv_bary)
    DO jb = i_startblk, i_endblk

     CALL get_indices_e(ptr_p, jb, i_startblk, i_endblk,        &
                        i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, nlev
        DO je = i_startidx, i_endidx

          !
          ! Calculate backward trajectories
          !

          ! position of barycenter in normal direction
          pos_barycenter(1) = - p_vn(je,jk,jb) * p_dthalf

          ! position of barycenter in tangential direction
          pos_barycenter(2) = - p_vt(je,jk,jb) * p_dthalf


          ! Determine the cell in which the barycenter is located (upwind cell)
          ! Cell indices are chosen such that the direction from cell 1 to cell 2
          ! is the positive direction of the normal vector N, therefore:
          IF ( p_vn(je,jk,jb) >= 0._wp ) THEN

            !! we are in cell 1 !!

            ! line and block indices of neighbor cell with barycenter
            p_cell_indices(je,jk,jb,1) = ptr_p%edges%cell_idx(je,jb,1)
            p_cell_indices(je,jk,jb,2) = ptr_p%edges%cell_blk(je,jb,1)


            ! Calculate the distance cell center --> barycenter for the cell,
            ! in which the barycenter is located. The distance vector points
            ! from the cell center to the barycenter.
            z_ntdistv_bary(1:2) = pos_barycenter(1:2)                      &
              &                 - ptr_int%pos_on_tplane_e(je,jb,1,1:2)


            ! In a last step, transform this distance vector into a rotated
            ! geographical coordinate system with its origin at the circumcenter
            ! of the upstream cell. Coordinate axes point to local East and local
            ! North.

            ! component in longitudinal direction
            p_distv_bary(je,jk,jb,1) =                                             &
              &    z_ntdistv_bary(1) * ptr_p%edges%primal_normal_cell(je,jb,1)%v1  &
              &  + z_ntdistv_bary(2) * ptr_p%edges%dual_normal_cell(je,jb,1)%v1

            ! component in latitudinal direction
            p_distv_bary(je,jk,jb,2) =                                             &
              &    z_ntdistv_bary(1) * ptr_p%edges%primal_normal_cell(je,jb,1)%v2  &
              &  + z_ntdistv_bary(2) * ptr_p%edges%dual_normal_cell(je,jb,1)%v2


          ELSE

            !! we are in cell 2 !!

            ! line and block indices of neighboring cell with barycenter
            p_cell_indices(je,jk,jb,1) = ptr_p%edges%cell_idx(je,jb,2)
            p_cell_indices(je,jk,jb,2) = ptr_p%edges%cell_blk(je,jb,2)


            ! Calculate the distance cell center --> barycenter for the cell,
            ! in which the barycenter is located. The distance vector points
            ! from the cell center to the barycenter.
            z_ntdistv_bary(1:2) = pos_barycenter(1:2)                      &
              &                 - ptr_int%pos_on_tplane_e(je,jb,2,1:2)


            ! In a last step, transform this distance vector into a rotated
            ! geographical coordinate system with its origin at the circumcenter
            ! of the upstream cell. Coordinate axes point to local East and local
            ! North.

            ! component in longitudinal direction
            p_distv_bary(je,jk,jb,1) =                                             &
              &    z_ntdistv_bary(1) * ptr_p%edges%primal_normal_cell(je,jb,2)%v1  &
              &  + z_ntdistv_bary(2) * ptr_p%edges%dual_normal_cell(je,jb,2)%v1

            ! component in latitudinal direction
            p_distv_bary(je,jk,jb,2) =                                             &
              &    z_ntdistv_bary(1) * ptr_p%edges%primal_normal_cell(je,jb,2)%v2  &
              &  + z_ntdistv_bary(2) * ptr_p%edges%dual_normal_cell(je,jb,2)%v2

          ENDIF


        ENDDO ! loop over edges
      ENDDO   ! loop over vertical levels
    END DO    ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE back_traj_o1



  !-------------------------------------------------------------------------
  !>
  !! Computation of first order backward trajectories for FFSL transport scheme
  !!
  !! Computes backward trajectories in order to determine an approximation to the
  !! departure region. Here, the departure region is approximated by a rhomboidal
  !! region, using a simple first order accurate trajectory computation. Computations
  !! are performed on a plane tangent to the edge midpoint. Coordinate axes point into
  !! the local normal and tangential direction.
  !! Once the 4 vertices of the departure region are known, the distance vector
  !! between the circumcenter of the upstream cell and the vertices is computed.
  !! In a final step, these vectors are transformed into a different coordinate system
  !! which has its origin at the circumcenter. The coordinate axes point to the local
  !! east and local north. This subroutine may be combined with any reconstruction method
  !! for the subgrid distribution.
  !! Note: Care has to be taken that the vertices of the departure region are stored in
  !! counterclockwise order. This is a requirement of the gaussian quadrature which
  !! follows lateron.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-05-12)
  !!
  SUBROUTINE back_traj_dreg_o1( ptr_p, ptr_int, p_vn, p_vt, p_dt, p_cell_indices, &
    &                     p_coords_dreg_v, opt_rlstart, opt_rlend )

    TYPE(t_patch), TARGET, INTENT(IN) ::     &  !< patch on which computation is performed
      &  ptr_p

    TYPE(t_int_state), TARGET, INTENT(IN) :: &  !< pointer to data structure for interpolation
      &  ptr_int

    REAL(wp), INTENT(IN)    ::  &  !< normal component of velocity vector at edge midpoints
      &  p_vn(:,:,:)

    REAL(wp), INTENT(IN)    ::  &  !< tangential component of velocity vector at
      &  p_vt(:,:,:)               !< edge midpoints

    REAL(wp), INTENT(IN)    ::  &  !< time step $\Delta t$
      &  p_dt

    REAL(wp), INTENT(OUT) ::    &  !< coordinates of departure region vertices. The origin
      &  p_coords_dreg_v(:,:,:,:,:)!< of the coordinate system is at the circumcenter of
                                   !< the upwind cell. Unit vectors point to local East
                                   !< and North. (geographical coordinates)
                                   !< dim: (nproma,nlev,ptr_p%nblks_e,4,2)

    INTEGER, INTENT(OUT)  ::    &  !< line and block indices of upwind cell
      &  p_cell_indices(:,:,:,:)   !< dim: (nproma,nlev,ptr_p%nblks_e,2)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    REAL(wp), POINTER  ::  &       !< pointer to coordinates of edge vertices
      &  ptr_ve(:,:,:,:)           !< on tangent plane

    REAL(wp) ::            &       !< coordinates of departure region vertices
      &  pos_dreg_vert_e(4,2)      !< in edge-based coordinate system

    REAL(wp) ::            &       !< coordinates of departure region vertices
      &  pos_dreg_vert_c(4,2)      !< as seen from translated coordinate system.
                                   !< origin at circumcenter of upwind cell

    INTEGER :: je, jk, jb          !< index of edge, vert level, block
    INTEGER :: nlev                !< number of full levels
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart, i_rlend, i_nchdom
  !-------------------------------------------------------------------------

    ! number of vertical levels
    nlev = ptr_p%nlev

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    ! number of child domains
    i_nchdom   = MAX(1,ptr_p%n_childdom)

    i_startblk = ptr_p%edges%start_blk(i_rlstart,1)
    i_endblk   = ptr_p%edges%end_blk(i_rlend,i_nchdom)

    ! pointer to the edge vertices
    ptr_ve => ptr_int%pos_on_tplane_e(:,:,7:8,:)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,pos_dreg_vert_e, &
!$OMP            pos_dreg_vert_c)
    DO jb = i_startblk, i_endblk

     CALL get_indices_e(ptr_p, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, nlev
        DO je = i_startidx, i_endidx



          ! departure region and correct counterclockwise numbering of vertices
          !--------------------------------------------------------------------
          !
          !  3\--------------\2       <- correct counterclockwise numbering
          !    \              \
          !     \      N       \      <- edge normal
          !      \      \       \
          !     4 \------\-------\1   <- triangle edge


          ! Determine the upwind cell
          ! Cell indices are chosen such that the direction from cell 1 to cell 2
          ! is the positive direction of the normal vector N.

          !
          ! Calculate backward trajectories, starting at the two edge vertices
          ! It is assumed the velocity vector is constant along the edge.
          !
          IF ( p_vn(je,jk,jb) >= 0._wp ) THEN

            !! we are in cell 1 !!

            ! Vertices of the departure cell
            ! Take care of correct counterclockwise numbering
            !
            ! position of vertex 4 in normal direction
            pos_dreg_vert_e(4,1) = ptr_ve(je,jb,1,1) - p_vn(je,jk,jb) * p_dt

            ! position of vertex 4 in tangential direction
            pos_dreg_vert_e(4,2) = ptr_ve(je,jb,1,2) - p_vt(je,jk,jb) * p_dt

            ! position of vertex 3 in normal direction
            pos_dreg_vert_e(3,1) = ptr_ve(je,jb,2,1) - p_vn(je,jk,jb) * p_dt

            ! position of vertex 3 in tangential direction
            pos_dreg_vert_e(3,2) = ptr_ve(je,jb,2,2) - p_vt(je,jk,jb) * p_dt

            ! position of vertex 1
            pos_dreg_vert_e(1,1:2) = ptr_ve(je,jb,1,1:2)

            ! position of vertex 2
            pos_dreg_vert_e(2,1:2) = ptr_ve(je,jb,2,1:2)


            ! line and block indices of upwind cell
            p_cell_indices(je,jk,jb,1) = ptr_p%edges%cell_idx(je,jb,1)
            p_cell_indices(je,jk,jb,2) = ptr_p%edges%cell_blk(je,jb,1)



            ! Calculate position of departure region vertices in a translated
            ! coordinate system. The origin is located at the circumcenter
            ! of the upwind cell. The distance vectors point from the cell center
            ! to the vertices.
            pos_dreg_vert_c(1,1:2) = pos_dreg_vert_e(1,1:2)                &
              &                     - ptr_int%pos_on_tplane_e(je,jb,1,1:2)
            pos_dreg_vert_c(2,1:2) = pos_dreg_vert_e(2,1:2)                &
              &                     - ptr_int%pos_on_tplane_e(je,jb,1,1:2)
            pos_dreg_vert_c(3,1:2) = pos_dreg_vert_e(3,1:2)                &
              &                     - ptr_int%pos_on_tplane_e(je,jb,1,1:2)
            pos_dreg_vert_c(4,1:2) = pos_dreg_vert_e(4,1:2)                &
              &                     - ptr_int%pos_on_tplane_e(je,jb,1,1:2)


            ! In a last step, these distance vectors are transformed into a rotated
            ! geographical coordinate system, which still has its origin at the circumcenter
            ! of the upwind cell. Now the coordinate axes point to local East and local
            ! North.

            ! components in longitudinal direction
            p_coords_dreg_v(je,jk,jb,1,1) =                                          &
              &    pos_dreg_vert_c(1,1) * ptr_p%edges%primal_normal_cell(je,jb,1)%v1 &
              &  + pos_dreg_vert_c(1,2) * ptr_p%edges%dual_normal_cell(je,jb,1)%v1

            p_coords_dreg_v(je,jk,jb,2,1) =                                          &
              &    pos_dreg_vert_c(2,1) * ptr_p%edges%primal_normal_cell(je,jb,1)%v1 &
              &  + pos_dreg_vert_c(2,2) * ptr_p%edges%dual_normal_cell(je,jb,1)%v1

            p_coords_dreg_v(je,jk,jb,3,1) =                                          &
              &    pos_dreg_vert_c(3,1) * ptr_p%edges%primal_normal_cell(je,jb,1)%v1 &
              &  + pos_dreg_vert_c(3,2) * ptr_p%edges%dual_normal_cell(je,jb,1)%v1

            p_coords_dreg_v(je,jk,jb,4,1) =                                          &
              &    pos_dreg_vert_c(4,1) * ptr_p%edges%primal_normal_cell(je,jb,1)%v1 &
              &  + pos_dreg_vert_c(4,2) * ptr_p%edges%dual_normal_cell(je,jb,1)%v1


            ! components in latitudinal direction
            p_coords_dreg_v(je,jk,jb,1,2) =                                          &
              &    pos_dreg_vert_c(1,1) * ptr_p%edges%primal_normal_cell(je,jb,1)%v2 &
              &  + pos_dreg_vert_c(1,2) * ptr_p%edges%dual_normal_cell(je,jb,1)%v2

            p_coords_dreg_v(je,jk,jb,2,2) =                                          &
              &    pos_dreg_vert_c(2,1) * ptr_p%edges%primal_normal_cell(je,jb,1)%v2 &
              &  + pos_dreg_vert_c(2,2) * ptr_p%edges%dual_normal_cell(je,jb,1)%v2

            p_coords_dreg_v(je,jk,jb,3,2) =                                          &
              &    pos_dreg_vert_c(3,1) * ptr_p%edges%primal_normal_cell(je,jb,1)%v2 &
              &  + pos_dreg_vert_c(3,2) * ptr_p%edges%dual_normal_cell(je,jb,1)%v2

            p_coords_dreg_v(je,jk,jb,4,2) =                                          &
              &    pos_dreg_vert_c(4,1) * ptr_p%edges%primal_normal_cell(je,jb,1)%v2 &
              &  + pos_dreg_vert_c(4,2) * ptr_p%edges%dual_normal_cell(je,jb,1)%v2


          ELSE

            !! we are in cell 2 !!

            ! Vertices of the departure cell
            ! Take care of the correct counterclockwise numbering
            !
            ! position of vertex 4 in normal direction
            pos_dreg_vert_e(2,1) = ptr_ve(je,jb,1,1) - p_vn(je,jk,jb) * p_dt

            ! position of vertex 4 in tangential direction
            pos_dreg_vert_e(2,2) = ptr_ve(je,jb,1,2) - p_vt(je,jk,jb) * p_dt

            ! position of vertex 3 in normal direction
            pos_dreg_vert_e(3,1) = ptr_ve(je,jb,2,1) - p_vn(je,jk,jb) * p_dt

            ! position of vertex 3 in tangential direction
            pos_dreg_vert_e(3,2) = ptr_ve(je,jb,2,2) - p_vt(je,jk,jb) * p_dt

            ! position of vertex 1
            pos_dreg_vert_e(1,1:2) = ptr_ve(je,jb,1,1:2)

            ! position of vertex 2
            pos_dreg_vert_e(4,1:2) = ptr_ve(je,jb,2,1:2)



            ! line and block indices of upwind cell
            p_cell_indices(je,jk,jb,1) = ptr_p%edges%cell_idx(je,jb,2)
            p_cell_indices(je,jk,jb,2) = ptr_p%edges%cell_blk(je,jb,2)



            ! Calculate position of departure region vertices in a translated
            ! coordinate system. The origin is located at the circumcenter
            ! of the upwind cell. The distance vectors point from the cell center
            ! to the vertices.
            pos_dreg_vert_c(1,1:2) = pos_dreg_vert_e(1,1:2)                &
              &                     - ptr_int%pos_on_tplane_e(je,jb,2,1:2)
            pos_dreg_vert_c(2,1:2) = pos_dreg_vert_e(2,1:2)                &
              &                     - ptr_int%pos_on_tplane_e(je,jb,2,1:2)
            pos_dreg_vert_c(3,1:2) = pos_dreg_vert_e(3,1:2)                &
              &                     - ptr_int%pos_on_tplane_e(je,jb,2,1:2)
            pos_dreg_vert_c(4,1:2) = pos_dreg_vert_e(4,1:2)                &
              &                     - ptr_int%pos_on_tplane_e(je,jb,2,1:2)


            ! In a last step, these distance vectors are transformed into a rotated
            ! geographical coordinate system, which still has its origin at the circumcenter
            ! of the upwind cell. But the coordinate axes now point to local East and local
            ! North.

            ! components in longitudinal direction
            p_coords_dreg_v(je,jk,jb,1,1) =                                          &
              &    pos_dreg_vert_c(1,1) * ptr_p%edges%primal_normal_cell(je,jb,2)%v1 &
              &  + pos_dreg_vert_c(1,2) * ptr_p%edges%dual_normal_cell(je,jb,2)%v1

            p_coords_dreg_v(je,jk,jb,2,1) =                                          &
              &    pos_dreg_vert_c(2,1) * ptr_p%edges%primal_normal_cell(je,jb,2)%v1 &
              &  + pos_dreg_vert_c(2,2) * ptr_p%edges%dual_normal_cell(je,jb,2)%v1

            p_coords_dreg_v(je,jk,jb,3,1) =                                          &
              &    pos_dreg_vert_c(3,1) * ptr_p%edges%primal_normal_cell(je,jb,2)%v1 &
              &  + pos_dreg_vert_c(3,2) * ptr_p%edges%dual_normal_cell(je,jb,2)%v1

            p_coords_dreg_v(je,jk,jb,4,1) =                                          &
              &    pos_dreg_vert_c(4,1) * ptr_p%edges%primal_normal_cell(je,jb,2)%v1 &
              &  + pos_dreg_vert_c(4,2) * ptr_p%edges%dual_normal_cell(je,jb,2)%v1


            ! components in latitudinal direction
            p_coords_dreg_v(je,jk,jb,1,2) =                                          &
              &    pos_dreg_vert_c(1,1) * ptr_p%edges%primal_normal_cell(je,jb,2)%v2 &
              &  + pos_dreg_vert_c(1,2) * ptr_p%edges%dual_normal_cell(je,jb,2)%v2

            p_coords_dreg_v(je,jk,jb,2,2) =                                          &
              &    pos_dreg_vert_c(2,1) * ptr_p%edges%primal_normal_cell(je,jb,2)%v2 &
              &  + pos_dreg_vert_c(2,2) * ptr_p%edges%dual_normal_cell(je,jb,2)%v2

            p_coords_dreg_v(je,jk,jb,3,2) =                                          &
              &    pos_dreg_vert_c(3,1) * ptr_p%edges%primal_normal_cell(je,jb,2)%v2 &
              &  + pos_dreg_vert_c(3,2) * ptr_p%edges%dual_normal_cell(je,jb,2)%v2

            p_coords_dreg_v(je,jk,jb,4,2) =                                          &
              &    pos_dreg_vert_c(4,1) * ptr_p%edges%primal_normal_cell(je,jb,2)%v2 &
              &  + pos_dreg_vert_c(4,2) * ptr_p%edges%dual_normal_cell(je,jb,2)%v2

          ENDIF


        ENDDO ! loop over edges
      ENDDO   ! loop over vertical levels
    END DO    ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE back_traj_dreg_o1



  !-------------------------------------------------------------------------
  !>
  !! Computation of second order backward trajectories for FFSL transport scheme
  !!
  !! Computes backward trajectories in order to determine the barycenter of the
  !! departure region. Here, an iterative second order method is used. Computations
  !! are performed on a plane tangent to the edge midpoint. Coordinate axes point
  !! into the local normal and tangential direction. A bilinear interpolation is
  !! applied using the velocity vectors at the three edge midpoints of the upwind cell
  !! in order to derive an improved estimate of the velocity.
  !!
  !! Once the barycenter of the departure region is known, the distance vector
  !! between the circumcenter of the upstream cell and the barycenter is computed.
  !! In a final step, this vector is transformed into a different coordinate system
  !! which has its origin at the circumcenter. The coordinate axes point to the local
  !! east and local north.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-03-24)
  !!
  SUBROUTINE back_traj_o2( ptr_p, ptr_int, p_vn, p_vt, p_dthalf, p_cell_indices,  &
    &                     p_distv_bary, opt_rlstart, opt_rlend )

    TYPE(t_patch), TARGET, INTENT(IN) ::      &  !< patch on which computation is performed
      &  ptr_p

    TYPE(t_int_state), TARGET, INTENT(IN) ::  &  !< pointer to data structure for interpolation
      &  ptr_int

    REAL(wp), INTENT(IN)    ::  &  !< normal component of velocity vector at edge midpoints
      &  p_vn(:,:,:)

    REAL(wp), INTENT(IN)    ::  &  !< tangential component of velocity vector at
      &  p_vt(:,:,:)               !< edge midpoints

    REAL(wp), INTENT(IN)  ::    &  !< $0.5 \Delta t$
      &  p_dthalf

    REAL(wp), INTENT(OUT) ::    &  !< distance vectors cell center --> barycenter of
      &  p_distv_bary(:,:,:,:)     !< departure region. (geographical coordinates)
                                   !< dim: (nproma,nlev,ptr_p%nblks_e,2)

    INTEGER, INTENT(OUT)  ::    &  !< line and block indices of cell centers in which the
      &  p_cell_indices(:,:,:,:)   !< computed barycenters are located.
                                   !< dim: (nproma,nlev,ptr_p%nblks_e,2)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    REAL(wp) :: pos_barycenter(2),  &   !< position of barycenter and distance vector
      &         z_ntdistv_bary(2)       !< cell center --> barycenter in 'normal' and
                                        !< 'tangential' coordinates.

    REAL(wp), POINTER  ::  &     !< pointer to positions of quadrilateral edge midpoints
      &  ptr_em(:,:,:,:)         !< on tangent plane

    REAL(wp)           ::  &     ! normal and tangential components of the velocity vectors
      &  z_vn_plane(nproma,ptr_p%nlev,ptr_p%nblks_e,4), &!< projected into local edge-based
      &  z_vt_plane(nproma,ptr_p%nlev,ptr_p%nblks_e,4)   !< coordinate system


    REAL(wp)  :: w1, w2, w3      !< weights for bilinear interpolation

    REAL(wp)  ::           &     !< resulting normal and tangential velocity component
      &  vn_new, vt_new          !< after bilinear interpolation

    INTEGER :: je, jk, jb        !< index of edge, vert level, block
    INTEGER :: nlev              !< number of full levels
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart, i_rlend, i_nchdom
    INTEGER :: nblks_e, npromz_e
    INTEGER :: zcell             !< determines whether the barycenter is located
                                 !< in cell 1 or 2
    INTEGER, POINTER ::    &     !< pointer for line and block indices of edge
      & iidx(:,:,:), iblk(:,:,:) !< midpoints for quadrilateral cell

!DR    REAL(wp) :: z_vabs_orig, z_vabs_new

  !-------------------------------------------------------------------------

    ! number of vertical levels
    nlev = ptr_p%nlev

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 1
    ENDIF

    ! number of child domains
    i_nchdom   = MAX(1,ptr_p%n_childdom)

    i_startblk = ptr_p%edges%start_blk(i_rlstart,1)
    i_endblk   = ptr_p%edges%end_blk(i_rlend,i_nchdom)


    ! For each edge the velocity vectors at the quadrilateral edges are
    ! transformed into the new edge-based coordinate system (coordinate axes
    ! normal and tangential to the inner edge).
    ! - projection of normal component into local normal direction
    ! - projection of tangential component into local normal direction
    ! - projection of normal component into local tangential direction
    ! - projection of tangential component into local tangential direction

    ! pointer to the quadrilateral edge midpoints
    ! Note that ptr_em(:,:,1:2,:) belong to cell number 1 and
    ! ptr_em(:,:,3:4,:) belong to cell number 2 as seen from the edge.
    ptr_em => ptr_int%pos_on_tplane_e(:,:,3:6,:)

    ! pointer to line and block indices of outer quadrilateral edge midpoints
    iidx => ptr_p%edges%quad_idx
    iblk => ptr_p%edges%quad_blk

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_p, jb, i_startblk, i_endblk,        &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif

            z_vn_plane(je,jk,jb,1) =                                                         &
              &   p_vn(iidx(je,jb,1),jk,iblk(je,jb,1)) * ptr_int%tplane_e_dotprod(je,jb,1,1) &
              & + p_vt(iidx(je,jb,1),jk,iblk(je,jb,1)) * ptr_int%tplane_e_dotprod(je,jb,1,2)

            z_vn_plane(je,jk,jb,2) =                                                         &
              &   p_vn(iidx(je,jb,2),jk,iblk(je,jb,2)) * ptr_int%tplane_e_dotprod(je,jb,2,1) &
              & + p_vt(iidx(je,jb,2),jk,iblk(je,jb,2)) * ptr_int%tplane_e_dotprod(je,jb,2,2)

            z_vn_plane(je,jk,jb,3) =                                                         &
              &   p_vn(iidx(je,jb,3),jk,iblk(je,jb,3)) * ptr_int%tplane_e_dotprod(je,jb,3,1) &
              & + p_vt(iidx(je,jb,3),jk,iblk(je,jb,3)) * ptr_int%tplane_e_dotprod(je,jb,3,2)

            z_vn_plane(je,jk,jb,4) =                                                         &
              &   p_vn(iidx(je,jb,4),jk,iblk(je,jb,4)) * ptr_int%tplane_e_dotprod(je,jb,4,1) &
              & + p_vt(iidx(je,jb,4),jk,iblk(je,jb,4)) * ptr_int%tplane_e_dotprod(je,jb,4,2)


            z_vt_plane(je,jk,jb,1) =                                                         &
              &   p_vn(iidx(je,jb,1),jk,iblk(je,jb,1)) * ptr_int%tplane_e_dotprod(je,jb,1,3) &
              & + p_vt(iidx(je,jb,1),jk,iblk(je,jb,1)) * ptr_int%tplane_e_dotprod(je,jb,1,4)

            z_vt_plane(je,jk,jb,2) =                                                         &
              &   p_vn(iidx(je,jb,2),jk,iblk(je,jb,2)) * ptr_int%tplane_e_dotprod(je,jb,2,3) &
              & + p_vt(iidx(je,jb,2),jk,iblk(je,jb,2)) * ptr_int%tplane_e_dotprod(je,jb,2,4)

            z_vt_plane(je,jk,jb,3) =                                                         &
              &   p_vn(iidx(je,jb,3),jk,iblk(je,jb,3)) * ptr_int%tplane_e_dotprod(je,jb,3,3) &
              & + p_vt(iidx(je,jb,3),jk,iblk(je,jb,3)) * ptr_int%tplane_e_dotprod(je,jb,3,4)

            z_vt_plane(je,jk,jb,4) =                                                         &
              &   p_vn(iidx(je,jb,4),jk,iblk(je,jb,4)) * ptr_int%tplane_e_dotprod(je,jb,4,3) &
              & + p_vt(iidx(je,jb,4),jk,iblk(je,jb,4)) * ptr_int%tplane_e_dotprod(je,jb,4,4)


!            ! Re-normalize the projected velocity vector
!            z_vabs_orig = SQRT( p_vn(ilq,jk,ibq)*p_vn(ilq,jk,ibq)                &
!              &         + p_vt(ilq,jk,ibq)*p_vt(ilq,jk,ibq) )
!            z_vabs_new  = SQRT( z_vn_plane(je,jk,jb,ne)*z_vn_plane(je,jk,jb,ne)  &
!              &         + z_vt_plane(je,jk,jb,ne)*z_vt_plane(je,jk,jb,ne) )

!            z_vn_plane(je,jk,jb,ne) = z_vn_plane(je,jk,jb,ne) * (z_vabs_orig/z_vabs_new)
!            z_vt_plane(je,jk,jb,ne) = z_vt_plane(je,jk,jb,ne) * (z_vabs_orig/z_vabs_new)


        ENDDO ! loop over edges
      ENDDO   ! loop over vertical levels
    END DO    ! loop over blocks
!$OMP END DO



!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,pos_barycenter,zcell,w1,w2,w3, &
!$OMP            vn_new,vt_new,z_ntdistv_bary)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_p, jb, i_startblk, i_endblk,        &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, nlev
        DO je = i_startidx, i_endidx


          !
          ! Calculate backward trajectories
          !

          ! First guess
          ! position of barycenter in normal direction
          pos_barycenter(1) = - p_vn(je,jk,jb) * p_dthalf

          ! position of barycenter in tangential direction
          pos_barycenter(2) = - p_vt(je,jk,jb) * p_dthalf



          IF (p_vn(je,jk,jb) >= 0._wp) THEN

            !! we are in cell 1 !!

            ! line and block indices of neighboring cell with barycenter
            p_cell_indices(je,jk,jb,1) = ptr_p%edges%cell_idx(je,jb,1)
            p_cell_indices(je,jk,jb,2) = ptr_p%edges%cell_blk(je,jb,1)

            zcell = 1

            ! calculate weights for bilinear interpolation of velocities onto the
            ! barycenter
            w3 = ( pos_barycenter(2)                                          &
              &   - pos_barycenter(1) * ptr_em(je,jb,1,2)/ptr_em(je,jb,1,1) ) &
              &   / ( ptr_em(je,jb,2,2)                                       &
              &   - ptr_em(je,jb,2,1) * ptr_em(je,jb,1,2)/ptr_em(je,jb,1,1) )

            w2 = (pos_barycenter(1) - w3 * ptr_em(je,jb,2,1))       &
              &  / ptr_em(je,jb,1,1)

            w1 = 1._wp - w2 - w3


            ! calculate improved normal and tangential velocity components using
            ! the bilinear interpolation weights.
            vn_new = w1 * p_vn(je,jk,jb) + w2 * z_vn_plane(je,jk,jb,1) &
              &    + w3 * z_vn_plane(je,jk,jb,2)
            vt_new = w1 * p_vt(je,jk,jb) + w2 * z_vt_plane(je,jk,jb,1) &
              &    + w3 * z_vt_plane(je,jk,jb,2)

          ELSE

            !! we are in cell 2 !!

            ! line and block indices of neighboring cell with barycenter
            p_cell_indices(je,jk,jb,1) = ptr_p%edges%cell_idx(je,jb,2)
            p_cell_indices(je,jk,jb,2) = ptr_p%edges%cell_blk(je,jb,2)

            zcell = 2


            ! calculate weights for bilinear interpolation of velocities onto the
            ! barycenter
            w3 = ( pos_barycenter(2)                                          &
              &   - pos_barycenter(1) * ptr_em(je,jb,3,2)/ptr_em(je,jb,3,1) ) &
              &   / ( ptr_em(je,jb,4,2)                                       &
              &   - ptr_em(je,jb,4,1) * ptr_em(je,jb,3,2)/ptr_em(je,jb,3,1) )

            w2 = (pos_barycenter(1) - w3 * ptr_em(je,jb,4,1))       &
              &  / ptr_em(je,jb,3,1)

            w1 = 1._wp - w2 - w3


            ! calculate improved normal and tangential velocity components using
            ! the bilinear interpolation weights.
            vn_new = w1 * p_vn(je,jk,jb) + w2 * z_vn_plane(je,jk,jb,3) &
              &    + w3 * z_vn_plane(je,jk,jb,4)
            vt_new = w1 * p_vt(je,jk,jb) + w2 * z_vt_plane(je,jk,jb,3) &
              &    + w3 * z_vt_plane(je,jk,jb,4)

          ENDIF


          ! Improved guess
          ! position of barycenter in normal direction
          pos_barycenter(1) = - vn_new * p_dthalf

          ! position of barycenter in tangential direction
          pos_barycenter(2) = - vt_new * p_dthalf


          ! Calculate the distance cell center --> barycenter for the cell,
          ! in which the barycenter is located. The distance vector points
          ! from the cell center to the barycenter.
          z_ntdistv_bary(1:2) = pos_barycenter(1:2)                      &
            &                 - ptr_int%pos_on_tplane_e(je,jb,zcell,1:2)


          ! In a last step, transform this distance vector into a rotated
          ! geographical coordinate system with its origin at the circumcenter
          ! of the upstream cell. Coordinate axes point to local East and local
          ! North.

          ! component in longitudinal direction
          p_distv_bary(je,jk,jb,1) =                                                 &
            &    z_ntdistv_bary(1) * ptr_p%edges%primal_normal_cell(je,jb,zcell)%v1  &
            &  + z_ntdistv_bary(2) * ptr_p%edges%dual_normal_cell(je,jb,zcell)%v1

          ! component in latitudinal direction
          p_distv_bary(je,jk,jb,2) =                                                 &
            &    z_ntdistv_bary(1) * ptr_p%edges%primal_normal_cell(je,jb,zcell)%v2  &
            &  + z_ntdistv_bary(2) * ptr_p%edges%dual_normal_cell(je,jb,zcell)%v2



        ENDDO ! loop over edges
      ENDDO   ! loop over vertical levels
    END DO    ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE back_traj_o2




  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Lax Friedrichs first order upwind flux,.
  !!
  !! Lax Friedrichs first order upwind flux,
  !! used in conservative advection routines.
  !! For passive advection, equivalent to
  !! any other first order upwind flux.
  !!
  !! @par Revision History
  !! Developed  by L.Bonaventura  (2004).
  !!
  FUNCTION laxfr_upflux( p_vn, p_psi1, p_psi2 )  RESULT(p_upflux)
    !

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: p_vn
    REAL(wp), INTENT(in) :: p_psi1, p_psi2

    REAL(wp) :: p_upflux

    !-----------------------------------------------------------------------
    p_upflux = 0.5_wp * (        p_vn  *( p_psi1 + p_psi2 )    &
      &                   - ABS( p_vn )*( p_psi2 - p_psi1 ) )

  END FUNCTION laxfr_upflux



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Lax Friedrichs first order upwind flux for vertical advection,.
  !!
  !! Generalized Lax Friedrichs first order upwind flux,
  !! used in conservative vertical advection routines.
  !! For passive advection, equivalent to any other first
  !! order upwind flux.
  !! Applicable to both pressure based and height based vertical
  !! coordinate systems. Depending on the coordinate system chosen,
  !! the sign of the second term in the flux equation changes.
  !! - (-) for pressure based vertical coordinate systems
  !! - (+) for height based coordinate systems
  !! In order to get the correct sign, the variable p_coeff_grid
  !! has been introduced which is =1 for pressure based and =-1
  !! for height based coordinate systems.
  !!
  !! @par Revision History
  !! Developed  by L.Bonaventura  (2004).
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized for p- and z-based vertical coordinate systems
  !!
  FUNCTION laxfr_upflux_v( p_vn, p_psi1, p_psi2, p_coeff_grid )  RESULT(p_upflux)
    !

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: p_vn
    REAL(wp), INTENT(in) :: p_psi1, p_psi2
    REAL(wp), INTENT(in) :: p_coeff_grid

    REAL(wp) :: p_upflux

    !-----------------------------------------------------------------------
    p_upflux = 0.5_wp * (                       p_vn  *( p_psi1 + p_psi2 )    &
      &                   - p_coeff_grid * ABS( p_vn )*( p_psi2 - p_psi1 ) )

  END FUNCTION laxfr_upflux_v




  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of quadratic tracer subgrid distribution
  !!
  !! Provides tracer-independent parts for a gauss quadrature of order 2.
  !! I.e. the 4 quadrature points in physical space and the product of weights
  !! and the determinant of the Jacobian for each quadrature point.
  !! This subroutine is specific to a reconstruction based on a quadratic
  !! polynomial. It needs to be called only once per time step, independent
  !! of the number of advected fields.
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-05-14)
  !!
  !!
  SUBROUTINE prep_gauss_quadrature_q( p_patch, p_coords_dreg_v,        &
    &                                 p_quad_vector_sum, p_rdreg_area, &
    &                                 opt_rlstart, opt_rlend )

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                           !< performed

    REAL(wp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:,:) !< in 2D cartesian coordinates
                                    !< dim: (nproma,nlev,nblks_e,4,2)

    REAL(wp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:) !< dim: (nproma,6,nlev,nblks_e)

    REAL(wp), INTENT(OUT) :: &      !< reciprocal area of departure region
      &  p_rdreg_area(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

   ! local variables
    REAL(wp) ::                &    !< coordinates of gaussian quadrature points
      &  z_gauss_pts(4,2)           !< in physical space

    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac(4)            !< each gaussian quadrature point.

    REAL(wp) ::                &    !< quadrature vector for single integration point
      &  z_quad_vector(4,6)

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: nlev                !< number of full levels
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom

  !-----------------------------------------------------------------------

    ! number of vertical levels
    nlev = p_patch%nlev

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,jg,i_startidx,i_endidx,z_gauss_pts,wgt_t_detjac,z_quad_vector)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, nlev

        DO je = i_startidx, i_endidx

          ! get coordinates of the quadrature points in physical space (mapping)
          z_gauss_pts(1,1) = DOT_PRODUCT(shape_func(1:4,1),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(1,2) = DOT_PRODUCT(shape_func(1:4,1),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(2,1) = DOT_PRODUCT(shape_func(1:4,2),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(2,2) = DOT_PRODUCT(shape_func(1:4,2),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(3,1) = DOT_PRODUCT(shape_func(1:4,3),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(3,2) = DOT_PRODUCT(shape_func(1:4,3),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(4,1) = DOT_PRODUCT(shape_func(1:4,4),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(4,2) = DOT_PRODUCT(shape_func(1:4,4),p_coords_dreg_v(je,jk,jb,1:4,2))


          ! get Jacobian determinant for each quadrature point and multiply with
          ! corresponding weights
          ! Note: dbl_eps is added, in order to have a meaningful 'edge value' 
          ! (better: area-average) even when the integration-area tends to zero.
          wgt_t_detjac(1) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(1),eta(1)) &
            &                      * wgt_zeta(1) *  wgt_eta(1) ) + dbl_eps
          wgt_t_detjac(2) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(2),eta(2)) &
            &                      * wgt_zeta(2) *  wgt_eta(2) ) + dbl_eps
          wgt_t_detjac(3) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(3),eta(3)) &
            &                      * wgt_zeta(3) *  wgt_eta(3) ) + dbl_eps
          wgt_t_detjac(4) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(4),eta(4)) &
            &                      * wgt_zeta(4) *  wgt_eta(4) ) + dbl_eps


          ! Get quadrature vector for each integration point and multiply by
          ! corresponding wgt_t_detjac
          DO jg=1, 4
            z_quad_vector(jg,1) = wgt_t_detjac(jg)
            z_quad_vector(jg,2) = wgt_t_detjac(jg) * z_gauss_pts(jg,1)
            z_quad_vector(jg,3) = wgt_t_detjac(jg) * z_gauss_pts(jg,2)
            z_quad_vector(jg,4) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1)**2)
            z_quad_vector(jg,5) = wgt_t_detjac(jg) * (z_gauss_pts(jg,2)**2)
            z_quad_vector(jg,6) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1) * z_gauss_pts(jg,2))
          ENDDO


          ! Sum quadrature vectors over all integration points
          p_quad_vector_sum(je,1,jk,jb) = SUM(z_quad_vector(:,1))
          p_quad_vector_sum(je,2,jk,jb) = SUM(z_quad_vector(:,2))
          p_quad_vector_sum(je,3,jk,jb) = SUM(z_quad_vector(:,3))
          p_quad_vector_sum(je,4,jk,jb) = SUM(z_quad_vector(:,4))
          p_quad_vector_sum(je,5,jk,jb) = SUM(z_quad_vector(:,5))
          p_quad_vector_sum(je,6,jk,jb) = SUM(z_quad_vector(:,6))


          ! reciprocal area of departure region
          p_rdreg_area(je,jk,jb) = 1._wp/SUM(wgt_t_detjac(1:4))

        ENDDO ! loop over edges

      ENDDO  ! loop over levels

    ENDDO  ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE prep_gauss_quadrature_q


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of cubic tracer subgrid distribution
  !!
  !! Provides tracer-independent parts for a gauss quadrature of order 2.
  !! I.e. the 4 quadrature points in physical space and the product of weights
  !! and the determinant of the Jacobian for each quadrature point.
  !! This subroutine is specific to a cubic reconstruction without third order
  !! cross derivatives.
  !! It needs to be called only once per time step, independent of the number
  !! of advected fields.
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-10-13)
  !!
  !!
  SUBROUTINE prep_gauss_quadrature_cpoor( p_patch, p_coords_dreg_v,        &
    &                                     p_quad_vector_sum, p_rdreg_area, &
    &                                     opt_rlstart, opt_rlend           )

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(IN) :: & !< patch on which computation is
      &  p_patch                           !< performed

    REAL(wp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:,:) !< in 2D cartesian coordinates
                                    !< dim: (nproma,nlev,nblks_e,4,2)

    REAL(wp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:) !< dim: (nproma,8,nlev,nblks_e)

    REAL(wp), INTENT(OUT) :: &      !< reciprocal area of departure region
      &  p_rdreg_area(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

   ! local variables
    REAL(wp) ::                &    !< coordinates of gaussian quadrature points
      &  z_gauss_pts(4,2)           !< in physical space

    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac(4)            !< each gaussian quadrature point.

    REAL(wp) ::                &    !< quadrature vector for single integration point
      &  z_quad_vector(4,8)

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: nlev                !< number of full levels
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom

  !-----------------------------------------------------------------------

    ! number of vertical levels
    nlev = p_patch%nlev

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,jg,i_startidx,i_endidx,z_gauss_pts,wgt_t_detjac,z_quad_vector)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, nlev

        DO je = i_startidx, i_endidx

          ! get coordinates of the quadrature points in physical space (mapping)
          z_gauss_pts(1,1) = DOT_PRODUCT(shape_func(1:4,1),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(1,2) = DOT_PRODUCT(shape_func(1:4,1),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(2,1) = DOT_PRODUCT(shape_func(1:4,2),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(2,2) = DOT_PRODUCT(shape_func(1:4,2),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(3,1) = DOT_PRODUCT(shape_func(1:4,3),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(3,2) = DOT_PRODUCT(shape_func(1:4,3),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(4,1) = DOT_PRODUCT(shape_func(1:4,4),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(4,2) = DOT_PRODUCT(shape_func(1:4,4),p_coords_dreg_v(je,jk,jb,1:4,2))


          ! get Jacobian determinant for each quadrature point and multiply with
          ! corresponding weights
          ! Note: dbl_eps is added, in order to have a meaningful 'edge value' 
          ! (better: area-average) even when the integration-area tends to zero.
          wgt_t_detjac(1) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(1),eta(1)) &
            &                      * wgt_zeta(1) *  wgt_eta(1) ) + dbl_eps
          wgt_t_detjac(2) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(2),eta(2)) &
            &                      * wgt_zeta(1) *  wgt_eta(2) ) + dbl_eps
          wgt_t_detjac(3) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(3),eta(3)) &
            &                      * wgt_zeta(2) *  wgt_eta(1) ) + dbl_eps
          wgt_t_detjac(4) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(4),eta(4)) &
            &                      * wgt_zeta(2) *  wgt_eta(2) ) + dbl_eps


          ! Get quadrature vector for each integration point and multiply by
          ! corresponding wgt_t_detjac
          DO jg=1, 4
            z_quad_vector(jg,1) = wgt_t_detjac(jg)
            z_quad_vector(jg,2) = wgt_t_detjac(jg) * z_gauss_pts(jg,1)
            z_quad_vector(jg,3) = wgt_t_detjac(jg) * z_gauss_pts(jg,2)
            z_quad_vector(jg,4) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1)**2)
            z_quad_vector(jg,5) = wgt_t_detjac(jg) * (z_gauss_pts(jg,2)**2)
            z_quad_vector(jg,6) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1) * z_gauss_pts(jg,2))
            z_quad_vector(jg,7) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1)**3)
            z_quad_vector(jg,8) = wgt_t_detjac(jg) * (z_gauss_pts(jg,2)**3)
          ENDDO


          ! Sum quadrature vectors over all integration points
          p_quad_vector_sum(je,1,jk,jb) = SUM(z_quad_vector(:,1))
          p_quad_vector_sum(je,2,jk,jb) = SUM(z_quad_vector(:,2))
          p_quad_vector_sum(je,3,jk,jb) = SUM(z_quad_vector(:,3))
          p_quad_vector_sum(je,4,jk,jb) = SUM(z_quad_vector(:,4))
          p_quad_vector_sum(je,5,jk,jb) = SUM(z_quad_vector(:,5))
          p_quad_vector_sum(je,6,jk,jb) = SUM(z_quad_vector(:,6))
          p_quad_vector_sum(je,7,jk,jb) = SUM(z_quad_vector(:,7))
          p_quad_vector_sum(je,8,jk,jb) = SUM(z_quad_vector(:,8))


          ! reciprocal area of departure region
          p_rdreg_area(je,jk,jb) = 1._wp/SUM(wgt_t_detjac(1:4))

        ENDDO ! loop over edges

      ENDDO  ! loop over levels

    ENDDO  ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE prep_gauss_quadrature_cpoor


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of cubic tracer subgrid distribution
  !!
  !! Provides tracer-independent parts for a gauss quadrature of order 2.
  !! I.e. the 4 quadrature points in physical space and the product of weights
  !! and the determinant of the Jacobian for each quadrature point.
  !! This subroutine is specific to a cubic reconstruction.
  !! It needs to be called only once per time step, independent of the number
  !! of advected fields.
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-10-13)
  !!
  !!
  SUBROUTINE prep_gauss_quadrature_c( p_patch, p_coords_dreg_v,        &
    &                                 p_quad_vector_sum, p_rdreg_area, &
    &                                 opt_rlstart, opt_rlend           )

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                           !< performed

    REAL(wp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:,:) !< in 2D cartesian coordinates
                                    !< dim: (nproma,nlev,nblks_e,4,2)

    REAL(wp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:) !< dim: (nproma,8,nlev,nblks_e)

    REAL(wp), INTENT(OUT) :: &      !< reciprocal area of departure region
      &  p_rdreg_area(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

   ! local variables
    REAL(wp) ::                &    !< coordinates of gaussian quadrature points
      &  z_gauss_pts(4,2)           !< in physical space

    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac(4)            !< each gaussian quadrature point.

    REAL(wp) ::                &    !< quadrature vector for single integration point
      &  z_quad_vector(4,10)

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: nlev                !< number of full levels
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom

  !-----------------------------------------------------------------------

    ! number of vertical levels
    nlev = p_patch%nlev

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,jg,i_startidx,i_endidx,z_gauss_pts,wgt_t_detjac,z_quad_vector)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, nlev

        DO je = i_startidx, i_endidx

          ! get coordinates of the quadrature points in physical space (mapping)
          z_gauss_pts(1,1) = DOT_PRODUCT(shape_func(1:4,1),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(1,2) = DOT_PRODUCT(shape_func(1:4,1),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(2,1) = DOT_PRODUCT(shape_func(1:4,2),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(2,2) = DOT_PRODUCT(shape_func(1:4,2),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(3,1) = DOT_PRODUCT(shape_func(1:4,3),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(3,2) = DOT_PRODUCT(shape_func(1:4,3),p_coords_dreg_v(je,jk,jb,1:4,2))
          z_gauss_pts(4,1) = DOT_PRODUCT(shape_func(1:4,4),p_coords_dreg_v(je,jk,jb,1:4,1))
          z_gauss_pts(4,2) = DOT_PRODUCT(shape_func(1:4,4),p_coords_dreg_v(je,jk,jb,1:4,2))


          ! get Jacobian determinant for each quadrature point and multiply with
          ! corresponding weights
          ! Note: dbl_eps is added, in order to have a meaningful 'edge value' 
          ! (better: area-average) even when the integration-area tends to zero.
          wgt_t_detjac(1) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(1),eta(1)) &
            &                      * wgt_zeta(1) *  wgt_eta(1) ) + dbl_eps
          wgt_t_detjac(2) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(2),eta(2)) &
            &                      * wgt_zeta(1) *  wgt_eta(2) ) + dbl_eps
          wgt_t_detjac(3) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(3),eta(3)) &
            &                      * wgt_zeta(2) *  wgt_eta(1) ) + dbl_eps
          wgt_t_detjac(4) = ( jac(p_coords_dreg_v(je,jk,jb,1:4,1),              &
            &                   p_coords_dreg_v(je,jk,jb,1:4,2),zeta(4),eta(4)) &
            &                      * wgt_zeta(2) *  wgt_eta(2) ) + dbl_eps


          ! Get quadrature vector for each integration point and multiply by
          ! corresponding wgt_t_detjac
          DO jg=1, 4
            z_quad_vector(jg,1) = wgt_t_detjac(jg)
            z_quad_vector(jg,2) = wgt_t_detjac(jg) * z_gauss_pts(jg,1)
            z_quad_vector(jg,3) = wgt_t_detjac(jg) * z_gauss_pts(jg,2)
            z_quad_vector(jg,4) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1)**2)
            z_quad_vector(jg,5) = wgt_t_detjac(jg) * (z_gauss_pts(jg,2)**2)
            z_quad_vector(jg,6) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1) * z_gauss_pts(jg,2))
            z_quad_vector(jg,7) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1)**3)
            z_quad_vector(jg,8) = wgt_t_detjac(jg) * (z_gauss_pts(jg,2)**3)
            z_quad_vector(jg,9) = wgt_t_detjac(jg) * (z_gauss_pts(jg,1)**2 * z_gauss_pts(jg,2))
            z_quad_vector(jg,10)= wgt_t_detjac(jg) * (z_gauss_pts(jg,1) * z_gauss_pts(jg,2)**2)
          ENDDO


          ! Sum quadrature vectors over all integration points
          p_quad_vector_sum(je,1,jk,jb) = SUM(z_quad_vector(:,1))
          p_quad_vector_sum(je,2,jk,jb) = SUM(z_quad_vector(:,2))
          p_quad_vector_sum(je,3,jk,jb) = SUM(z_quad_vector(:,3))
          p_quad_vector_sum(je,4,jk,jb) = SUM(z_quad_vector(:,4))
          p_quad_vector_sum(je,5,jk,jb) = SUM(z_quad_vector(:,5))
          p_quad_vector_sum(je,6,jk,jb) = SUM(z_quad_vector(:,6))
          p_quad_vector_sum(je,7,jk,jb) = SUM(z_quad_vector(:,7))
          p_quad_vector_sum(je,8,jk,jb) = SUM(z_quad_vector(:,8))
          p_quad_vector_sum(je,9,jk,jb) = SUM(z_quad_vector(:,9))
          p_quad_vector_sum(je,10,jk,jb)= SUM(z_quad_vector(:,10))


          ! reciprocal area of departure region
          p_rdreg_area(je,jk,jb) = 1._wp/SUM(wgt_t_detjac(1:4))

        ENDDO ! loop over edges

      ENDDO  ! loop over levels

    ENDDO  ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE prep_gauss_quadrature_c



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Computes Jacobian determinant for gaussian quadrature
  !!
  !! Computes Jacobian determinant for gaussian quadrature
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-05-14)
  !!
  !!
  FUNCTION jac(x, y, zeta, eta)  RESULT(det_jac)

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: x(1:4), y(1:4)  !< coordinates of vertices in x-y-system
    REAL(wp), INTENT(IN) :: zeta, eta       !< integration point in \zeta,\eta-system

    ! RETURN VALUE:
    REAL(wp) :: det_jac

    REAL(wp), DIMENSION(2,2) :: jacob

  !-----------------------------------------------------------------------

    jacob(1,1) = -(1._wp - eta) * x(1) + (1._wp - eta) * x(2)  &
      &        +  (1._wp + eta) * x(3) - (1._wp + eta) * x(4)
    jacob(1,2) = -(1._wp - eta) * y(1) + (1._wp - eta) * y(2)  &
      &        +  (1._wp + eta) * y(3) - (1._wp + eta) * y(4)
    jacob(2,1) = -(1._wp - zeta)* x(1) - (1._wp + zeta)* x(2)  &
      &        +  (1._wp + zeta)* x(3) + (1._wp - zeta)* x(4)
    jacob(2,2) = -(1._wp - zeta)* y(1) - (1._wp + zeta)* y(2)  &
      &        +  (1._wp + zeta)* y(3) + (1._wp - zeta)* y(4)

    det_jac = 0.0625_wp * (jacob(1,1)*jacob(2,2) - jacob(1,2)*jacob(2,1))

  END FUNCTION jac


END MODULE mo_advection_utils

