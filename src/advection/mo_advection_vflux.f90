!>
!! Computation of vertical tracer flux
!!
!! Vertical upwind fluxes are calculated at triangle centers on half-levels.
!! Possible options for vertical flux calculation include
!! - first order Godunov method (UP1)
!! - second order MUSCL method with or without CFL restriction
!! - third order PPM method with or without CFL restriction
!!
!! Semi-monotone and monotone limiters are available for MUSCL and PPM
!!
!! These routines compute only the correct half level value of
!! 'c*weta'. The vertical divergence is computed in step_advection.
!!
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Jochen Foerstner, DWD (2008-05-15)
!! Modification by Daniel Reinert, DWD (2009-08-06)
!! - code restructured
!! Modification by Daniel Reinert, DWD (2009-08-12)
!! - included piecewise parabolic method (PPM)
!! Modification by Daniel Reinert, DWD (2009-08-12)
!! - recoding of MUSCL in order to account for time and space
!!   dependent surface pressure and layer thickness
!! Modification by Daniel Reinert, DWD (2010-01-22)
!! - modified MUSCL scheme which handles CFL>1 (see Lin and Rood (1996))
!!   added optional semi-monotone limiter.
!! Modification by Daniel Reinert, DWD (2010-01-09)
!! - transferred vertical flux calculation (UP, MUSCL, PPM) to this
!!   new module
!! Modification by Daniel Reinert, DWD (2010-04-23)
!! - moved slope limiter to new module mo_advection_limiter
!!
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
MODULE mo_advection_vflux

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, SUCCESS, min_rlcell_int
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_model_domain,        ONLY: t_patch
  USE mo_run_nml,             ONLY: nproma, ntracer, iequations, msg_level,   &
    &                               lvert_nest
  USE mo_advection_utils,     ONLY: laxfr_upflux_v, laxfr_upflux, coeff_grid, &
    &                               iup_v, imuscl_v, imuscl_vcfl, ippm_v,     &
    &                               ippm_vcfl, islopel_vsm, islopel_vm,       &
    &                               ifluxl_vpd, ino_flx, izero_grad,          &
    &                               iparent_flx, lcompute, lcleanup
  USE mo_advection_limiter,   ONLY: v_muscl_slimiter_mo, v_muscl_slimiter_sm, &
   &                                v_ppm_slimiter_mo, v_ppm_slimiter_sm,     &
   &                                vflx_limiter_pd, vflx_limiter_pd_ha
  USE mo_loopindices,         ONLY: get_indices_c

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: vert_upwind_flux, upwind_vflux_up, upwind_vflux_muscl,  &
    &       upwind_vflux_ppm, upwind_vflux_ppm_cfl

  !-------------------------------------------------------------------------

CONTAINS


  !-------------------------------------------------------------------------
  !
  !

  !>
  !! Calculation of vertical upwind flux at triangle centers on half levels
  !!
  !! Calculation of vertical upwind flux at triangle centers on half levels
  !! using either
  !! - the first order Godunov method (UP1)
  !! - the second order MUSCL method
  !! - the third order PPM method
  !!
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-02-10)
  !! Modification by Daniel Reinert (2010-11-05)
  !! - tracer loop moved from step_advection to vert_upwind_flux
  !! Modification by Daniel Reinert (2010-11-08)
  !! - lcompute and lcleanup are precomputed in setup_transport.
  !!   This was necessary to allow for different flux-methods
  !!   for different tracers.
  !! Modification by Daniel Reinert, DWD (2011-02-15)
  !! - new field providing the upper margin tracer flux (required)
  !!
  !
  ! !LITERATURE
  ! MUSCL: Ahmad et al. (2006), Int. J. Num. Meth. Fluids, 50, 1247-1268
  !        Lin et al. (1994), MWR, 122, 1575-1593
  ! PPM  : Colella and Woodward (1984), JCP, 54, 174-201
  !        Carpenter et al. (1989), MWR, 118, 586-612
  !        Lin and Rood (1996), MWR, 124, 2046-2070 (see also for CFL-
  !                                                  independent versions)
  !
  SUBROUTINE vert_upwind_flux( p_patch, p_cc, p_mflx_contra_v, p_w_contra,    &
    &                      p_dtime, p_pres_ic, p_pres_mc, p_cellhgt_mc_now,   &
    &                      p_rcellhgt_mc_now, p_cellmass_now,                 &
    &                      p_ivadv_tracer, p_itype_vlimit, p_iubc_adv,        &
    &                      p_iadv_slev, p_upflux, opt_topflx_tra, opt_q_int,  &
    &                      opt_rho_ic, opt_rlstart, opt_rlend )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is 
      &  p_patch                             !< performed

    REAL(wp), INTENT(IN) ::  &      !< advected cell centered variable
      &  p_cc(:,:,:,:)              !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp), INTENT(INOUT) ::  &   !< contravariant vertical mass flux
      &  p_mflx_contra_v(:,:,:)     !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::  &      !< contravariant vertical velocity
      &  p_w_contra(:,:,:)          !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) :: p_dtime !< time step

    REAL(wp), INTENT(IN) ::  &      !< half level pressure at time n
      &  p_pres_ic(:,:,:)           !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::  &      !< full level pressure at time n
      &  p_pres_mc(:,:,:)           !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::  &      !< cell height defined at full levels for
      &  p_cellhgt_mc_now(:,:,:)    !< time step n (either \Delta p or \Delta z)
                                    !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::  &      !< reciprocal of cell height at
      &  p_rcellhgt_mc_now(:,:,:)   !< full levels at time step n

    REAL(wp), TARGET, INTENT(IN)::& !< NH: density weighted cell height at full levels
      &  p_cellmass_now(:,:,:)      !< at time step n [kg/m**2]
                                    !< HA: pressure thickness for full levels 
                                    !< at time step n [Pa]
                                    !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN) ::   &      !< parameter to select numerical
      &  p_ivadv_tracer(:)          !< scheme for vertical transport
                                    !< dim: (ntracer)

    INTEGER, INTENT(IN) ::   &      !< parameter to select the limiter
      &  p_itype_vlimit(:)          !< for vertical transport
                                    !< dim: (ntracer)

    INTEGER, INTENT(IN) ::   &      !< vertical start level for transport
      &  p_iadv_slev(:)             !< dim: (ntracer)

    INTEGER, INTENT(IN) ::   &      !< selects upper boundary condition
      &  p_iubc_adv

    REAL(wp), INTENT(OUT) :: &      !< variable in which the upwind flux is stored
      &  p_upflux(:,:,:,:)          !< dim: (nproma,nlevp1,nblks_c,ntracer)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< vertical tracer flux at upper boundary 
      &  opt_topflx_tra(:,:,:)          !< NH: [kg/m**2/s]
                                        !< HA: [Pa/s]
                                        !< dim: (nproma,nblks_c,ntracer)

    REAL(wp), INTENT(OUT), OPTIONAL :: & !< tracer value at upper boundary of child nest 
      &  opt_q_int(:,:,:)               !< NH: [kg/kg]
                                        !< HA: [kg/kg]
                                        !< dim: (nproma,nblks_c,ntracer)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< half level density at n+1/2 (only necessary
      &  opt_rho_ic(:,:,:)              !< for the NH-core, when muscl_cfl 
                                        !< is applied in vertical direction
                                        !< dim: (nproma,nlevp1,nblks_c)
    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    INTEGER :: jt                   !< tracer loop index
    INTEGER :: jc, jb               !< cell and block loop index
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart_c, i_rlend_c, i_nchdom
    !-----------------------------------------------------------------------
    !
    ! Loop over different tracers
    !
    ! Note (DR): Since we define p_ivadv_tracer as a 1D array of dimension
    ! ntracer, we can select different flux calculation methods for
    ! different tracers. We may also decide not to advect special
    ! tracers. Furthermore, options regarding the desired limiter can be passed 
    ! via the argument list. The same is true for precomputed lcompute/lcleanup
    ! values.
    IF (PRESENT(opt_topflx_tra)) THEN

      DO jt = 1, ntracer

        ! Select desired flux calculation method
        SELECT  CASE( p_ivadv_tracer(jt) )

        CASE( iup_v )
          ! CALL first order upwind
          CALL upwind_vflux_up( p_patch, p_cc(:,:,:,jt), p_iubc_adv,   &! in
            &                   p_mflx_contra_v, p_upflux(:,:,:,jt),   &! in,out
            &                   opt_topflx_tra=opt_topflx_tra(:,:,jt), &! in
            &                   opt_slev=p_iadv_slev(jt),              &! in
            &                   opt_rlstart=opt_rlstart,               &! in
            &                   opt_rlend=opt_rlend                    )! in


        CASE( imuscl_v, imuscl_vcfl )
          ! CALL second order MUSCL
          CALL upwind_vflux_muscl( p_patch, p_cc(:,:,:,jt), p_iubc_adv,      &! in
            &                      p_mflx_contra_v, p_w_contra, p_dtime,     &! in
            &                      lcompute%muscl_v(jt),                     &! in
            &                      lcleanup%muscl_v(jt), p_ivadv_tracer(jt), &! in
            &                      p_itype_vlimit(jt), p_pres_ic, p_pres_mc, &! in
            &                      p_cellhgt_mc_now, p_rcellhgt_mc_now,      &! in
            &                      p_upflux(:,:,:,jt), opt_rho_ic=opt_rho_ic,&! out,in
            &                      opt_topflx_tra=opt_topflx_tra(:,:,jt),    &! in
            &                      opt_slev=p_iadv_slev(jt),                 &! in
            &                      opt_rlstart=opt_rlstart,                  &! in
            &                      opt_rlend=opt_rlend                       )! in
        CASE( ippm_v )
          ! CALL third order PPM
          CALL upwind_vflux_ppm( p_patch, p_cc(:,:,:,jt), p_iubc_adv,  &! in
            &                  p_mflx_contra_v, p_w_contra, p_dtime,   &! in
            &                  p_itype_vlimit(jt), p_cellhgt_mc_now,   &! in
            &                  p_rcellhgt_mc_now, p_upflux(:,:,:,jt),  &! in,out
            &                  opt_topflx_tra=opt_topflx_tra(:,:,jt),  &! in
            &                  opt_slev=p_iadv_slev(jt),               &! in
            &                  opt_rlstart=opt_rlstart,                &! in
            &                  opt_rlend=opt_rlend                     )! in

        CASE( ippm_vcfl )
          ! CALL third order PPM which handles long time steps (i.e. CFL>1)
          CALL upwind_vflux_ppm_cfl( p_patch, p_cc(:,:,:,jt), p_iubc_adv,    &! in
            &                  p_mflx_contra_v, p_dtime, lcompute%ppm_v(jt), &! in
            &                  lcleanup%ppm_v(jt), p_itype_vlimit(jt),       &! in
            &                  p_cellhgt_mc_now, p_cellmass_now,             &! in
            &                  p_upflux(:,:,:,jt),                           &! out
            &                  opt_topflx_tra=opt_topflx_tra(:,:,jt),        &! in
            &                  opt_slev=p_iadv_slev(jt),                     &! in
            &                  opt_rlstart=opt_rlstart,                      &! in
            &                  opt_rlend=opt_rlend                           )! in
        END SELECT
      END DO  ! Tracer loop


      !
      ! get face value "q_int" at vertical boundary of child nest
      !
      ! determine if upper boundary values are needed
      IF (lvert_nest .AND. (p_patch%nshift_child > 0)) THEN 

        ! refinement control start/end level for cells
        i_rlstart_c = 1
        i_rlend_c   = min_rlcell_int

        ! number of child domains
        i_nchdom = MAX(1,p_patch%n_childdom)

        i_startblk = p_patch%cells%start_blk(i_rlstart_c,1)
        i_endblk   = p_patch%cells%end_blk(i_rlend_c,i_nchdom)

        DO jb = i_startblk, i_endblk
          CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
            &                 i_startidx, i_endidx, i_rlstart_c, i_rlend_c )

          DO jt = 1, ntracer
            DO jc = i_startidx, i_endidx
              opt_q_int(jc,jb,jt) = p_upflux(jc,p_patch%nshift_child,jb,jt) &
                &        /(p_mflx_contra_v(jc,p_patch%nshift_child,jb) + dbl_eps)
            ENDDO
          ENDDO
        ENDDO

      ENDIF


    ELSE ! opt_topflx_tra not present (i.e. for hydrostatic model)

      DO jt = 1, ntracer

        ! Select desired flux calculation method
        SELECT  CASE( p_ivadv_tracer(jt) )

        CASE( iup_v )
          ! CALL first order upwind
          CALL upwind_vflux_up( p_patch, p_cc(:,:,:,jt), p_iubc_adv,   &! in
            &                   p_mflx_contra_v, p_upflux(:,:,:,jt),   &! in,out
            &                   opt_slev=p_iadv_slev(jt),              &! in
            &                   opt_rlstart=opt_rlstart,               &! in
            &                   opt_rlend=opt_rlend                    )! in

        CASE( imuscl_v, imuscl_vcfl )
          ! CALL second order MUSCL
          CALL upwind_vflux_muscl( p_patch, p_cc(:,:,:,jt), p_iubc_adv,      &! in
            &                      p_mflx_contra_v, p_w_contra, p_dtime,     &! in
            &                      lcompute%muscl_v(jt),                     &! in
            &                      lcleanup%muscl_v(jt), p_ivadv_tracer(jt), &! in
            &                      p_itype_vlimit(jt), p_pres_ic, p_pres_mc, &! in
            &                      p_cellhgt_mc_now, p_rcellhgt_mc_now,      &! in
            &                      p_upflux(:,:,:,jt), opt_rho_ic=opt_rho_ic,&! out,in
            &                      opt_slev=p_iadv_slev(jt),                 &! in
            &                      opt_rlstart=opt_rlstart,                  &! in
            &                      opt_rlend=opt_rlend                       )! in

        CASE( ippm_v )
          ! CALL third order PPM
          CALL upwind_vflux_ppm( p_patch, p_cc(:,:,:,jt), p_iubc_adv,  &! in
            &                  p_mflx_contra_v, p_w_contra, p_dtime,   &! in
            &                  p_itype_vlimit(jt), p_cellhgt_mc_now,   &! in
            &                  p_rcellhgt_mc_now, p_upflux(:,:,:,jt),  &! in,out
            &                  opt_slev=p_iadv_slev(jt),               &! in
            &                  opt_rlstart=opt_rlstart,                &! in
            &                  opt_rlend=opt_rlend                     )! in

        CASE( ippm_vcfl )
          ! CALL third order PPM which handles long time steps (i.e. CFL>1)
          CALL upwind_vflux_ppm_cfl( p_patch, p_cc(:,:,:,jt), p_iubc_adv,    &! in
            &                  p_mflx_contra_v, p_dtime, lcompute%ppm_v(jt), &! in
            &                  lcleanup%ppm_v(jt), p_itype_vlimit(jt),       &! in
            &                  p_cellhgt_mc_now, p_cellmass_now,             &! in
            &                  p_upflux(:,:,:,jt),                           &! out
            &                  opt_slev=p_iadv_slev(jt),                     &! in
            &                  opt_rlstart=opt_rlstart,                      &! in
            &                  opt_rlend=opt_rlend                           )! in

        END SELECT
      END DO  ! Tracer loop

    END IF ! PRESENT(opt_topflx_tra)

  END SUBROUTINE vert_upwind_flux




  !-------------------------------------------------------------------------
  !>
  !! The first order Godunov method
  !!
  !! Calculation of time averaged vertical tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !! Initial revision by Jochen Foerstner, DWD (2008-05-15)
  !! Modification by Daniel Reinert, DWD (2010-02-09)
  !! - transferred to separate subroutine
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized to height based vertical coordinate systems
  !!
  SUBROUTINE upwind_vflux_up( p_patch, p_cc, p_iubc_adv, p_mflx_contra_v, &
    &                         p_upflux, opt_topflx_tra, opt_slev,         &
    &                         opt_rlstart, opt_rlend )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_vflux: upwind_vflux_up'

    TYPE(t_patch), TARGET, INTENT(IN) ::  & !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::   &   !< advected cell centered variable
      &  p_cc(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN)  ::   &   !< selects upper boundary condition
      &  p_iubc_adv

    REAL(wp), INTENT(IN) ::   &   !< contravariant vertical mass flux
      &  p_mflx_contra_v(:,:,:)   !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(OUT) ::  &   !< vertical tracer flux at half levels
      &  p_upflux(:,:,:)          !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< vertical tracer flux at upper boundary 
      &  opt_topflx_tra(:,:)            !< dim: (nproma,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    REAL(wp) ::  &                             !< necessary, to make this routine
     &  zparent_topflx(nproma,p_patch%nblks_c) !< compatible to the hydrost. core 
                                       
    INTEGER  :: slev                   !< vertical start level
    INTEGER  :: nlev, nlevp1           !< number of full and half levels
    INTEGER  :: jc, jk, jb             !< index of cell, vertical level and block
    INTEGER  :: jkm1                   !< jk - 1
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    !-------------------------------------------------------------------------

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_topflx_tra) ) THEN
      zparent_topflx(:,:) = opt_topflx_tra(:,:)
    ELSE
      zparent_topflx(:,:) = 0._wp
    ENDIF

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

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    !
    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jkm1,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )


      DO jk = slev+1, nlev

        ! index of top half level
        jkm1 = jk - 1

        DO jc = i_startidx, i_endidx

          ! calculate vertical tracer flux
          p_upflux(jc,jk,jb) =                                  &
            &  laxfr_upflux_v( p_mflx_contra_v(jc,jk,jb),       &
            &                p_cc(jc,jkm1,jb), p_cc(jc,jk,jb),  &
            &                coeff_grid )

        END DO ! end loop over cells

      ENDDO ! end loop over vertical levels


      !
      ! set upper and lower boundary condition
      !
      CALL set_bc_vadv(p_upflux(:,slev+1,jb),            &! in
        &              p_mflx_contra_v(:,slev+1,jb),     &! in
        &              p_mflx_contra_v(:,slev  ,jb),     &! in
        &              p_iubc_adv, i_startidx, i_endidx, &! in
        &              zparent_topflx(:,jb),             &! in
        &              p_upflux(:,slev,jb),              &! out
        &              p_upflux(:,nlevp1,jb)             )! out

    ENDDO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE upwind_vflux_up




  !-------------------------------------------------------------------------
  !>
  !! The second order MUSCL scheme
  !!
  !! Calculation of time averaged vertical tracer fluxes using the second
  !! order MUSCL scheme.
  !!
  !! @par Revision History
  !! Initial revision by Jochen Foerstner, DWD (2008-05-15)
  !! Modification by Daniel Reinert, DWD (2009-08-06)
  !! - code restructured
  !! Modification by Daniel Reinert, DWD (2009-08-12)
  !! - recoding of MUSCL in order to account for time and space
  !!   dependent surface pressure and thus layer thickness
  !! Modification by Daniel Reinert, DWD (2009-08-12)
  !! - modified MUSCL scheme which handles CFL>1 (see Lin and Rood (1996))
  !!   and an optional positive definite limiter.
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized to height based vertical coordinate systems. Included
  !!   parameter coeff_grid, in order to apply the same code to either a
  !!   pressure based or height based vertical coordinate system.
  !!
  ! !LITERATURE
  ! MUSCL: Ahmad et al. (2006), Int. J. Num. Meth. Fluids, 50, 1247-1268
  !        Lin et al. (1994), MWR, 122, 1575-1593
  !
  SUBROUTINE upwind_vflux_muscl( p_patch, p_cc, p_iubc_adv, p_mflx_contra_v,  &
    &                        p_w_contra, p_dtime, ld_compute, ld_cleanup,     &
    &                        p_ivadv_tracer, p_itype_vlimit, p_pres_ic,       &
    &                        p_pres_mc, p_cellhgt_mc_now, p_rcellhgt_mc_now,  &
    &                        p_upflux, opt_rho_ic, opt_topflx_tra, opt_slev,  &
    &                        opt_rlstart, opt_rlend )


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_vflux: upwind_vflux_muscl'

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::  &    !< advected cell centered variable
      &  p_cc(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN)  ::   &   !< selects upper boundary condition
      &  p_iubc_adv

    REAL(wp), INTENT(IN) ::  &    !< contravariant vertical mass flux
      &  p_mflx_contra_v(:,:,:)   !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< contravariant vertical velocity
      &  p_w_contra(:,:,:)        !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) :: p_dtime  !< time step

    REAL(wp), INTENT(IN) ::  &    !< half level pressure at time n
      &  p_pres_ic(:,:,:)         !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< full level pressure at time n
      &  p_pres_mc(:,:,:)         !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< layer thickness at cell center at time n
      &  p_cellhgt_mc_now(:,:,:)  !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< reciprocal of layer thickness at cell center
      &  p_rcellhgt_mc_now(:,:,:) !< at time n; dim: (nproma,nlev,nblks_c)

    LOGICAL, INTENT(IN) ::   &    !< flag, if .TRUE. compute geometrical terms
      &  ld_compute

    LOGICAL, INTENT(IN) ::   &    !< flag, if .TRUE. clean up geometrical terms
      &  ld_cleanup

    INTEGER, INTENT(IN) ::   &    !< parameter to select numerical
      &  p_ivadv_tracer           !< scheme for vertical transport

    INTEGER, INTENT(IN) :: p_itype_vlimit    !< parameter to select the limiter
                                             !< for vertical transport

    REAL(wp), INTENT(OUT) ::  &   !< variable in which the upwind flux is stored
      &  p_upflux(:,:,:)          !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< half level density at n+1/2 (only necessary
      &  opt_rho_ic(:,:,:)              !< for the NH-core, when muscl_cfl or ppm_cfl
                                        !< is applied in vertical direction
                                        !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< vertical tracer flux at upper boundary 
      &  opt_topflx_tra(:,:)            !< dim: (nproma,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    REAL(wp), ALLOCATABLE, SAVE :: &
      &  z_delp_c1(:,:,:),  &        !< pressure difference between half level and
      &  z_delp_c2(:,:,:)            !< corresponding upper (z_delp_c1) and lower
                                     !< (z_delp_c2) full level.
                                     !< (nproma,nlevp1,p_patch%nblks_c)

    REAL(wp), ALLOCATABLE, SAVE :: &
      &  z_alpha_1(:,:,:)            !< weight for linear interpolation from full
                                     !< level to half level (MUSCL)
                                     !< (nproma,nlevp1,p_patch%nblks_c)

    REAL(wp), ALLOCATABLE, SAVE :: &
      &  z_gmoment(:,:,:)            !< geometrical moment, which accounts for the fact,
                                     !< that the mass point is not located at the center
                                     !< of mass of the vertical grid cell for the
                                     !< hydrostatic model version
                                     !< (nproma,nlev,p_patch%nblks_c)

    REAL(wp), ALLOCATABLE, SAVE :: & !< CFL number (weta>0, w<0)
      &  z_cfl_m0(:,:,:)             !< (nproma,nlevp1,p_patch%nblks_c)

    REAL(wp), ALLOCATABLE, SAVE :: & !< CFL number (weta<0, w>0)
      &  z_cfl_p0(:,:,:)             !< (nproma,nlevp1,p_patch%nblks_c)

    REAL(wp) :: &                    !< face values of transported field
      &  z_face(nproma,p_patch%nlevp1,p_patch%nblks_c)
   
    REAL(wp) :: &                    !< gradient at cell center
      &  z_grad(nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp) :: &                    !< linear extrapolation value 1
      &  z_lext_1(nproma,p_patch%nlevp1)

    REAL(wp) :: &                    !< linear extrapolation value 2
      &  z_lext_2(nproma,p_patch%nlevp1)

    REAL(wp) :: &                    !< integer fluxes for weta>0, w<0
      &  z_iflx_m(nproma,p_patch%nlevp1)

    REAL(wp) :: &                    !< integer fluxes for weta<0, w>0
      &  z_iflx_p(nproma,p_patch%nlevp1)

    REAL(wp) :: &                    !< copy of z_cfl_m0 for each tracer, to avoid loss
      &  z_cfl_m(nproma,p_patch%nlevp1,p_patch%nblks_c) !< of original values
                                                        
    REAL(wp) :: &                    !< dito for z_cfl_p0
      &  z_cfl_p(nproma,p_patch%nlevp1,p_patch%nblks_c)

    REAL(wp) ::  &                             !< necessary, to make this routine
     &  zparent_topflx(nproma,p_patch%nblks_c) !< compatible to the hydrost. core 

    INTEGER  :: slev, slevp1             !< vertical start level and start level +1
    INTEGER  :: nlev, nlevp1             !< number of full and half levels
    INTEGER  :: jkm1, ikp1, ikp1_ic      !< vertical level minus and plus one
    INTEGER  :: cfl_counter              !< checks, wheter we encountered CFL>1
    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: jkkp, jkkm
    INTEGER  :: ist                      !< status variable
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER, DIMENSION(nproma)::   &     !< similar to jkm and jk, but number of
      &   jkm_int, jk_int                !< cells to skip for CFL>1 have been
                                         !< substracted/added

!DR    REAL(wp) :: cfl_m_max, cfl_p_max !< maximum Courant number from both sides

    !-------------------------------------------------------------------------

    !
    ! advection is done with an upwind scheme and a piecwise linear
    ! approx. of the cell centered values at the edges
    ! is used.

    ! 3 options:  standard without limiter
    !             standard with positive definite or monotone limiter
    !             special version with limiter which can handle CFL >1
    !

    ! check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev   = opt_slev
      slevp1 = opt_slev + 1
    ELSE
      slev   = 1
      slevp1 = 2
    END IF

    IF ( PRESENT(opt_topflx_tra) ) THEN
      zparent_topflx(:,:) = opt_topflx_tra(:,:)
    ELSE
      zparent_topflx(:,:) = 0._wp
    ENDIF

    IF (.NOT. PRESENT(opt_rho_ic) .AND. iequations ==3 .AND. &
      &   p_ivadv_tracer == imuscl_vcfl ) THEN
        CALL finish ( TRIM(routine),                         &
          &    'optional argument opt_rho_ic is missing' )
    ENDIF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c-1
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell_int
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)


    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1


    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

    IF ( ld_compute ) THEN
      ! allocate temporary arrays for layer thickness, weights and
      ! Courant numbers
      ALLOCATE( z_delp_c1(nproma,nlevp1,p_patch%nblks_c),    &
        &       z_delp_c2(nproma,nlevp1,p_patch%nblks_c),    &
        &       z_alpha_1(nproma,nlevp1,p_patch%nblks_c),    &
        &       z_gmoment(nproma,nlev,p_patch%nblks_c),      &
        &       z_cfl_m0(nproma,nlevp1,p_patch%nblks_c),     &
        &       z_cfl_p0(nproma,nlevp1,p_patch%nblks_c),     &
        &       STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                         &
          &  'allocation for z_delp_c1, z_delp_c2, '     //  &
          &  'z_alpha_1, z_gmoment, z_cfl_m0, z_cfl_p0 failed' )
      ENDIF
    END IF

    !
    ! 0. precalculate time dependent vertical distances, weights and
    !    geometrical moment in terms of pressure (at time step n) for
    !    vertical differencing
    !
    !    Calculate Courant number at cell interfaces for the cases
    !    weta >0 and weta <0
    !
!$OMP PARALLEL
    IF (ld_compute) THEN
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,jkm1)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend )


        ! boundary conditions for uppermost and lowermost half level
        z_delp_c1(i_startidx:i_endidx,1,jb) = 0._wp
        z_delp_c2(i_startidx:i_endidx,1,jb) = p_pres_ic(i_startidx:i_endidx,1,jb) &
          &                                 - p_pres_mc(i_startidx:i_endidx,1,jb)
        z_delp_c1(i_startidx:i_endidx,nlevp1,jb) = p_pres_ic(i_startidx:i_endidx,nlevp1,jb) &
          &                                      - p_pres_mc(i_startidx:i_endidx,nlev,jb)
        z_delp_c2(i_startidx:i_endidx,nlevp1,jb) = 0._wp

        z_alpha_1(i_startidx:i_endidx,1,jb)      = 0._wp
        z_alpha_1(i_startidx:i_endidx,nlevp1,jb) = 1._wp

        ! geometrical moment
        z_gmoment(i_startidx:i_endidx,1,jb) =                                   &
          &                   0.5_wp*p_cellhgt_mc_now(i_startidx:i_endidx,1,jb) &
          &                   + coeff_grid * z_delp_c2(i_startidx:i_endidx,1,jb)

        ! Courant number at top and bottom
        z_cfl_m0(i_startidx:i_endidx,1,jb)      = 0._wp
        z_cfl_p0(i_startidx:i_endidx,1,jb)      = 0._wp
        z_cfl_m0(i_startidx:i_endidx,nlevp1,jb) = 0._wp
        z_cfl_p0(i_startidx:i_endidx,nlevp1,jb) = 0._wp

        DO jk = 2, nlev

          ! index of top half level
          jkm1 = jk - 1

          DO jc = i_startidx, i_endidx

            ! pressure difference between half level and corresponding upper
            ! and lower full level
            ! z_delp_c1: half level - upper full level (>0)
            ! z_delp_c2: half level - lower full level (<0)
            z_delp_c1(jc,jk,jb) = p_pres_ic(jc,jk,jb) - p_pres_mc(jc,jkm1,jb)
            z_delp_c2(jc,jk,jb) = p_pres_ic(jc,jk,jb) - p_pres_mc(jc,  jk,jb)


            ! weight for linear interpolation of tracer values
            ! the second weight is not stored. It's simply (1-z_alpha_1)
            z_alpha_1(jc,jk,jb) = z_delp_c1(jc,jk,jb)                           &
              &                 / ( p_pres_mc(jc,jk,jb) - p_pres_mc(jc,jkm1,jb) )


            ! geometrical moment which forces the piecewise linear
            ! reconstruction to be conservative
            z_gmoment(jc,jk,jb) = 0.5_wp*p_cellhgt_mc_now(jc,jk,jb)       &
              &                 + coeff_grid * z_delp_c2(jc,jk,jb)


            ! Calculate local Courant number at half levels
            ! z_cfl_m0 for the case weta >0 (w <0)
            ! z_cfl_p0 for the case weta <0 (w >0)
            z_cfl_m0(jc,jk,jb) = ABS(p_w_contra(jc,jk,jb))               &
              &                * p_dtime * p_rcellhgt_mc_now(jc,jkm1,jb)

            z_cfl_p0(jc,jk,jb) = ABS(p_w_contra(jc,jk,jb))               &
              &                * p_dtime * p_rcellhgt_mc_now(jc,jk,jb)

          END DO
        END DO
      ENDDO
!$OMP END DO

    ENDIF


    !
    ! 1. reconstruction of (unlimited) cell based vertical gradient
    !    and specification of boundary values for z_face and z_grad
    !
    !    In addition, copy z_cfl_[m,p]0 to work-array z_cfl_[m,p]
    !
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,ikp1_ic,ikp1)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

      ! set value of tracer at top vertical face
      z_face(i_startidx:i_endidx,slev,jb) = p_cc(i_startidx:i_endidx,slev,jb)

      ! calculate values of tracer at second vertical face
      z_face(i_startidx:i_endidx,slevp1,jb) = z_alpha_1(i_startidx:i_endidx,slevp1,jb) &
        &                    * p_cc(i_startidx:i_endidx,slevp1,jb)                     &
        &                    + ( 1._wp - z_alpha_1(i_startidx:i_endidx,slevp1,jb) )    &
        &                    * p_cc(i_startidx:i_endidx,slev,jb)

      ! vertical gradient at top
      z_grad(i_startidx:i_endidx,slev,jb) = coeff_grid                                 &
        &                      * ( z_face(i_startidx:i_endidx,slevp1,jb)               &
        &                      - z_face(i_startidx:i_endidx,slev,jb) )                 &
        &                      * p_rcellhgt_mc_now(i_startidx:i_endidx,slev,jb)

      ! local Courant number at half levels
      ! z_cfl_m for the case weta >0 (w <0)
      ! z_cfl_p for the case weta <0 (w >0)
      z_cfl_m(i_startidx:i_endidx,slev:nlevp1,jb) =                              &
        &                          z_cfl_m0(i_startidx:i_endidx,slev:nlevp1,jb)
      z_cfl_p(i_startidx:i_endidx,slev:nlevp1,jb) =                              &
        &                          z_cfl_p0(i_startidx:i_endidx,slev:nlevp1,jb)

      DO jk = slevp1, nlev

        ! index of bottom half level
        ikp1_ic = jk + 1
        ikp1 = MIN( ikp1_ic, nlev )

        DO jc = i_startidx, i_endidx

          ! calculate tracer values at vertical face
          z_face(jc,ikp1_ic,jb) =                                         &
            &    z_alpha_1(jc,ikp1_ic,jb) * p_cc(jc,ikp1,jb)              &
            &   + ( 1._wp - z_alpha_1(jc,ikp1_ic,jb) ) * p_cc(jc,jk,jb)

          ! vertical gradient
          z_grad(jc,jk,jb) = coeff_grid                                   &
            &              * ( z_face(jc,ikp1_ic,jb) - z_face(jc,jk,jb) ) &
            &              * p_rcellhgt_mc_now(jc,jk,jb)

        END DO ! end loop over cells

      ENDDO ! end loop over vertical levels

    ENDDO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL


!DR  cfl_p_max=MAXVAL(z_cfl_p(:,34:54,:))
!DR  cfl_m_max=MAXVAL(z_cfl_m(:,34:54,:))
!DR  WRITE(message_text,'(a,e16.8)') 'CFL_P MAX', cfl_p_max
!DR  CALL message('vertical_upwind_flux',message_text)
!DR  WRITE(message_text,'(a,e16.8)') 'CFL_M MAX', cfl_m_max
!DR  CALL message('vertical_upwind_flux',message_text)




    ! 2. If desired, slopes are limited using either a monotonous or
    !    a semi-monotonous version of the Barth-Jespersen slope limiter.
    !
    IF (p_itype_vlimit == islopel_vsm) THEN
      ! semi-monotonic (sm) slope limiter
      CALL v_muscl_slimiter_sm( p_patch, p_cc, z_gmoment, z_delp_c1,      &! in
        &                       z_delp_c2, z_grad, opt_rlstart=i_rlstart, &! in,inout,in
        &                       opt_rlend=i_rlend, opt_slev=slev          )! in

    ELSE IF (p_itype_vlimit == islopel_vm) THEN
      ! monotonous (mo) slope limiter
      CALL v_muscl_slimiter_mo( p_patch, p_cc, z_gmoment, z_delp_c1,      &! in
        &                       z_delp_c2, z_grad, opt_rlstart=i_rlstart, &! in,inout,in
        &                       opt_rlend=i_rlend, opt_slev=slev          )! in
    ENDIF



    SELECT  CASE( p_ivadv_tracer )
    CASE( imuscl_vcfl )           ! works for CFL >1

    !
    ! 3. calculation of upwind fluxes. IF CFL > 1, the fluxes are computed as
    !    the sum of integer-fluxes and one fractional flux. IF CFL <1 the
    !    fluxes are only comprised of the fractional flux. The fractional
    !    flux is calculated assuming a piecewise linear approx. for the
    !    subgrid distribution.
    !

    ! Comment: In principle it is possible to generalize this code snippet, too.
    ! Unfortunately this would lead to some computational overhead which can be avoided
    ! (at least partly), by writing distinct code snippet for the hydrostatic and
    ! non-hydrostatic dynamical core. Note that this is only necessary for the Courant-
    ! number independent version.

      IF (iequations /= 3) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_iflx_m,z_iflx_p,z_lext_1, &
!$OMP            z_lext_2,jkm1,cfl_counter,jkm_int,jk_int,jkkm,jkkp)
        DO jb = i_startblk, i_endblk

          CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
            &                 i_startidx, i_endidx, i_rlstart, i_rlend )


          z_iflx_m(i_startidx:i_endidx,slev:nlevp1) = 0._wp
          z_iflx_p(i_startidx:i_endidx,slev:nlevp1) = 0._wp
          z_lext_1(i_startidx:i_endidx,slev)        = p_cc(i_startidx:i_endidx,slev,jb)
          z_lext_2(i_startidx:i_endidx,slev)        = p_cc(i_startidx:i_endidx,slev,jb)


          DO jk = slevp1, nlev

            jkm1 = jk - 1
            jkkm = jk
            jkkp = jk

            jk_int(i_startidx:i_endidx)  = jk     ! Initialization
            jkm_int(i_startidx:i_endidx) = jkm1

            cfl_counter = 0


            ! Loop for calculation of Integer fluxes
            DO WHILE ( (MAXVAL(z_cfl_m(i_startidx:i_endidx,jk,jb)) > 1._wp  &
                 & .OR. MAXVAL(z_cfl_p(i_startidx:i_endidx,jk,jb)) > 1._wp) &
                 & .AND. jkkm > 2 .AND. jkkp < nlev )


              cfl_counter = cfl_counter + 1  ! Checks, whether we stepped into this loop, or not
              jkkm = jkkm - 1
              jkkp = jkkp + 1

              DO jc = i_startidx, i_endidx

                ! Integer fluxes of cell above (case of weta > 0; w < 0)
                IF ( z_cfl_m(jc,jk,jb) > 1._wp ) THEN


                  ! Integer flux (division by p_dtime is done at the end)
                  z_iflx_m(jc,jk) = z_iflx_m(jc,jk) + p_cc(jc,jkm_int(jc),jb) &
                    &              * p_cellhgt_mc_now(jc,jkm_int(jc),jb)


                  ! Courant number - 1
                  z_cfl_m(jc,jk,jb) = ( z_cfl_m(jc,jk,jb) - 1._wp )           &
                    &               * p_cellhgt_mc_now(jc,jkm_int(jc),jb)     &
                    &               * p_rcellhgt_mc_now(jc,jkm_int(jc)-1,jb)


                  ! Account for number of cells which we need to step upwind
                  ! for calculating the fractional flux
                  jkm_int(jc) = jkm_int(jc) - 1
                ENDIF

                ! Integer fluxes of cell below (case of weta < 0; w > 0)
                IF ( z_cfl_p(jc,jk,jb) > 1._wp ) THEN

                  ! Integer flux (division by p_dtime is done at the end)
                  z_iflx_p(jc,jk) = z_iflx_p(jc,jk) - p_cc(jc,jk_int(jc),jb)  &
                    &              * p_cellhgt_mc_now(jc,jk_int(jc),jb)


                  ! Courant number - 1
                  z_cfl_p(jc,jk,jb) = ( z_cfl_p(jc,jk,jb) - 1._wp )           &
                    &                * p_cellhgt_mc_now(jc,jk_int(jc),jb)     &
                    &                * p_rcellhgt_mc_now(jc,jk_int(jc)+1,jb)

                  ! Account for number of cells which we need to step upwind
                  ! for calculating the fractional flux
                  jk_int(jc) = jk_int(jc) + 1
                ENDIF

              ENDDO  ! end loop over cells

            ENDDO  ! end while loop



            ! if CFL < 1, everywhere in (1:nlen,jk,jb) then use faster standard
            ! calculation of (fractional) upwind fluxes
            IF (cfl_counter == 0) THEN

              DO jc = i_startidx, i_endidx

                ! linear extrapolated values
                ! first (of cell above) (case of weta > 0; w < 0; physical downwelling)
                z_lext_1(jc,jk) = p_cc(jc,jkm1,jb) + 0.5_wp                         &
                  &             * z_grad(jc,jkm1,jb) * p_cellhgt_mc_now(jc,jkm1,jb) &
                  &             * (1._wp - MIN(1._wp,z_cfl_m(jc,jk,jb)))

                ! second (of cell below) (case of weta < 0; w > 0; physical upwelling)
                z_lext_2(jc,jk) = p_cc(jc,jk,jb) - 0.5_wp                         &
                  &             * z_grad(jc,jk,jb) *  p_cellhgt_mc_now(jc,jk,jb)  &
                  &             * (1._wp - MIN(1._wp,z_cfl_p(jc,jk,jb)))

                !
                ! calculate vertical tracer flux
                !
                p_upflux(jc,jk,jb) =                              &
                  &  laxfr_upflux( p_mflx_contra_v(jc,jk,jb),     &
                  &              z_lext_1(jc,jk), z_lext_2(jc,jk))

              END DO ! end loop over cells

            ELSE

              ! if CFL > 1 anywhere in (i_startidx:i_endidx,jk,jb), then the following 
              ! somewhat more expensive calculation of the sum of integer and fractional 
              ! upwind fluxes is needed.
              DO jc = i_startidx, i_endidx
                IF (p_mflx_contra_v(jc,jk,jb) > 0._wp ) THEN

                  ! Sum of Integer fluxes and fractional flux
                  ! if weta(jc,jk,jb) > 0._wp; w < 0
                  p_upflux(jc,jk,jb)= ( z_iflx_m(jc,jk)                           &
                    &             + p_cellhgt_mc_now(jc,jkm_int(jc),jb)           &
                    &             * z_cfl_m(jc,jk,jb) * ( p_cc(jc,jkm_int(jc),jb) &
                    &             + 0.5_wp * z_grad(jc,jkm_int(jc),jb)            &
                    &             * p_cellhgt_mc_now(jc,jkm_int(jc),jb)           &
                    &             * (1._wp - MIN(1._wp,z_cfl_m(jc,jk,jb)))) )     &
                    &             / p_dtime

                ELSE

                  ! Sum of Integer fluxes and fractional flux
                  ! if weta(jc,jk,jb) <= 0._wp; w >= 0
                  p_upflux(jc,jk,jb) = ( z_iflx_p(jc,jk)                          &
                    &             - p_cellhgt_mc_now(jc,jk_int(jc),jb)            &
                    &             * z_cfl_p(jc,jk,jb) * ( p_cc(jc,jk_int(jc),jb)  &
                    &             - 0.5_wp * z_grad(jc,jk_int(jc),jb)             &
                    &             * p_cellhgt_mc_now(jc,jk_int(jc),jb)            &
                    &             * (1._wp - MIN(1._wp,z_cfl_p(jc,jk,jb)))) )     &
                    &             / p_dtime

                ENDIF

              ENDDO ! end loop over cells

            ENDIF

          ENDDO ! end loop over vertical levels

          !
          ! set upper and lower boundary condition
          !
          CALL set_bc_vadv(p_upflux(:,slev+1,jb),            &! in
            &              p_mflx_contra_v(:,slev+1,jb),     &! in
            &              p_mflx_contra_v(:,slev  ,jb),     &! in
            &              p_iubc_adv, i_startidx, i_endidx, &! in
            &              zparent_topflx(:,jb),             &! in
            &              p_upflux(:,slev,jb),              &! out
            &              p_upflux(:,nlevp1,jb)             )! out

        ENDDO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

      ELSE  ! IF (iequations = 3), i.e. non-hydrostatic core

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_iflx_m,z_iflx_p,z_lext_1, &
!$OMP            z_lext_2,jkm1,cfl_counter,jkm_int,jk_int,jkkm,jkkp)
        DO jb = i_startblk, i_endblk

          CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
            &                 i_startidx, i_endidx, i_rlstart, i_rlend )


          z_iflx_m(i_startidx:i_endidx,slev:nlevp1) = 0._wp
          z_iflx_p(i_startidx:i_endidx,slev:nlevp1) = 0._wp
          z_lext_1(i_startidx:i_endidx,slev)        = p_cc(i_startidx:i_endidx,slev,jb)
          z_lext_2(i_startidx:i_endidx,slev)        = p_cc(i_startidx:i_endidx,slev,jb)


          DO jk = slevp1, nlev

            jkm1 = jk - 1
            jkkm = jk
            jkkp = jk

            jk_int(i_startidx:i_endidx)  = jk     ! Initialization
            jkm_int(i_startidx:i_endidx) = jkm1

            cfl_counter = 0


            ! Loop for calculation of Integer fluxes
            DO WHILE ( (MAXVAL(z_cfl_m(i_startidx:i_endidx,jk,jb)) > 1._wp  &
                 & .OR. MAXVAL(z_cfl_p(i_startidx:i_endidx,jk,jb)) > 1._wp) &
                 & .AND. jkkm > 2 .AND. jkkp < nlev )


              cfl_counter = cfl_counter + 1  ! Checks, whether we stepped into this loop, or not
              jkkm = jkkm - 1
              jkkp = jkkp + 1

              DO jc = i_startidx, i_endidx

                ! Integer fluxes of cell above (case of weta > 0; w < 0)
                IF ( z_cfl_m(jc,jk,jb) > 1._wp ) THEN


                  ! Integer flux (division by p_dtime is done at the end)
                  z_iflx_m(jc,jk) = z_iflx_m(jc,jk) - p_cc(jc,jkm_int(jc),jb) &
                    &              * p_cellhgt_mc_now(jc,jkm_int(jc),jb)      &
                    &              * opt_rho_ic(jc,jk,jb)


                  ! Courant number - 1
                  z_cfl_m(jc,jk,jb) = ( z_cfl_m(jc,jk,jb) - 1._wp )           &
                    &               * p_cellhgt_mc_now(jc,jkm_int(jc),jb)     &
                    &               * p_rcellhgt_mc_now(jc,jkm_int(jc)-1,jb)


                  ! Account for number of cells which we need to step upwind
                  ! for calculating the fractional flux
                  jkm_int(jc) = jkm_int(jc) - 1
                ENDIF

                ! Integer fluxes of cell below (case of weta < 0; w > 0)
                IF ( z_cfl_p(jc,jk,jb) > 1._wp ) THEN

                  ! Integer flux (division by p_dtime is done at the end)
                  z_iflx_p(jc,jk) = z_iflx_p(jc,jk) + p_cc(jc,jk_int(jc),jb)  &
                    &              * p_cellhgt_mc_now(jc,jk_int(jc),jb)       &
                    &              * opt_rho_ic(jc,jk,jb)


                  ! Courant number - 1
                  z_cfl_p(jc,jk,jb) = ( z_cfl_p(jc,jk,jb) - 1._wp )           &
                    &                * p_cellhgt_mc_now(jc,jk_int(jc),jb)     &
                    &                * p_rcellhgt_mc_now(jc,jk_int(jc)+1,jb)

                  ! Account for number of cells which we need to step upwind
                  ! for calculating the fractional flux
                  jk_int(jc) = jk_int(jc) + 1
                ENDIF

              ENDDO  ! end loop over cells

            ENDDO  ! end while loop



            ! if CFL < 1, everywhere in (1:nlen,jk,jb) then use faster standard
            ! calculation of (fractional) upwind fluxes
            IF (cfl_counter == 0) THEN

              DO jc = i_startidx, i_endidx

                ! linear extrapolated values
                ! first (of cell above) (case of weta > 0; w < 0; physical downwelling)
                z_lext_1(jc,jk) = p_cc(jc,jkm1,jb) - 0.5_wp                         &
                  &             * z_grad(jc,jkm1,jb) * p_cellhgt_mc_now(jc,jkm1,jb) &
                  &             * (1._wp - MIN(1._wp,z_cfl_m(jc,jk,jb)))

                ! second (of cell below) (case of weta < 0; w > 0; physical upwelling)
                z_lext_2(jc,jk) = p_cc(jc,jk,jb) + 0.5_wp                         &
                  &             * z_grad(jc,jk,jb) *  p_cellhgt_mc_now(jc,jk,jb)  &
                  &             * (1._wp - MIN(1._wp,z_cfl_p(jc,jk,jb)))

                !
                ! calculate vertical tracer flux
                !
                p_upflux(jc,jk,jb) =                                 &
                  &  laxfr_upflux( p_mflx_contra_v(jc,jk,jb),        &
                  &                z_lext_2(jc,jk), z_lext_1(jc,jk) )

              END DO ! end loop over cells

            ELSE

              ! if CFL > 1 anywhere in (1:nlen,jk,jb), then the following somewhat
              ! more expensive calculation of the sum of integer and fractional upwind
              ! fluxes is needed.
              DO jc = i_startidx, i_endidx
                IF (p_mflx_contra_v(jc,jk,jb) <= 0._wp ) THEN

                  ! Sum of Integer fluxes and fractional flux
                  ! if w(jc,jk,jb) <= 0._wp; physical downwelling
                  p_upflux(jc,jk,jb)= ( z_iflx_m(jc,jk)                           &
                    &             - opt_rho_ic(jc,jk,jb)                          &
                    &             * p_cellhgt_mc_now(jc,jkm_int(jc),jb)           &
                    &             * z_cfl_m(jc,jk,jb) * ( p_cc(jc,jkm_int(jc),jb) &
                    &             - 0.5_wp * z_grad(jc,jkm_int(jc),jb)            &
                    &             * p_cellhgt_mc_now(jc,jkm_int(jc),jb)           &
                    &             * (1._wp - MIN(1._wp,z_cfl_m(jc,jk,jb)))) )     &
                    &             / p_dtime

                ELSE

                  ! Sum of Integer fluxes and fractional flux
                  ! if w(jc,jk,jb) > 0._wp; physical upwelling
                  p_upflux(jc,jk,jb) = ( z_iflx_p(jc,jk)                          &
                    &             + opt_rho_ic(jc,jk,jb)                          &
                    &             * p_cellhgt_mc_now(jc,jk_int(jc),jb)            &
                    &             * z_cfl_p(jc,jk,jb) * ( p_cc(jc,jk_int(jc),jb)  &
                    &             + 0.5_wp * z_grad(jc,jk_int(jc),jb)             &
                    &             * p_cellhgt_mc_now(jc,jk_int(jc),jb)            &
                    &             * (1._wp - MIN(1._wp,z_cfl_p(jc,jk,jb)))) )     &
                    &             / p_dtime

                ENDIF

              ENDDO ! end loop over cells

            ENDIF

          ENDDO ! end loop over vertical levels

          !
          ! set upper and lower boundary condition
          !
          CALL set_bc_vadv(p_upflux(:,slev+1,jb),            &! in
            &              p_mflx_contra_v(:,slev+1,jb),     &! in
            &              p_mflx_contra_v(:,slev  ,jb),     &! in
            &              p_iubc_adv, i_startidx, i_endidx, &! in
            &              zparent_topflx(:,jb),             &! in
            &              p_upflux(:,slev,jb),              &! out
            &              p_upflux(:,nlevp1,jb)             )! out

        ENDDO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

      ENDIF ! IF (iequations /= 3)



    CASE( imuscl_v )
    !
    ! 3. extrapolation using piecewise linear approx. of the transported
    ! quantity to the edge and finally, calculation of the upwind fluxes
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_lext_1,z_lext_2,jkm1)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend )


        z_lext_1(i_startidx:i_endidx,slev) = p_cc(i_startidx:i_endidx,slev,jb)
        z_lext_2(i_startidx:i_endidx,slev) = p_cc(i_startidx:i_endidx,slev,jb)

        DO jk = slevp1, nlev

          ! index of top half level
          jkm1 = jk - 1

          DO jc = i_startidx, i_endidx

            ! linear extrapolated values
            ! first (of cell above) (case of weta > 0; w < 0; physical downwelling)
            z_lext_1(jc,jk) = p_cc(jc,jkm1,jb) + coeff_grid * 0.5_wp            &
              &             * z_grad(jc,jkm1,jb) * p_cellhgt_mc_now(jc,jkm1,jb) &
              &             * (1._wp - z_cfl_m(jc,jk,jb))


            ! second (of cell below) (case of weta < 0; w > 0; physical upwelling)
            z_lext_2(jc,jk) = p_cc(jc,jk,jb) - coeff_grid * 0.5_wp           &
              &             * z_grad(jc,jk,jb) * p_cellhgt_mc_now(jc,jk,jb)  &
              &             * (1._wp - z_cfl_p(jc,jk,jb))


            !
            ! calculate vertical tracer flux
            !
            p_upflux(jc,jk,jb) =                                  &
              &  laxfr_upflux_v( p_mflx_contra_v(jc,jk,jb),       &
              &                z_lext_1(jc,jk), z_lext_2(jc,jk),  &
              &                coeff_grid )

          END DO ! end loop over cells

        ENDDO ! end loop over vertical levels


        !
        ! set upper and lower boundary condition
        !
        CALL set_bc_vadv(p_upflux(:,slev+1,jb),            &! in
          &              p_mflx_contra_v(:,slev+1,jb),     &! in
          &              p_mflx_contra_v(:,slev  ,jb),     &! in
          &              p_iubc_adv, i_startidx, i_endidx, &! in
          &              zparent_topflx(:,jb),             &! in
          &              p_upflux(:,slev,jb),              &! out
          &              p_upflux(:,nlevp1,jb)             )! out

      ENDDO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

      !
      ! 6. If desired, apply a flux limiter to limit computed fluxes.
      !    These flux limiters are based on work by Zalesak (1979)
      !
      IF (iequations /= 3) THEN
        IF (p_itype_vlimit == ifluxl_vpd) THEN
          ! positive-definite (pd) limiter
          CALL vflx_limiter_pd_ha( p_patch, p_dtime, p_cc, p_upflux,     & !in,inout
            &                   opt_rlstart=i_rlstart, opt_rlend=i_rlend,& !in
            &                   opt_slev=slev                            ) !in 
        ENDIF
      ELSE
        IF (p_itype_vlimit == ifluxl_vpd) THEN
          ! positive-definite (pd) limiter
          CALL vflx_limiter_pd( p_patch, p_dtime, p_cc, p_upflux,     & !in,inout
            &                opt_rlstart=i_rlstart, opt_rlend=i_rlend,& !in
            &                opt_slev=slev                            ) !in
        ENDIF
      ENDIF

    END SELECT


    IF ( ld_cleanup ) THEN
      ! deallocate temporary arrays for layer thickness, weights and
      ! Courant numbers
      DEALLOCATE( z_delp_c1, z_delp_c2, z_alpha_1,       &
      &           z_gmoment, z_cfl_m0, z_cfl_p0, STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                     &
          &  'deallocation for z_delp_c1, z_delp_c2, '// &
          &  'z_alpha_1, z_gmoment, z_cfl_m0, z_cfl_p0 failed' )
      ENDIF
    END IF


  END SUBROUTINE upwind_vflux_muscl



  !-------------------------------------------------------------------------
  !>
  !! The third order PPM scheme
  !!
  !! Calculation of time averaged vertical tracer fluxes using the third
  !! order PPM scheme.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2009-08-12)
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized to height based vertical coordinate systems. Included
  !!   parameter coeff_grid, in order to apply the same code to either a
  !!   pressure based or height based vertical coordinate system.
  !! Modification by Daniel Reinert, DWD (2011-01-17)
  !! - added optional parameter opt_lout_edge which will provide the 
  !!   reconstructed 'edge' value.
  !!
  !
  ! !LITERATURE
  ! - Colella and Woodward (1984), JCP, 54, 174-201
  ! - Carpenter et al. (1989), MWR, 118, 586-612
  ! - Lin and Rood (1996), MWR, 124, 2046-2070
  !
  SUBROUTINE upwind_vflux_ppm( p_patch, p_cc, p_iubc_adv, p_mflx_contra_v,  &
    &                      p_w_contra, p_dtime, p_itype_vlimit,             &
    &                      p_cellhgt_mc_now, p_rcellhgt_mc_now, p_upflux,   &
    &                      opt_lout_edge, opt_topflx_tra, opt_slev,         &
    &                      opt_rlstart, opt_rlend )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_vflux: upwind_vflux_ppm'

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::  &    !< advected cell centered variable
      &  p_cc(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN)  ::   &   !< selects upper boundary condition
      &  p_iubc_adv

    REAL(wp), INTENT(IN) ::  &    !< contravariant vertical mass flux
      &  p_mflx_contra_v(:,:,:)   !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< contravariant vertical velocity
      &  p_w_contra(:,:,:)        !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) :: p_dtime  !< time step


    REAL(wp), INTENT(IN) ::  &    !< layer thickness at cell center at time n
      &  p_cellhgt_mc_now(:,:,:)  !< dim: (nproma,nlev,nblks_c)


    REAL(wp), INTENT(IN) ::  &    !< reciprocal of layer thickness at cell center
      &  p_rcellhgt_mc_now(:,:,:) !< at time n; dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(OUT) :: &    !< output field, containing the tracer mass flux
      &  p_upflux(:,:,:)          !< or the reconstructed edge value
                                  !< dim: (nproma,nlevp1,nblks_c)

    INTEGER, INTENT(IN) :: p_itype_vlimit  !< parameter to select the limiter for
                                           !< vertical transport

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: output edge value (.TRUE.),
     & opt_lout_edge                   !< or the flux across the edge 
                                       !< (.FALSE./not specified)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< vertical tracer flux at upper boundary 
      &  opt_topflx_tra(:,:)            !< dim: (nproma,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    LOGICAL  :: l_out_edgeval          !< corresponding local variable; default 
                                       !< .FALSE. i.e. output flux across the edge

    REAL(wp) :: &                      !< face values of transported field
      &  z_face(nproma,p_patch%nlevp1,p_patch%nblks_c)

    REAL(wp) :: &                      !< face value (upper face)
      &  z_face_up(nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp) :: &                      !< face value (lower face)
      &  z_face_low(nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp) :: &                      !< linear extrapolation value 1
      &  z_lext_1(nproma,p_patch%nlevp1)
 
    REAL(wp) :: &                      !< linear extrapolation value 2
      &  z_lext_2(nproma,p_patch%nlevp1)
 
    REAL(wp) :: &                      !< CFL number (weta>0, w<0)
      &  z_cfl_m(nproma,p_patch%nlevp1,p_patch%nblks_c) !< (nproma,nlevp1,p_patch%nblks_c)
                                                         
    REAL(wp) :: &                      !< CFL number (weta<0, w>0)
      &  z_cfl_p(nproma,p_patch%nlevp1,p_patch%nblks_c) !< (nproma,nlevp1,p_patch%nblks_c)
                                                       
    REAL(wp) :: &                      !< monotonized slope
      &  z_slope(nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp) ::  &                             !< necessary, to make this routine
     &  zparent_topflx(nproma,p_patch%nblks_c) !< compatible to the hydrost. core

    REAL(wp) :: z_slope_u, z_slope_l   !< one-sided slopes
    REAL(wp) :: z_delta_m, z_delta_p   !< difference between lower and upper face value
                                       !< for weta >0 and weta <0
    REAL(wp) :: z_a11, z_a12           !< 1/6 * a6,i (see Colella and Woodward (1984))
    REAL(wp) :: z_weta_dt              !< weta times p_dtime

    INTEGER  :: slev, slevp1           !< vertical start level and start level +1
    INTEGER  :: nlev, nlevp1           !< number of full and half levels
    INTEGER  :: ikm1, ikp1, ikp1_ic, & !< vertical level minus and plus one, plus two
      &  ikp2
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: jc, jk, jb              !< index of cell, vertical level and block

!DR    REAL(wp) :: cfl_m_max, cfl_p_max !< maximum Courant number from both sides
    !-----------------------------------------------------------------------

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev  = opt_slev
      slevp1= opt_slev + 1
    ELSE
      slev  = 1
      slevp1= 2
    END IF

    IF ( PRESENT(opt_lout_edge) ) THEN
      l_out_edgeval = opt_lout_edge
    ELSE
      l_out_edgeval = .FALSE.
    ENDIF

    IF ( PRESENT(opt_topflx_tra) ) THEN
      zparent_topflx(:,:) = opt_topflx_tra(:,:)
    ELSE
      zparent_topflx(:,:) = 0._wp
    ENDIF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c-1
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell_int
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)


    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1


    !
    ! advection is done with an upwind scheme and a piecwise parabolic
    ! approx. of the subgrid distribution is used.
    !
    ! 3 options:  standard without limiter
    !             standard with semi-monotone or monotone limiter
    !             special version with limiter which handles CFL >1

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

    !
    ! 1. Calculate Courant number for weta>0 (w<0) and weta<0 (w>0)
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,ikm1,z_weta_dt,ikp1_ic,ikp1, &
!$OMP            z_slope_u,z_slope_l,ikp2)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

      ! Courant number at top
      z_cfl_p(i_startidx:i_endidx,slev,jb) = 0._wp
      z_cfl_m(i_startidx:i_endidx,slev,jb) = 0._wp
      ! Courant number at bottom
      z_cfl_p(i_startidx:i_endidx,nlevp1,jb) = 0._wp
      z_cfl_m(i_startidx:i_endidx,nlevp1,jb) = 0._wp

      DO jk = slevp1, nlev

        ! index of top half level
        ikm1 = jk - 1

        DO jc = i_startidx, i_endidx

          ! Calculate local Courant number at half levels
          ! z_cfl_m for weta >0 (w <0)
          ! z_cfl_p for weta <0 (w >0)
          z_weta_dt = ABS(p_w_contra(jc,jk,jb)) * p_dtime

          z_cfl_m(jc,jk,jb) = z_weta_dt * p_rcellhgt_mc_now(jc,ikm1,jb)

          z_cfl_p(jc,jk,jb) = z_weta_dt * p_rcellhgt_mc_now(jc,jk,jb)

        END DO ! end loop over cells

      ENDDO ! end loop over vertical levels


      !
      ! 2. Calculate monotonized slope
      !
      z_slope(i_startidx:i_endidx,slev,jb) = 0._wp

      DO jk = slevp1, nlev

        ! index of top half level
        ikm1    = jk - 1
        ! index of bottom half level
        ikp1_ic = jk + 1
        ikp1    = MIN( ikp1_ic, nlev )

        DO jc = i_startidx, i_endidx

          z_slope_u = 2._wp * (p_cc(jc,jk,jb) - p_cc(jc,ikm1,jb))
          z_slope_l = 2._wp * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))

          IF ((z_slope_u * z_slope_l) .GT. 0._wp) THEN

            z_slope(jc,jk,jb) = ( p_cellhgt_mc_now(jc,jk,jb)                             &
              &  / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)            &
              &  + p_cellhgt_mc_now(jc,ikp1,jb)) )                                       &
              &  * ( (2._wp * p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)) &
              &  / (p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,jk,jb))           &
              &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))                                   &
              &  + (p_cellhgt_mc_now(jc,jk,jb) + 2._wp * p_cellhgt_mc_now(jc,ikp1,jb))   &
              &  / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb))           &
              &  * (p_cc(jc,jk,jb) - p_cc(jc,ikm1,jb)) )

            z_slope(jc,jk,jb) = SIGN(                                            &
              &  MIN( ABS(z_slope(jc,jk,jb)), ABS(z_slope_u), ABS(z_slope_l) ),  &
              &    z_slope(jc,jk,jb))

          ELSE

            z_slope(jc,jk,jb) = 0._wp

          ENDIF

        END DO ! end loop over cells

      END DO ! end loop over vertical levels

      !
      ! 3. reconstruct face values at vertical half-levels
      !

      ! Boundary values for two highest and lowest half-levels
      !
      ! for faces k=slevp1 and k=nlevp1-1 reconstructed face values are calculated by
      ! interpolating a quadratic (instead of quartic) polynomial through 3
      ! values of the indefinite integral A=\int_{\eta_{0}}^{\eta}q\,\mathrm{d}\eta
      !
      ! for faces k=slev and k=nlevp1 a zero gradient condition is assumed and the
      ! face values are set to the tracer values of the corresponding cell centers
      !
      DO jc = i_startidx, i_endidx

        z_face(jc,slevp1,jb) = p_cc(jc,slev,jb)*(1._wp - (p_cellhgt_mc_now(jc,slev,jb)&
          &       / p_cellhgt_mc_now(jc,slevp1,jb))) + (p_cellhgt_mc_now(jc,slev,jb)  &
          &       /(p_cellhgt_mc_now(jc,slev,jb) + p_cellhgt_mc_now(jc,slevp1,jb)))   &
          &       * ((p_cellhgt_mc_now(jc,slev,jb) / p_cellhgt_mc_now(jc,slevp1,jb))  &
          &       * p_cc(jc,slev,jb) + p_cc(jc,slevp1,jb))

        z_face(jc,nlev,jb) = p_cc(jc,nlev-1,jb)*( 1._wp                               &
          &       - (p_cellhgt_mc_now(jc,nlev-1,jb) / p_cellhgt_mc_now(jc,nlev,jb)))  &
          &       + (p_cellhgt_mc_now(jc,nlev-1,jb)/(p_cellhgt_mc_now(jc,nlev-1,jb)   &
          &       + p_cellhgt_mc_now(jc,nlev,jb))) * ((p_cellhgt_mc_now(jc,nlev-1,jb) &
          &       / p_cellhgt_mc_now(jc,nlev,jb)) * p_cc(jc,nlev-1,jb)                &
          &       + p_cc(jc,nlev,jb))

        z_face(jc,slev,jb)   = p_cc(jc,slev,jb)
        z_face(jc,nlevp1,jb) = p_cc(jc,nlev,jb)

      ENDDO


      DO jk = slevp1, nlev-2

        ! index of top half level
        ikm1 = jk - 1
        ! index of bottom half level
        ikp1 = jk + 1
        ikp2 = jk + 2

        DO jc = i_startidx, i_endidx

          z_face(jc,ikp1,jb) = p_cc(jc,jk,jb) + (p_cellhgt_mc_now(jc,jk,jb)           &
            &  / (p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb)))         &
            &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))                                  &
            &  + (1._wp/(p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)    &
            &  + p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,ikp2,jb)))        &
            &  * ( (2._wp * p_cellhgt_mc_now(jc,ikp1,jb) * p_cellhgt_mc_now(jc,jk,jb) &
            &  / (p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb)))         &
            &  * ( (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb))        &
            &  / (2._wp*p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb))    &
            &  - (p_cellhgt_mc_now(jc,ikp2,jb) + p_cellhgt_mc_now(jc,ikp1,jb))        &
            &  / (2._wp*p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,jk,jb)) )  &
            &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb)) - p_cellhgt_mc_now(jc,jk,jb)     &
            &  * z_slope(jc,ikp1,jb) * (p_cellhgt_mc_now(jc,ikm1,jb)                  &
            &  + p_cellhgt_mc_now(jc,jk,jb)) / (2._wp*p_cellhgt_mc_now(jc,jk,jb)      &
            &  + p_cellhgt_mc_now(jc,ikp1,jb)) + p_cellhgt_mc_now(jc,ikp1,jb)         &
            &  * z_slope(jc,jk,jb) * (p_cellhgt_mc_now(jc,ikp1,jb)                    &
            &  + p_cellhgt_mc_now(jc,ikp2,jb)) / (p_cellhgt_mc_now(jc,jk,jb)          &
            &  + 2._wp*p_cellhgt_mc_now(jc,ikp1,jb)) )

        END DO ! end loop over cells

      END DO ! end loop over vertical levels

    END DO
!$OMP END DO
!$OMP END PARALLEL

!DR  cfl_p_max=MAXVAL(z_cfl_p(:,34:54,:))
!DR  cfl_m_max=MAXVAL(z_cfl_m(:,34:54,:))
!DR  WRITE(message_text,'(a,e16.8)') 'CFL_P MAX', cfl_p_max
!DR  CALL message('vertical_upwind_flux',message_text)
!DR  WRITE(message_text,'(a,e16.8)') 'CFL_M MAX', cfl_m_max
!DR  CALL message('vertical_upwind_flux',message_text)


    !
    ! 4. Limitation of first guess parabola (which is based on z_face)
    ! Note that z_face_up(k) does not need to equal z_face_low(k-1) after
    ! the limitation procedure.
    ! Therefore 2 additional fields z_face_up and z_face_low are
    ! introduced.
    !
    IF (p_itype_vlimit == islopel_vsm) THEN
      ! semi-monotonic (sm) limiter
      CALL v_ppm_slimiter_sm( p_patch, p_cc, z_face,                  & !in
        &                   z_face_up, z_face_low,                    & !inout
        &                   opt_rlstart=i_rlstart, opt_rlend=i_rlend, & !in
        &                   opt_slev=slev                             ) !in
    ELSE IF (p_itype_vlimit == islopel_vm) THEN
      ! monotonic (mo) limiter
      CALL v_ppm_slimiter_mo( p_patch, p_cc, z_face, z_slope,         & !in
        &                   z_face_up, z_face_low,                    & !inout
        &                   opt_rlstart=i_rlstart, opt_rlend=i_rlend, & !in
        &                   opt_slev=slev                             ) !in
    ELSE
      ! simply copy face values to 'face_up' and 'face_low' arrays
!$OMP PARALLEL
!$OMP DO PRIVATE(jk,ikp1,jb,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend )

        DO jk = slev, nlev
          ! index of bottom half level
          ikp1 = jk + 1
          z_face_up(i_startidx:i_endidx,jk,jb)  = z_face(i_startidx:i_endidx,jk,jb)
          z_face_low(i_startidx:i_endidx,jk,jb) = z_face(i_startidx:i_endidx,ikp1,jb)
        ENDDO
      ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
    ENDIF



!$OMP PARALLEL
    IF ( l_out_edgeval ) THEN

      ! 5a. Compute edge value of advected quantity by applying a 
      !     piecewise parabolic reconstruction

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_lext_1,z_lext_2,ikm1,z_delta_m, &
!$OMP            z_delta_p,z_a11,z_a12)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend )

        ! for the time being, when computing the edge value, we set 
        ! the top and bottom values to 0
        p_upflux(i_startidx:i_endidx,  slev,jb) = 0._wp
        p_upflux(i_startidx:i_endidx,nlevp1,jb) = 0._wp

        z_lext_1(i_startidx:i_endidx,slev)   = p_cc(i_startidx:i_endidx,slev,jb)
        z_lext_2(i_startidx:i_endidx,slev)   = p_cc(i_startidx:i_endidx,slev,jb)
        z_lext_1(i_startidx:i_endidx,nlevp1) = p_cc(i_startidx:i_endidx,nlev,jb)
        z_lext_2(i_startidx:i_endidx,nlevp1) = p_cc(i_startidx:i_endidx,nlev,jb)

        DO jk = slevp1, nlev

          ! index of top half level
          ikm1 = jk -1


          DO jc = i_startidx, i_endidx

            ! linear extrapolated values
            ! for the height based coordinate system multiplication by coeff_grid
            ! is not necessary due to compensating (-) signs.
            ! first (of cell above) (case of w < 0; weta > 0)
            z_delta_m = z_face_low(jc,ikm1,jb) - z_face_up(jc,ikm1,jb)
            z_a11     = p_cc(jc,ikm1,jb)                                  &
              &       - 0.5_wp * (z_face_low(jc,ikm1,jb) + z_face_up(jc,ikm1,jb))

            z_lext_1(jc,jk) = p_cc(jc,ikm1,jb)                            &
              &  + (0.5_wp * z_delta_m * (1._wp - z_cfl_m(jc,jk,jb)))     &
              &  - z_a11*(1._wp - 3._wp*z_cfl_m(jc,jk,jb)                 &
              &  + 2._wp*z_cfl_m(jc,jk,jb)*z_cfl_m(jc,jk,jb))


            ! second (of cell below) (case of w > 0; weta < 0)
            z_delta_p = z_face_low(jc,jk,jb) - z_face_up(jc,jk,jb)
            z_a12     = p_cc(jc,jk,jb)                                    &
              &       - 0.5_wp * (z_face_low(jc,jk,jb) + z_face_up(jc,jk,jb))

            z_lext_2(jc,jk) = p_cc(jc,jk,jb)                              &
              &  - (0.5_wp * z_delta_p * (1._wp - z_cfl_p(jc,jk,jb)))     &
              &  - z_a12*(1._wp - 3._wp*z_cfl_p(jc,jk,jb)                 &
              &  + 2._wp*z_cfl_p(jc,jk,jb)*z_cfl_p(jc,jk,jb))

            !
            ! calculate 'edge value' of advected quantity
            !
            p_upflux(jc,jk,jb) =                                         &
              &  laxfr_upflux_v( SIGN(1._wp, p_mflx_contra_v(jc,jk,jb)), &
              &                z_lext_1(jc,jk), z_lext_2(jc,jk),         &
              &                -1.0_wp )  ! for test purposes only in nh model

            ! sign of the edge value
            p_upflux(jc,jk,jb) =  p_upflux(jc,jk,jb)*SIGN(1._wp, p_mflx_contra_v(jc,jk,jb))

          END DO ! end loop over cells

        ENDDO ! end loop over vertical levels

      ENDDO ! end loop over blocks
!$OMP END DO

    ELSE
    !
    ! 5b. extrapolation using piecewise parabolic approx. of the transported
    ! quantity to the edge and finally, calculation of the upwind fluxes
    !

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_lext_1,z_lext_2,ikm1,z_delta_m, &
!$OMP            z_delta_p,z_a11,z_a12)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend )

        z_lext_1(i_startidx:i_endidx,slev)   = p_cc(i_startidx:i_endidx,slev,jb)
        z_lext_2(i_startidx:i_endidx,slev)   = p_cc(i_startidx:i_endidx,slev,jb)
        z_lext_1(i_startidx:i_endidx,nlevp1) = p_cc(i_startidx:i_endidx,nlev,jb)
        z_lext_2(i_startidx:i_endidx,nlevp1) = p_cc(i_startidx:i_endidx,nlev,jb)

        DO jk = slevp1, nlev

          ! index of top half level
          ikm1 = jk -1

          DO jc = i_startidx, i_endidx

            ! linear extrapolated values
            ! for the height based coordinate system multiplication by coeff_grid
            ! is not necessary due to compensating (-) signs.
            ! first (of cell above) (case of w < 0; weta > 0)
            z_delta_m = z_face_low(jc,ikm1,jb) - z_face_up(jc,ikm1,jb)
            z_a11     = p_cc(jc,ikm1,jb)                                  &
              &       - 0.5_wp * (z_face_low(jc,ikm1,jb) + z_face_up(jc,ikm1,jb))

            z_lext_1(jc,jk) = p_cc(jc,ikm1,jb)                            &
              &  + (0.5_wp * z_delta_m * (1._wp - z_cfl_m(jc,jk,jb)))     &
              &  - z_a11*(1._wp - 3._wp*z_cfl_m(jc,jk,jb)                 &
              &  + 2._wp*z_cfl_m(jc,jk,jb)*z_cfl_m(jc,jk,jb))


            ! second (of cell below) (case of w > 0; weta < 0)
            z_delta_p = z_face_low(jc,jk,jb) - z_face_up(jc,jk,jb)
            z_a12     = p_cc(jc,jk,jb)                                    &
              &       - 0.5_wp * (z_face_low(jc,jk,jb) + z_face_up(jc,jk,jb))

            z_lext_2(jc,jk) = p_cc(jc,jk,jb)                              &
              &  - (0.5_wp * z_delta_p * (1._wp - z_cfl_p(jc,jk,jb)))     &
              &  - z_a12*(1._wp - 3._wp*z_cfl_p(jc,jk,jb)                 &
              &  + 2._wp*z_cfl_p(jc,jk,jb)*z_cfl_p(jc,jk,jb))

            !
            ! calculate vertical tracer flux
            !
            p_upflux(jc,jk,jb) =                                  &
              &  laxfr_upflux_v( p_mflx_contra_v(jc,jk,jb),       &
              &                z_lext_1(jc,jk), z_lext_2(jc,jk),  &
              &                coeff_grid )

          END DO ! end loop over cells

        ENDDO ! end loop over vertical levels


        !
        ! set upper and lower boundary condition
        !
        CALL set_bc_vadv(p_upflux(:,slev+1,jb),            &! in
          &              p_mflx_contra_v(:,slev+1,jb),     &! in
          &              p_mflx_contra_v(:,slev  ,jb),     &! in
          &              p_iubc_adv, i_startidx, i_endidx, &! in
          &              zparent_topflx(:,jb),             &! in
          &              p_upflux(:,slev,jb),              &! out
          &              p_upflux(:,nlevp1,jb)             )! out

      ENDDO ! end loop over blocks
!$OMP END DO

      ENDIF
!$OMP END PARALLEL

      !
      ! 6. If desired, apply a flux limiter to limit computed fluxes.
      !    These flux limiters are based on work by Zalesak (1979)
      !
      IF (iequations /= 3) THEN
        IF (p_itype_vlimit == ifluxl_vpd) THEN
          ! positive-definite (pd) limiter
          CALL vflx_limiter_pd_ha( p_patch, p_dtime, p_cc, p_upflux,    & !in,inout
            &                 opt_rlstart=i_rlstart, opt_rlend=i_rlend, & !in
            &                 opt_slev=slev                             ) !in
        ENDIF
      ELSE
        IF (p_itype_vlimit == ifluxl_vpd) THEN
          ! positive-definite (pd) limiter
          CALL vflx_limiter_pd( p_patch, p_dtime, p_cc, p_upflux,       & !in,inout
            &                opt_rlstart=i_rlstart, opt_rlend=i_rlend,  & !in
            &                opt_slev=slev                              ) !in
        ENDIF
      ENDIF


  END SUBROUTINE upwind_vflux_ppm



  !-------------------------------------------------------------------------
  !>
  !! The third order PPM scheme for large time steps (CFL>1)
  !!
  !! Calculation of time averaged vertical tracer fluxes or tracer edge 
  !! values using the third order PPM scheme. This scheme can handle 
  !! large time steps (i.e. CFL>1)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-01-14)
  !!
  !
  ! !LITERATURE
  ! - Colella and Woodward (1984), JCP, 54, 174-201
  ! - Carpenter et al. (1989), MWR, 118, 586-612
  ! - Lin and Rood (1996), MWR, 124, 2046-2070 (CFL-independent version)
  !
  SUBROUTINE upwind_vflux_ppm_cfl( p_patch, p_cc, p_iubc_adv, p_mflx_contra_v, &
    &                      p_dtime,  ld_compute, ld_cleanup, p_itype_vlimit,   &
    &                      p_cellhgt_mc_now, p_cellmass_now, p_upflux,         &
    &                      opt_lout_edge, opt_topflx_tra, opt_slev,            &
    &                      opt_rlstart, opt_rlend )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_vflux: upwind_vflux_ppm_cfl'

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::  &    !< advected cell centered variable
      &  p_cc(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN)  ::   &   !< selects upper boundary condition
      &  p_iubc_adv

    REAL(wp), INTENT(INOUT) ::  & !< contravariant vertical mass flux
      &  p_mflx_contra_v(:,:,:)   !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< time step
      &  p_dtime

    LOGICAL, INTENT(IN)  ::  &    !< flag, if .TRUE. compute geometrical terms
      &  ld_compute

    LOGICAL, INTENT(IN)  ::  &    !< flag, if .TRUE. clean up geometrical terms
      &  ld_cleanup

    INTEGER, INTENT(IN)  ::  &    !< parameter to select the limiter for
      &  p_itype_vlimit           !< vertical transport

    REAL(wp), INTENT(IN) ::  &    !< layer thickness at cell center at time n
      &  p_cellhgt_mc_now(:,:,:)  !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< NH: density weighted cell height at full levels
      &  p_cellmass_now(:,:,:)    !< at time step n [kg/m**2]
                                  !< HA: pressure thickness for full levels 
                                  !< at time step n [Pa]
                                  !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(OUT) :: &    !< output field, containing the tracer mass flux
      &  p_upflux(:,:,:)          !< or the reconstructed edge value
                                  !< dim: (nproma,nlevp1,nblks_c)

    LOGICAL, INTENT(IN), OPTIONAL ::  & !< optional: output edge value (.TRUE.),
     & opt_lout_edge                    !< or the flux across the edge 
                                        !< (.FALSE./not specified)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< vertical tracer flux at upper boundary 
      &  opt_topflx_tra(:,:)            !< dim: (nproma,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    LOGICAL  :: l_out_edgeval     !< corresponding local variable; default 
                                  !< .FALSE. i.e. output flux across the edge

    REAL(wp) :: &                 !< face values of transported field
      &  z_face(nproma,p_patch%nlevp1,p_patch%nblks_c)

    REAL(wp) :: &                 !< face value (upper face)
      &  z_face_up(nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp) :: &                 !< face value (lower face)
      &  z_face_low(nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp) :: &                 !< integer fluxes for w>0 (weta<0)
      &  z_iflx_p(nproma,p_patch%nlevp1)

    REAL(wp) :: &                 !< integer fluxes for w<0 (weta>0)
      &  z_iflx_m(nproma,p_patch%nlevp1)

    REAL(wp) :: &                 !< monotonized slope
      &  z_slope(nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp) :: z_slope_u, z_slope_l     !< one-sided slopes
    REAL(wp) :: z_delta_p, z_delta_m     !< difference between upper and lower face value
                                         !< for w>0 and w<0
    REAL(wp) :: z_a11, z_a12             !< 1/6 * a6,i (see Colella and Woodward (1984))

    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: ikm1, ikp1, ikp1_ic, &   !< vertical level minus and plus one, plus two
      &  ikp2
    INTEGER  :: slev, slevp1             !< vertical start level and start level +1
    INTEGER  :: nlev, nlevp1             !< number of full and half levels

    INTEGER  :: ji_p, ji_m               !< loop variable for index list
    INTEGER  :: ist                      !< status variable
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom

    INTEGER, PARAMETER :: nlist_max=3    !< maximum number of index lists
    INTEGER  :: nlist_p, nlist_m         !< list number
    INTEGER  :: nlist                    !< list loop variable

    REAL(wp), ALLOCATABLE, SAVE  ::   &  !< fractional (mass weighted) Courant number 
      &  z_cflfrac_p(:,:,:)              !< for w>0

    REAL(wp), ALLOCATABLE, SAVE  ::   &  !< fractional (mass weighted) Courant number 
      &  z_cflfrac_m(:,:,:)              !< for w<0

    INTEGER, ALLOCATABLE, SAVE  ::    &  !< Index and level lists and list dimensions 
      &  i_indlist_p(:,:,:),  &          !< for points with CFL>1,2,3
      &  i_indlist_m(:,:,:),  &
      &  i_levlist_p(:,:,:),  &
      &  i_levlist_m(:,:,:),  &
      &  i_listdim_p(:,:),    &
      &  i_listdim_m(:,:)

    INTEGER, ALLOCATABLE, SAVE  ::    &  !< shifted indices
      &  jk_int_p(:,:,:),             &  ! jk+s
      &  jk_int_m(:,:,:)                 ! jk-s, with shift index s

    INTEGER  :: jk_shift
    INTEGER  :: counter_p, counter_m, &  !< check whether any of the points has 
      &  counter_jip, counter_jim        !< CFL>nlist_p/m

    REAL(wp) ::   &                      !< Integer flux
      &  z_flx_int(nproma,p_patch%nlevp1,p_patch%nblks_c)

    REAL(wp) ::   &                      !< high order flux
      &  z_flx_frac_high(nproma,p_patch%nlevp1,p_patch%nblks_c)

!DR    REAL(wp)  ::   &                     !< low order flux (will be necessary if 
!DR      &  z_flx_frac_low(nproma,p_patch%nlevp1,p_patch%nblks_c) 
!DR                                         !< we implement the FCT version

    REAL(wp) ::   &                      !< maximum vertical Courant number
      &  max_cfl(nproma,p_patch%nlevp1,p_patch%nblks_c)

    REAL(wp) ::  &                              !< necessary, to make this routine
      &  zparent_topflx(nproma,p_patch%nblks_c) !< compatible to the hydrost. core 

    REAL(wp) ::   &                      !< dummy variable
      &  z_dummy

    !-----------------------------------------------------------------------


    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev  = opt_slev
      slevp1= opt_slev + 1
    ELSE
      slev  = 1
      slevp1= 2
    END IF

    IF ( PRESENT(opt_lout_edge) ) THEN
      l_out_edgeval = opt_lout_edge
    ELSE
      l_out_edgeval = .FALSE.
    ENDIF

    IF ( PRESENT(opt_topflx_tra) ) THEN
      zparent_topflx(:,:) = opt_topflx_tra(:,:)
    ELSE
      zparent_topflx(:,:) = 0._wp
    ENDIF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c-1
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell_int
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)


    !
    ! advection is done with an upwind scheme where a piecwise parabolic
    ! approx. of the subgrid distribution is used.
    !
    ! 3 options:  standard without limiter
    !             standard with semi-monotone or monotone limiter
    !             special version with limiter which handles CFL >1

    IF ( ld_compute ) THEN
      ! allocate temporary arrays 
      ALLOCATE( i_indlist_p(nproma*nlevp1,nlist_max,p_patch%nblks_c),  &
        &       i_indlist_m(nproma*nlevp1,nlist_max,p_patch%nblks_c),  &
        &       i_levlist_p(nproma*nlevp1,nlist_max,p_patch%nblks_c),  &
        &       i_levlist_m(nproma*nlevp1,nlist_max,p_patch%nblks_c),  &
        &       i_listdim_p(nlist_max,p_patch%nblks_c),                &
        &       i_listdim_m(nlist_max,p_patch%nblks_c),                &
        &       jk_int_p(nproma,nlev,p_patch%nblks_c),                 &
        &       jk_int_m(nproma,nlev,p_patch%nblks_c),                 &
        &       z_cflfrac_p(nproma,nlevp1,p_patch%nblks_c),            &
        &       z_cflfrac_m(nproma,nlevp1,p_patch%nblks_c), STAT=ist   )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                         &
          &  'allocation for i_indlist_p, i_indlist_m, i_levlist_p, '  //  &
          &  'i_levlist_m, i_listdim_p, i_listdim_m, jk_int_p, '       //  &
          &  'jk_int_m, z_cflfrac_p, z_cflfrac_m failed '  )
      ENDIF

      IF (msg_level >= 11) THEN
        ! otherwise an error occurs when running in debug mode
        max_cfl(:,:,:) = 0._wp
      ENDIF
    END IF


!$OMP PARALLEL
    !
    ! The contravariant mass flux should never exactly vanish
    !
    IF (l_out_edgeval) THEN
      coeff_grid=-1 ! assume that this is only used in nh model
!$OMP DO PRIVATE(jb,jk,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend )

        DO jk = slevp1, nlev
          p_mflx_contra_v(i_startidx:i_endidx,jk,jb) =                            &
          &              p_mflx_contra_v(i_startidx:i_endidx,jk,jb)               &
          &              + SIGN(dbl_eps,p_mflx_contra_v(i_startidx:i_endidx,jk,jb))
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF



    !
    ! 1. Compute density weighted (fractional) Courant number 
    !    for w<0 (weta>0) and w>0 (weta<0) and integer shift s
    !
    IF (ld_compute) THEN
!$OMP DO PRIVATE(jb,jk,jc,ikm1,i_startidx,i_endidx,z_dummy,nlist_p,nlist_m, &
!$OMP            counter_p,counter_m,counter_jip,counter_jim)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

      ! set start values for index list dimension
      i_listdim_p(1:nlist_max,jb) = 0
      i_listdim_m(1:nlist_max,jb) = 0

      ! (fractional) Courant number at top
      z_cflfrac_p(i_startidx:i_endidx,slev,jb) = 0._wp
      z_cflfrac_m(i_startidx:i_endidx,slev,jb) = 0._wp
      ! (fractional) Courant number at bottom
      z_cflfrac_p(i_startidx:i_endidx,nlevp1,jb) = 0._wp
      z_cflfrac_m(i_startidx:i_endidx,nlevp1,jb) = 0._wp


      !
      ! compute (fractional) Courant number
      !
      DO jk = slevp1, nlev

        ikm1 = jk-1

        DO jc = i_startidx, i_endidx

          ! initialize shift index
          jk_int_p(jc,jk,jb) = jk
          jk_int_m(jc,jk,jb) = ikm1

          z_dummy = p_dtime * ABS(p_mflx_contra_v(jc,jk,jb))

          ! compute Courant number
          z_cflfrac_p(jc,jk,jb) = z_dummy / p_cellmass_now(jc,jk,jb)
          z_cflfrac_m(jc,jk,jb) = z_dummy / p_cellmass_now(jc,ikm1,jb)

          max_cfl(jc,jk,jb) = MAX(z_cflfrac_p(jc,jk,jb),z_cflfrac_m(jc,jk,jb))
        ENDDO

      ENDDO


      ! If CFL>1 then split the CFL number into the fractional CFL number 
      ! and the index shift s.
      IF ( MAXVAL(max_cfl(i_startidx:i_endidx,slevp1:nlev,jb)) > 1._wp ) THEN

        DO jk = slevp1, nlev

          ! start construction of fractional Courant number
          z_cflfrac_p(i_startidx:i_endidx,jk,jb) = p_dtime                  &
             &              * ABS(p_mflx_contra_v(i_startidx:i_endidx,jk,jb))
          z_cflfrac_m(i_startidx:i_endidx,jk,jb) = p_dtime                  &
             &              * ABS(p_mflx_contra_v(i_startidx:i_endidx,jk,jb))

          ! initialize list number
          nlist_p = 0
          nlist_m = 0

          ! checks whether there exists any point with 'large CFL number'
          counter_p = 1
          counter_m = 1

          ! loop until no point has CFL > nlist_p, or nlist_p > nlist_max
          DO WHILE(counter_p > 0 .AND. nlist_p < nlist_max )

            ! get number of current list
            nlist_p     = nlist_p + 1
            ! re-initialize counter for CFL>nlist_p
            counter_p   = 0

            ! copy value from counter in vector form to counter in scalar form, 
            ! since otherwise the following DO-Loop will not vectorize.
            counter_jip = i_listdim_p(nlist_p,jb)

            DO jc = i_startidx, i_endidx

              IF ( z_cflfrac_p(jc,jk,jb) > p_cellmass_now(jc,jk_int_p(jc,jk,jb),jb) &
                &  .AND. jk_int_p(jc,jk,jb) <= nlev-1 ) THEN

                z_cflfrac_p(jc,jk,jb) = z_cflfrac_p(jc,jk,jb)                   &
                  &                   - p_cellmass_now(jc,jk_int_p(jc,jk,jb),jb)

                ! update shift index
                jk_int_p(jc,jk,jb) = jk_int_p(jc,jk,jb) + 1

                ! tests whether we need to loop once again
                counter_p = counter_p + 1

                ! Fill index lists with those points that need index shifts
                ! Note that we have to use a scalar counter instead of a vector, like
                ! i_listdim_p(nlist_p,jb). Otherwise this loop will not vectorize.  
                counter_jip = counter_jip + 1
                i_indlist_p(counter_jip,nlist_p,jb) = jc
                i_levlist_p(counter_jip,nlist_p,jb) = jk

              ELSE
                ! compute fractional Courant number
                z_cflfrac_p(jc,jk,jb) = z_cflfrac_p(jc,jk,jb)                    &
                 &                    / p_cellmass_now(jc,jk_int_p(jc,jk,jb),jb)

              ENDIF

            END DO ! end loop over cells

            ! store index of current index list
            ! after the last jk-loop this will be the list dimension
            i_listdim_p(nlist_p,jb) = counter_jip

          ENDDO  ! DO WHILE loop
 

          ! loop until no point has CFL>nlist_m, or nlist_m > nlist_max
          DO WHILE(counter_m > 0 .AND. nlist_m < nlist_max )

            ! set number of current list
            nlist_m     = nlist_m + 1
            ! re-initialize counter for CFL>nlist_m
            counter_m   = 0

            ! copy value from counter in vector form to counter in scalar form, 
            ! since otherwise the following DO-Loop will not vectorize.
            counter_jim = i_listdim_m(nlist_m,jb)


            DO jc = i_startidx, i_endidx

              IF ( z_cflfrac_m(jc,jk,jb) > p_cellmass_now(jc,jk_int_m(jc,jk,jb),jb) &
                & .AND. jk_int_m(jc,jk,jb) >= slevp1) THEN

                 z_cflfrac_m(jc,jk,jb) = z_cflfrac_m(jc,jk,jb)                   &
                  &                    - p_cellmass_now(jc,jk_int_m(jc,jk,jb),jb)

                ! Index shift
                jk_int_m(jc,jk,jb) = jk_int_m(jc,jk,jb) - 1

                ! tests whether we need to loop once again
                counter_m = counter_m + 1

                ! Fill index lists with those points that need index shifts
                ! Note that we have to use a scalar instead of a vector, like
                ! i_listdim_m(nlist_m,jb). Otherwise it will not vectorize  
                counter_jim = counter_jim + 1
                i_indlist_m(counter_jim,nlist_m,jb) = jc
                i_levlist_m(counter_jim,nlist_m,jb) = jk

              ELSE
                ! compute fractional Courant number
                z_cflfrac_m(jc,jk,jb) = z_cflfrac_m(jc,jk,jb)                   &
                  &                   / p_cellmass_now(jc,jk_int_m(jc,jk,jb),jb)
              ENDIF

            END DO ! end loop over cells

            ! store index of current index list.
            ! after the last jk-loop this will be the list dimension
            i_listdim_m(nlist_m,jb) = counter_jim

          ENDDO  ! DO WHILE loop


        ENDDO ! end loop over vertical levels

      END IF

    ENDDO ! end loop over blocks
!$OMP END DO NOWAIT

!$OMP SINGLE
    IF (msg_level >= 11) THEN ! print maximum vertical CFL number
      WRITE(message_text,'(e16.8)') MAXVAL(max_cfl(:,:,:))
      CALL message('maximum vertical CFL:',message_text)
    ENDIF ! msg_level >= 11
!$OMP END SINGLE NOWAIT
    END IF


    !
    ! 2. Compute monotonized slope
    !
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,ikm1,ikp1_ic,ikp1,z_slope_u, &
!$OMP            z_slope_l,ikp2)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

      z_slope(i_startidx:i_endidx,slev,jb) = 0._wp

      DO jk = slevp1, nlev

        ! index of top half level
        ikm1    = jk - 1
        ! index of bottom half level
        ikp1_ic = jk + 1
        ikp1    = MIN( ikp1_ic, nlev )

        DO jc = i_startidx, i_endidx

          z_slope_u = 2._wp * (p_cc(jc,jk,jb) - p_cc(jc,ikm1,jb))
          z_slope_l = 2._wp * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))

          IF ((z_slope_u * z_slope_l) .GT. 0._wp) THEN

            z_slope(jc,jk,jb) = ( p_cellhgt_mc_now(jc,jk,jb)                             &
              &  / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)            &
              &  + p_cellhgt_mc_now(jc,ikp1,jb)) )                                       &
              &  * ( (2._wp * p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)) &
              &  / (p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,jk,jb))           &
              &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))                                   &
              &  + (p_cellhgt_mc_now(jc,jk,jb) + 2._wp * p_cellhgt_mc_now(jc,ikp1,jb))   &
              &  / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb))           &
              &  * (p_cc(jc,jk,jb) - p_cc(jc,ikm1,jb)) )

            z_slope(jc,jk,jb) = SIGN(                                            &
              &  MIN( ABS(z_slope(jc,jk,jb)), ABS(z_slope_u), ABS(z_slope_l) ),  &
              &    z_slope(jc,jk,jb))

          ELSE

            z_slope(jc,jk,jb) = 0._wp

          ENDIF

        END DO

      END DO


      !
      ! 3. reconstruct face values at vertical half-levels
      !

      ! Boundary values for two highest and lowest half-levels
      !
      ! for faces k=slevp1 and k=nlevp1-1 reconstructed face values are calculated by
      ! interpolating a quadratic (instead of quartic) polynomial through 3
      ! values of the indefinite integral A=\int_{\eta_{0}}^{\eta}q\,\mathrm{d}\eta
      !
      ! for faces k=slev and k=nlevp1 a zero gradient condition is assumed and the
      ! face values are set to the tracer values of the corresponding cell centers
      !
      DO jc = i_startidx, i_endidx

        z_face(jc,slevp1,jb) = p_cc(jc,slev,jb)*(1._wp - (p_cellhgt_mc_now(jc,slev,jb)&
          &       / p_cellhgt_mc_now(jc,slevp1,jb))) + (p_cellhgt_mc_now(jc,slev,jb)  &
          &       /(p_cellhgt_mc_now(jc,slev,jb) + p_cellhgt_mc_now(jc,slevp1,jb)))   &
          &       * ((p_cellhgt_mc_now(jc,slev,jb) / p_cellhgt_mc_now(jc,slevp1,jb))  &
          &       * p_cc(jc,slev,jb) + p_cc(jc,slevp1,jb))

        z_face(jc,nlev,jb) = p_cc(jc,nlev-1,jb)*( 1._wp                               &
          &       - (p_cellhgt_mc_now(jc,nlev-1,jb) / p_cellhgt_mc_now(jc,nlev,jb)))  &
          &       + (p_cellhgt_mc_now(jc,nlev-1,jb)/(p_cellhgt_mc_now(jc,nlev-1,jb)   &
          &       + p_cellhgt_mc_now(jc,nlev,jb))) * ((p_cellhgt_mc_now(jc,nlev-1,jb) &
          &       / p_cellhgt_mc_now(jc,nlev,jb)) * p_cc(jc,nlev-1,jb)                &
          &       + p_cc(jc,nlev,jb))

        z_face(jc,slev,jb)   = p_cc(jc,slev,jb)
        z_face(jc,nlevp1,jb) = p_cc(jc,nlev,jb)

      ENDDO


      DO jk = slevp1, nlev-2

        ! index of top half level
        ikm1 = jk - 1
        ! index of bottom half level
        ikp1 = jk + 1
        ikp2 = jk + 2

        DO jc = i_startidx, i_endidx

          z_face(jc,ikp1,jb) = p_cc(jc,jk,jb) + (p_cellhgt_mc_now(jc,jk,jb)           &
            &  / (p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb)))         &
            &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))                                  &
            &  + (1._wp/(p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)    &
            &  + p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,ikp2,jb)))        &
            &  * ( (2._wp * p_cellhgt_mc_now(jc,ikp1,jb) * p_cellhgt_mc_now(jc,jk,jb) &
            &  / (p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb)))         &
            &  * ( (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb))        &
            &  / (2._wp*p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb))    &
            &  - (p_cellhgt_mc_now(jc,ikp2,jb) + p_cellhgt_mc_now(jc,ikp1,jb))        &
            &  / (2._wp*p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,jk,jb)) )  &
            &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb)) - p_cellhgt_mc_now(jc,jk,jb)     &
            &  * z_slope(jc,ikp1,jb) * (p_cellhgt_mc_now(jc,ikm1,jb)                  &
            &  + p_cellhgt_mc_now(jc,jk,jb)) / (2._wp*p_cellhgt_mc_now(jc,jk,jb)      &
            &  + p_cellhgt_mc_now(jc,ikp1,jb)) + p_cellhgt_mc_now(jc,ikp1,jb)         &
            &  * z_slope(jc,jk,jb) * (p_cellhgt_mc_now(jc,ikp1,jb)                    &
            &  + p_cellhgt_mc_now(jc,ikp2,jb)) / (p_cellhgt_mc_now(jc,jk,jb)          &
            &  + 2._wp*p_cellhgt_mc_now(jc,ikp1,jb)) )

        END DO

      END DO

    END DO
!$OMP END DO
!$OMP END PARALLEL



    !
    ! 4. Limitation of first guess parabola (which is based on z_face)
    ! Note that z_face_up(k) does not need to equal z_face_low(k-1) after
    ! the limitation procedure.
    ! Therefore 2 additional fields z_face_up and z_face_low are
    ! introduced.
    !
    IF (p_itype_vlimit == islopel_vsm) THEN
      ! semi-monotonic (sm) limiter
      CALL v_ppm_slimiter_sm( p_patch, p_cc, z_face,              & !in
        &                     z_face_up, z_face_low,              & !inout
        &                     opt_rlstart=i_rlstart,              & !in
        &                     opt_rlend=i_rlend, opt_slev=slev    ) !in
    ELSE IF (p_itype_vlimit == islopel_vm) THEN
      ! monotonic (mo) limiter
      CALL v_ppm_slimiter_mo( p_patch, p_cc, z_face, z_slope,     & !in
        &                     z_face_up, z_face_low,              & !inout
        &                     opt_rlstart=i_rlstart,              & !in
        &                     opt_rlend=i_rlend, opt_slev=slev    ) !in
    ELSE
      ! simply copy face values to 'face_up' and 'face_low' arrays
!$OMP PARALLEL
!$OMP DO PRIVATE(jk,ikp1,jb,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend )

        DO jk = slev, nlev
          ! index of bottom half level
          ikp1 = jk + 1
          z_face_up(i_startidx:i_endidx,jk,jb)  = z_face(i_startidx:i_endidx,jk,jb)
          z_face_low(i_startidx:i_endidx,jk,jb) = z_face(i_startidx:i_endidx,ikp1,jb)
        ENDDO
      ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
    ENDIF



    !
    ! 3. calculation of upwind fluxes. IF CFL > 1, the fluxes are the sum of
    !    integer-fluxes and a fractional flux. IF CFL <1 the fluxes are only
    !    comprised of the fractional flux. The fractional flux is calculated
    !    assuming a piecewise parabolic approx. for the subgrid distribution.
    !

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlist,ji_p,ji_m,i_startidx,i_endidx,jk_shift,z_iflx_m, &
!$OMP            z_iflx_p,z_delta_m,z_delta_p,z_a11,z_a12)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )


      z_iflx_m(i_startidx:i_endidx,slev:nlevp1) = 0._wp
      z_iflx_p(i_startidx:i_endidx,slev:nlevp1) = 0._wp

      ! Loop over all lists (nlist will serve as shift index)
      DO nlist = 1, nlist_max

        IF (i_listdim_p(nlist,jb) == 0) CYCLE

        !
        ! loop over all cells in i_indlist_p
        !
        ! integer fluxes for w>0 (weta<0)
        !
!CDIR NODEP
        DO ji_p=1,i_listdim_p(nlist,jb)

          ! get jc and jk index from precomputed list
          jc = i_indlist_p(ji_p,nlist,jb)
          jk = i_levlist_p(ji_p,nlist,jb)

          ! integer shift (depends on the applied list)
          jk_shift = jk + nlist - 1

          ! Integer flux (division by p_dtime is done at the end)
          z_iflx_p(jc,jk) = z_iflx_p(jc,jk) - coeff_grid * p_cc(jc,jk_shift,jb) &
            &              * p_cellmass_now(jc,jk_shift,jb)

        ENDDO  ! loop over cells in i_indlist_p


        IF (i_listdim_m(nlist,jb) == 0) CYCLE

        !
        ! loop over all cells in i_indlist_m
        !
        ! integer fluxes for w<0 (weta>0)
        !
!CDIR NODEP
        DO ji_m=1,i_listdim_m(nlist,jb)

          ! get jc and jk index from precomputed list
          jc = i_indlist_m(ji_m,nlist,jb)
          jk = i_levlist_m(ji_m,nlist,jb)

          ! integer shift (depends on the applied list)
          jk_shift = jk - nlist

          ! Integer flux (division by p_dtime is done at the end)
          z_iflx_m(jc,jk) = z_iflx_m(jc,jk) + coeff_grid*p_cc(jc,jk_shift,jb) &
            &              * p_cellmass_now(jc,jk_shift,jb)

        ENDDO  ! loop over cells in i_indlist_m

      ENDDO  ! loop over index lists



      ! compute fractional flux as well as total flux (sum of integer- and 
      ! fractional flux)
      DO jk = slevp1, nlev

        DO jc = i_startidx, i_endidx

          ! Note that currently coeff_grid=-1 for NH-version
          IF (coeff_grid * p_mflx_contra_v(jc,jk,jb) < 0._wp ) THEN

            ! fractional flux
            ! if w > 0 , weta < 0 (physical upwelling)
            ! note that the second coeff_grid factor in front of z_delta_p 
            ! is obsolete due to a compensating (-) sign emerging from the 
            ! computation of z_delta_p.
            z_delta_p = z_face_up(jc,jk_int_p(jc,jk,jb),jb)               &
              &         - z_face_low(jc,jk_int_p(jc,jk,jb),jb)
            z_a12     = p_cc(jc,jk_int_p(jc,jk,jb),jb)                     &
              &         - 0.5_wp * (z_face_low(jc,jk_int_p(jc,jk,jb),jb)   &
              &         + z_face_up(jc,jk_int_p(jc,jk,jb),jb))


            ! fractional high order flux   
            z_flx_frac_high(jc,jk,jb) = ( - coeff_grid                              &
              &         * p_cellmass_now(jc,jk_int_p(jc,jk,jb),jb)                  &
              &         * z_cflfrac_p(jc,jk,jb) *( p_cc(jc,jk_int_p(jc,jk,jb),jb)   &
              &         + (0.5_wp * z_delta_p * (1._wp - z_cflfrac_p(jc,jk,jb)))    &
              &         - z_a12*(1._wp - 3._wp*z_cflfrac_p(jc,jk,jb)                &
              &         + 2._wp*z_cflfrac_p(jc,jk,jb)**2) ) )                       &
              &         / p_dtime


            ! integer flux
            z_flx_int(jc,jk,jb)= z_iflx_p(jc,jk)/p_dtime


!DR            ! fractional low order flux
!DR            z_flx_frac_low(jc,jk,jb) = ( - coeff_grid                            &
!DR              &         * p_cellmass_now(jc,jk_int(jc,jk,jb),jb)                 &
!DR              &         * z_cflfrac_p(jc,jk,jb) * p_cc(jc,jk_int(jc,jk,jb),jb) ) &
!DR              &         / p_dtime

          ELSE

            ! fractional flux
            ! if w < 0 , weta > 0 (physical downwelling)
            ! note that the second coeff_grid factor in front of z_delta_p 
            ! is obsolete due to a compensating (-) sign emerging from the 
            ! computation of z_delta_p.
            z_delta_m = z_face_up(jc,jk_int_m(jc,jk,jb),jb)               &
              &         - z_face_low(jc,jk_int_m(jc,jk,jb),jb)
            z_a11     = p_cc(jc,jk_int_m(jc,jk,jb),jb)                     &
              &         - 0.5_wp * (z_face_low(jc,jk_int_m(jc,jk,jb),jb)   &
              &         + z_face_up(jc,jk_int_m(jc,jk,jb),jb))


            ! fractional high order flux           
            z_flx_frac_high(jc,jk,jb)= ( coeff_grid                                 &
              &         * p_cellmass_now(jc,jk_int_m(jc,jk,jb),jb)                  &
              &         * z_cflfrac_m(jc,jk,jb) * ( p_cc(jc,jk_int_m(jc,jk,jb),jb)  &
              &         - (0.5_wp * z_delta_m * (1._wp - z_cflfrac_m(jc,jk,jb)))    &
              &         - z_a11*(1._wp - 3._wp*z_cflfrac_m(jc,jk,jb)                &
              &         + 2._wp*z_cflfrac_m(jc,jk,jb)**2) ) )                       &
              &         / p_dtime


            ! integer flux
            z_flx_int(jc,jk,jb)= z_iflx_m(jc,jk)/p_dtime

!DR            ! fractional low order flux
!DR            z_flx_frac_low(jc,jk,jb)= ( coeff_grid                                 &
!DR              &         * p_cellmass_now(jc,jk_int_m(jc,jk,jb),jb)                 &
!DR              &         * z_cflfrac_m(jc,jk,jb) * p_cc(jc,jk_int_m(jc,jk,jb),jb) ) &
!DR              &         / p_dtime


          ENDIF

          ! full flux (integer- plus fractional flux)
          p_upflux(jc,jk,jb) = z_flx_int(jc,jk,jb) + z_flx_frac_high(jc,jk,jb)

        ENDDO ! end loop over cells

      ENDDO ! end loop over vertical levels


      !
      ! set upper and lower boundary condition
      !
      CALL set_bc_vadv(p_upflux(:,slev+1,jb),            &! in
        &              p_mflx_contra_v(:,slev+1,jb),     &! in
        &              p_mflx_contra_v(:,slev  ,jb),     &! in
        &              p_iubc_adv, i_startidx, i_endidx, &! in
        &              zparent_topflx(:,jb),             &! in
        &              p_upflux(:,slev,jb),              &! out
        &              p_upflux(:,nlevp1,jb)             )! out

    ENDDO ! end loop over blocks
!$OMP END DO


    ! If desired, get edge value of advected quantity 
    IF ( l_out_edgeval ) THEN
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend )

        DO jk = slevp1, nlev

          DO jc = i_startidx, i_endidx
            p_upflux(jc,jk,jb) = p_upflux(jc,jk,jb)/p_mflx_contra_v(jc,jk,jb)
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO
    ENDIF
!$OMP END PARALLEL


    !
    ! 6. If desired, apply a flux limiter to limit computed fluxes.
    !    These flux limiters are based on work by Zalesak (1979)
    !
    IF ( iequations /= 3 ) THEN
      IF (p_itype_vlimit == ifluxl_vpd) THEN
        ! positive-definite (pd) limiter
        CALL vflx_limiter_pd_ha( p_patch, p_dtime, p_cc, p_upflux,  & !in,inout
          &                   opt_rlstart=i_rlstart,                & !in
          &                   opt_rlend=i_rlend, opt_slev=slev      ) !in
      ENDIF
    ELSE
      IF (p_itype_vlimit == ifluxl_vpd) THEN
        ! positive-definite (pd) limiter
        CALL vflx_limiter_pd( p_patch, p_dtime, p_cc, p_upflux,     & !in,inout
          &                   opt_rlstart=i_rlstart,                & !in
          &                   opt_rlend=i_rlend, opt_slev=slev      ) !in
      ENDIF
    ENDIF


    IF ( ld_cleanup ) THEN
      ! deallocate temporary arrays 
      DEALLOCATE( i_indlist_p, i_indlist_m, i_levlist_p,           &
      &           i_levlist_m, i_listdim_p, i_listdim_m, jk_int_p, &
      &           jk_int_m, z_cflfrac_p, z_cflfrac_m,STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                     &
          &  'deallocation for i_indlist_p, i_indlist_m, i_levlist_p, ' //  &
          &  'i_levlist_m, i_listdim_p, i_listdim_m, jk_int_p, '        //  &
          &  'jk_int_m, z_cflfrac_p, z_cflfrac_m failed '  )
      ENDIF

    END IF


  END SUBROUTINE upwind_vflux_ppm_cfl



  !-------------------------------------------------------------------------
  !>
  !! Set upper and lower boundary condition for vertical transport
  !!
  !! Set upper and lower boundary condition for vertical transport.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-04-12)
  !!
  !
  SUBROUTINE set_bc_vadv(upflx_top_p1, mflx_top_p1, mflx_top, iubc_adv, &
    &                    i_start, i_end, parent_topflx, upflx_top,      &
    &                    upflx_bottom )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_vflux: set_ubc_adv'

    REAL(wp), INTENT(IN)     :: & !< computed tracer flux at second half level
      &  upflx_top_p1(nproma)
    REAL(wp), INTENT(IN)     :: & !< mass flux at second half level
      &  mflx_top_p1(nproma)
    REAL(wp), INTENT(IN)     :: & !< mass flux at upper boundary
      &  mflx_top(nproma)
    INTEGER, INTENT(IN)      :: & !< selects upper boundary condition
      &  iubc_adv
    INTEGER, INTENT(IN)      :: & !< start and end index
      & i_start, i_end
    REAL(wp), INTENT(IN)     :: & !< tracer flux at upper boundary, 
      &  parent_topflx(nproma)    !< interpolated from parent grid
    REAL(wp), INTENT(OUT)    :: & !< upper boundary condition
      &  upflx_top(nproma)
    REAL(wp), INTENT(OUT)    :: & !< lower boundary condition
      &  upflx_bottom(nproma)


    !-------------------------------------------------------------------------

    ! 
    ! flux at top boundary
    ! 
    SELECT CASE (iubc_adv)
      CASE ( ino_flx )     ! no flux
        upflx_top(i_start:i_end) = 0._wp
 
      CASE ( izero_grad )  ! zero gradient
        upflx_top(i_start:i_end) = upflx_top_p1(i_start:i_end)      &
            &           * mflx_top(i_start:i_end)                   &
            &           / ( mflx_top_p1(i_start:i_end)              &
            &           + SIGN(dbl_eps, mflx_top_p1(i_start:i_end)))

      CASE ( iparent_flx ) ! interpolated flux from parent grid
        upflx_top(i_start:i_end) = parent_topflx(i_start:i_end)
    END SELECT

    !
    ! flux at bottom boundary
    !
    upflx_bottom(i_start:i_end) = 0._wp

  END SUBROUTINE set_bc_vadv



END MODULE mo_advection_vflux

