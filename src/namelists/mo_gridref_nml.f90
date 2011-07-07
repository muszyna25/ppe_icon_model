!>
!! Namelist for the configuration of the grid refinement.
!!
!! This subroutine is called by control_model
!!
!! @par Revision History
!! Revision history in mo_grf_intp_stat (r4300)
!! Modification by Constantin Junk, MPI-M (2011-05-05)
!! - some updates on the structure
!! - moved gridref namelist variables from mo_grf_intp_data_strc to
!!   mo_gridref_nml
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
MODULE mo_gridref_nml

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_exception,           ONLY: message, finish
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH
  USE mo_master_nml,          ONLY: lrestart
  USE mo_model_domain_import, ONLY: n_dom
  USE mo_namelist,            ONLY: position_nml, POSITIONED
  USE mo_gridref_config,      ONLY: gridref_config
  USE mo_run_nml,             ONLY: lshallow_water
  USE mo_mpi,                 ONLY: p_pe, p_io
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !----------------------------------!
  ! gridref_ctl namelist variables   !
  !----------------------------------!

  INTEGER  :: rbf_vec_kern_grf_e ! rbf kernel for vector interpolation

  ! scale factors for rbf grid refinement interpolation
  REAL(wp) :: rbf_scale_grf_e(max_dom)

  INTEGER  :: grf_intmethod_c,  &  ! switch for type of grid refinement interpolation
    &         grf_intmethod_ct, &  ! (see below for explanation of options)
    &         grf_intmethod_e

  INTEGER  :: grf_velfbk     ! switch for velocity feedback method
                             ! 1 = averaging over child edges 1 and 2;
                             ! 2 = 2nd-order method using RBF reconstruction to child vertices
  
  INTEGER  :: grf_scalfbk    ! switch for feedback method of scalar dynamical variables
                             ! 1 = area-weighted averaging
                             ! 2 = bilinear interpolation

  INTEGER  :: grf_tracfbk    ! switch for feedback method of passive tracer variables
                             ! 1 = area-weighted averaging
                             ! 2 = bilinear interpolation

  ! Exponents for IDW interpolation in idw_compute_coeff_grf
  REAL(wp) :: grf_idw_exp_e12, grf_idw_exp_e34

  ! Denominators of normalized diffusion coefficients for boundary diffusion
  REAL(wp) :: denom_diffu_v, denom_diffu_t

  NAMELIST/gridref_ctl/  rbf_vec_kern_grf_e, rbf_scale_grf_e,             &
    &                    grf_velfbk, grf_scalfbk, grf_tracfbk,            &
    &                    grf_idw_exp_e12, grf_idw_exp_e34,                &
    &                    grf_intmethod_c, grf_intmethod_e,                &
    &                    grf_intmethod_ct, denom_diffu_v, denom_diffu_t

  PUBLIC :: rbf_vec_kern_grf_e, rbf_scale_grf_e,                          &
    &                    grf_velfbk, grf_scalfbk, grf_tracfbk,            &
    &                    grf_idw_exp_e12, grf_idw_exp_e34,                &
    &                    grf_intmethod_c, grf_intmethod_e,                &
    &                    grf_intmethod_ct, denom_diffu_v, denom_diffu_t

  PUBLIC :: gridref_nml_setup, read_gridref_namelist

CONTAINS

!>
!!   Set up the configuration for grid refinement.
!!
!!
!! @par Revision History
!!  Created by Guenther Zaengl, MPI-M (2009-02-09).
!!  Modification by Constantin Junk, MPI-M (2011-05-05)
!!  - renamed setup_gridref to gridref_nml_setup
!!
SUBROUTINE gridref_nml_setup

  INTEGER :: istat, funit

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = 'mo_gridref_nml:setup_gridref'

  !------------------------------------------------------------
  ! set up the default values for gridref_ctl
  !------------------------------------------------------------

  CALL message(TRIM(routine),'default settings')

  ! Switch for interpolation method used for cell-based dynamical
  grf_intmethod_c   = 2         ! 1: copying, 2: gradient-based interpolation
  ! Switch for interpolation method used for tracer variables
  grf_intmethod_ct  = 2         ! 1: copying, 2: gradient-based interpolation
  ! Currently, grf_intmethod_c is used for temperature only; other variables are copied
  grf_intmethod_e   = 4         ! 1: IDW, 2: RBF, 3: IDW/gradient-based, 4: RBF/gradient-based
  ! Switch for velocity feedback method.
  grf_velfbk      = 1         ! 1: average over child edges 1 and 2
                              ! 2: 2nd-order method using RBF reconstruction to child vertices
  ! Switch for feedback method for scalar dynamical variables
  grf_scalfbk     = 2         ! 1: area-weighted averaging
                              ! 2: bilinear interpolation
  ! Switch for feedback method for passive tracer variables
  grf_tracfbk     = 2         ! 1: area-weighted averaging
                              ! 2: bilinear interpolation
  
  ! Exponents for IDW interpolation function
  grf_idw_exp_e12  = 1.2_wp   ! child edges 1 and 2
  grf_idw_exp_e34  = 1.7_wp   ! child edges 3 and 4

  ! RBF kernels for grid refinement interpolation
  rbf_vec_kern_grf_e = 1        ! 1: Gaussian, 2: 1/(1+r**2), 3: inverse multiquadric

  ! zero whole arrays
  rbf_scale_grf_e(:) = 0.0_wp

  ! Initialize namelist fields for scaling factors (dimension 1:n_dom); used part only
  rbf_scale_grf_e(1:n_dom) = 0.5_wp  ! default setting for vector grf interpolation

  ! Denominator for temperature boundary diffusion
  denom_diffu_t = 135._wp

  ! Denominator for velocity boundary diffusion
  IF (lshallow_water) THEN
    denom_diffu_v = 250._wp
  ELSE
    denom_diffu_v = 200._wp
  ENDIF

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values in the previous integration.
    !----------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('gridref_ctl')
      READ(funit,NML=gridref_ctl)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL position_nml ('gridref_ctl', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, gridref_ctl)
    END SELECT

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=gridref_ctl)                                                             
    CALL store_and_close_namelist(funit, 'gridref_ctl')                                      

  ! write the contents of the namelist to an ASCII file
  
  IF(p_pe == p_io) WRITE(nnml_output,nml=gridref_ctl)


END SUBROUTINE gridref_nml_setup


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for NWP physics. 
  !!
  !! This subroutine 
  !! - reads the Namelist for NWP physics
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-06-07)
  !!
  SUBROUTINE read_gridref_namelist
    !
    INTEGER :: istat, funit
    INTEGER :: jg           ! loop index

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_gridref_nml: read_gridref_namelist'

    !-----------------------------------------------------------------------

    !-----------------------!
    ! 1. default settings   !
    !-----------------------!

    ! Switch for interpolation method used for cell-based dynamical
    grf_intmethod_c   = 2     ! 1: copying, 2: gradient-based interpolation
    ! Switch for interpolation method used for tracer variables
    grf_intmethod_ct  = 2     ! 1: copying, 2: gradient-based interpolation
    ! Currently, grf_intmethod_c is used for temperature only; other variables are copied
    grf_intmethod_e   = 4     ! 1: IDW, 2: RBF, 3: IDW/gradient-based, 
                              ! 4: RBF/gradient-based

    ! Switch for velocity feedback method.
    grf_velfbk      = 1       ! 1: average over child edges 1 and 2
                              ! 2: 2nd-order method using RBF reconstruction to child vertices
    ! Switch for feedback method for scalar dynamical variables
    grf_scalfbk     = 2       ! 1: area-weighted averaging
                              ! 2: bilinear interpolation
    ! Switch for feedback method for passive tracer variables
    grf_tracfbk     = 2       ! 1: area-weighted averaging
                              ! 2: bilinear interpolation
  
    ! Exponents for IDW interpolation function
    grf_idw_exp_e12  = 1.2_wp ! child edges 1 and 2
    grf_idw_exp_e34  = 1.7_wp ! child edges 3 and 4

    ! RBF kernels for grid refinement interpolation
    rbf_vec_kern_grf_e = 1    ! 1: Gaussian, 2: 1/(1+r**2), 3: inverse multiquadric

    ! zero whole arrays
    rbf_scale_grf_e(:) = 0.0_wp

    ! Initialize namelist fields for scaling factors (dimension 1:n_dom); used part only
    rbf_scale_grf_e(1:n_dom) = 0.5_wp  ! default setting for vector grf interpolation

    ! Denominator for temperature boundary diffusion
    denom_diffu_t = 135._wp

    ! Denominator for velocity boundary diffusion
    IF (lshallow_water) THEN
      denom_diffu_v = 250._wp
    ELSE
      denom_diffu_v = 200._wp
    ENDIF


    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('gridref_ctl')
      READ(funit,NML=gridref_ctl)
      CALL close_tmpfile(funit)
    END IF


    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL position_nml ('gridref_ctl', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, gridref_ctl)
    END SELECT


    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
    DO jg = 1,max_dom
      gridref_config(jg)%rbf_vec_kern_grf_e = rbf_vec_kern_grf_e
      gridref_config(jg)%rbf_scale_grf_e = rbf_scale_grf_e
      gridref_config(jg)%grf_velfbk = grf_velfbk
      gridref_config(jg)%grf_scalfbk = grf_scalfbk
      gridref_config(jg)%grf_tracfbk = grf_tracfbk
      gridref_config(jg)%grf_idw_exp_e12 = grf_idw_exp_e12
      gridref_config(jg)%grf_idw_exp_e34 = grf_idw_exp_e34
      gridref_config(jg)%grf_intmethod_c = grf_intmethod_c
      gridref_config(jg)%grf_intmethod_e = grf_intmethod_e
      gridref_config(jg)%grf_intmethod_ct = grf_intmethod_ct
      gridref_config(jg)%denom_diffu_v = denom_diffu_v
      gridref_config(jg)%denom_diffu_t = denom_diffu_t
    ENDDO

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=gridref_ctl)                    
    CALL store_and_close_namelist(funit, 'gridref_ctl') 


    ! 6. write the contents of the namelist to an ASCII file
    !
    IF(p_pe == p_io) WRITE(nnml_output,nml=gridref_ctl)


  END SUBROUTINE read_gridref_namelist

END MODULE mo_gridref_nml
