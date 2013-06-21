!>
!! Namelist for configuration of grid refinement algorithms.
!!
!! @par Revision History
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
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio 
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile

  USE mo_gridref_config,      ONLY:   &
    &                            config_rbf_vec_kern_grf_e => rbf_vec_kern_grf_e,& 
    &                            config_rbf_scale_grf_e    => rbf_scale_grf_e,&
    &                            config_grf_velfbk         => grf_velfbk,&
    &                            config_grf_scalfbk        => grf_scalfbk,&
    &                            config_grf_tracfbk        => grf_tracfbk,&
    &                            config_grf_idw_exp_e12    => grf_idw_exp_e12,&
    &                            config_grf_idw_exp_e34    => grf_idw_exp_e34,&
    &                            config_grf_intmethod_c    => grf_intmethod_c,&
    &                            config_grf_intmethod_e    => grf_intmethod_e,&
    &                            config_grf_intmethod_ct   => grf_intmethod_ct,&
    &                            config_l_mass_consvcorr   => l_mass_consvcorr,&
    &                            config_l_density_nudging  => l_density_nudging,&
    &                            config_denom_diffu_v      => denom_diffu_v,&
    &                            config_denom_diffu_t      => denom_diffu_t
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings, log_nml_settings 


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_gridref_namelist

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !---------------------
  ! namelist variables
  !---------------------

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

  LOGICAL  :: l_mass_consvcorr  ! .true.: apply mass conservation correction
  LOGICAL  :: l_density_nudging ! .true.: apply density nudging near lateral nest boundaries if feedback is turned on
                                ! (in case of one-way nesting, all prognostic variables are nudged irrespective of this switch)

  ! Exponents for IDW interpolation in idw_compute_coeff_grf
  REAL(wp) :: grf_idw_exp_e12, grf_idw_exp_e34

  ! Denominators of normalized diffusion coefficients for boundary diffusion
  REAL(wp) :: denom_diffu_v, denom_diffu_t

  NAMELIST/gridref_nml/  rbf_vec_kern_grf_e, rbf_scale_grf_e,             &
    &                    grf_velfbk, grf_scalfbk, grf_tracfbk,            &
    &                    grf_idw_exp_e12, grf_idw_exp_e34,                &
    &                    grf_intmethod_c, grf_intmethod_e,                &
    &                    grf_intmethod_ct, denom_diffu_v, denom_diffu_t,  &
    &                    l_mass_consvcorr, l_density_nudging

CONTAINS
  !-------------------------------------------------------------------------
  !>
  !! This subroutine 
  !! - reads the Namelist for local grid refinement 
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-07-06)
  !!
  SUBROUTINE read_gridref_namelist( filename )

    CHARACTER(LEN=*),INTENT(IN) :: filename
    INTEGER :: istat, funit
    !0!CHARACTER(len=*),PARAMETER :: routine = 'mo_gridref_nml:read_gridref_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

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

    ! Initialize namelist fields for scaling factors (dimension 1:max_dom); used part only
    rbf_scale_grf_e(1:max_dom) = 0.5_wp  ! default setting for vector grf interpolation

    ! Denominator for temperature boundary diffusion
    denom_diffu_t = 135._wp

    ! Denominator for velocity boundary diffusion
    denom_diffu_v = 200._wp

    ! Mass conservation correction turned on by default
    l_mass_consvcorr = .TRUE. 

    ! Density nudging near nest boundaries turned on by default
    l_density_nudging = .TRUE.

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('gridref_nml')
      READ(funit,NML=gridref_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('gridref_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      WRITE(temp_defaults(), gridref_nml)                     ! write defaults to temporary text file
      READ (nnml, gridref_nml, iostat=istat)                  ! overwrite default settings
      WRITE(temp_settings(), gridref_nml)                     ! write settings to temporary text file
      CALL log_nml_settings("nml.log")
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
      config_rbf_vec_kern_grf_e = rbf_vec_kern_grf_e
      config_rbf_scale_grf_e = rbf_scale_grf_e
      config_grf_velfbk = grf_velfbk
      config_grf_scalfbk = grf_scalfbk
      config_grf_tracfbk = grf_tracfbk
      config_grf_idw_exp_e12 = grf_idw_exp_e12
      config_grf_idw_exp_e34 = grf_idw_exp_e34
      config_grf_intmethod_c = grf_intmethod_c
      config_grf_intmethod_e = grf_intmethod_e
      config_grf_intmethod_ct = grf_intmethod_ct
      config_denom_diffu_v = denom_diffu_v
      config_denom_diffu_t = denom_diffu_t
      config_l_mass_consvcorr = l_mass_consvcorr
      config_l_density_nudging = l_density_nudging

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=gridref_nml)                    
      CALL store_and_close_namelist(funit, 'gridref_nml') 
    ENDIF
    ! 6. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=gridref_nml)

  END SUBROUTINE read_gridref_namelist

END MODULE mo_gridref_nml
