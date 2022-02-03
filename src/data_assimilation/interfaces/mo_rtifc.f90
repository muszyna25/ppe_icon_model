!
!+ Interface for the RTTOV library (version 10 and later)
!
MODULE mo_rtifc
!
! Description:
!   This module and the associated mo_rtifc_*.f90 modules contain subroutines 
!   to work with RTTOV. The purpose of this module is to make work with 
!   different RTTOV versions more user friendly. The interfaces and options of 
!   this module should not change from one RTTOV version to the next - in 
!   contrast to the original RTTOV code.
!   The following is provided:
!    - routines to modify/check/print RTTOV configuration/options
!    - an initialization routine which initializes the RTTOV modules
!      and reads the required instrument specific coefficients.
!      If _RTIFC_DISTRIBCOEF is set during compilation, this routine has the 
!      option to read the coefficients on one PE and distribute them to the
!      others.
!    - Routines to fill the RTTOV profiles structure.
!    - Routines to call rttov in forward and k mode
!    - Some auxiliary routines
!   This module (mo_rtifc.f90) is just a wrapper module, that uses routines, 
!   variables and parameters from associated modules (depending on the Macros that
!   were set during compilation).
!   Basic routines and stuff, that should not change from one RTTOV version to
!   another (e.g. DWD specific stuff) is located in mo_rtifc_base.f90.
!   All stuff that depends somehow on the RTTOV version is located in version
!   specific modules mo_rtifc_${rttov_version}.f90 (rttov_version=nort,10,12,13...).
!
! Current Code Owner: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    fax:   +49 69 8062 3721
!    email:  robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_8         2009/12/09 Marc Schwaerz
!  interface for the rttov (currently version 9) library
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Detlef Pingel
!  update to current version
! V1_19        2012-04-16 Andreas Messer
!  RTTOV10 interface, Optimizations
! V1_20        2012-06-18 Andreas Messer
!  extend to RTTOV10; optimize
! V1_22        2013-02-13 Andreas Messer
!  Merged some changes from Alexandre Lanciani into rtifc.
!  Implementation of vectorized K-mode (Robin Faulwetter).
!  add namelist parameter for setting the rttov levels.
!  move module mo_rtifc.f90 from directory basic to oo-model.
! V1_23        2013-03-26 Andreas Rhodin
!  new parameters: ctp_k, cfraction_k (adjoint cloud top height and fraction)
! V1_26        2013/06/27 Robin Faulwetter
!  Introduced a check on the influence of the surface onto radiances.
!  Introduced USE_MWSURF bit. Corrected the usage of other USE_* bits.
! V1_27        2013-11-08 Robin Faulwetter
!  Implemented FOV-dependent obserrors for TOVS.
!  Fixed undesired behaviour of thinning
! V1_28        2014/02/26 Andreas Rhodin
!  pass cloud-top-height and cloud-fraction to RTTOV routines by process_tovs_mult
! V1_31        2014-08-21 Robin Faulwetter
!  Unify mo_rtifc with COSMO.
!  New, much faster write_rttov_prof routine.
! V1_35        2014-11-07 Andreas Rhodin
!  correct printout in case of error
! V1_42        2015-06-08 Robin Faulwetter
!  Unified modules for radiance processing with COSMO
! V1_43        2015-08-19 Robin Faulwetter
!  Added features required for high peaking radiances over land/clouds
! V1_47        2016-06-06 Robin Faulwetter
!  Many improvements for radiances.
! V1_48        2016-10-06 Robin Faulwetter
!  Implemented RTTOV12
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Marc Schwaerz  EUMETSAT/DWD 2008-08-20 Initial release
! Detlef Pingel           DWD 2009-06-18 Provision of gradient matrix interface
! Marc Schwaerz  EUMETSAT/DWD 2009-07-30 Provision of tangent linear interface
!                                        Implementation of mpi-distributed reading of coeff files
!                                        Multiple profile processing enabled
! Marc Schwaerz  EUMETSAT/DWD 2009-09-30 Possibility to turn off rttov support via a single
!                                        pre-processor directive called NO_RTTOV
! Robin Faulwetter        DWD 2020-05    Split into multiple modules and major cleanup/
!                                        reorganisation
!=======================================================================

#include "mo_rtifc_macros.incf"

  
  !-------------
  ! Modules used
  !-------------

  use mo_rtifc_base

#if (_RTTOV_VERSION == 13) 
  use mo_rtifc_13,       only: rtifc_vers,             &! RTTOV version number this interface is compiled for
                               rtifc_version,          &! version string
                               rtifc_init,             &! Initialize RTTOV, i.e. read coeffs
                               rtifc_coef_index,       &! Return index of coeffs for given satid/instr
                               rtifc_cleanup,          &! Cleanup RTTOV and mo_rtifc structures
                               rtifc_get_preslev,      &! Get RTTOV levels
                               rtifc_fill_input,       &! Fill rttov profile structure
                               rtifc_direct,           &! rttov_direct call
                               rtifc_k,                &! rttov_k call
                               rtifc_dealloc_profiles, &! deallocate rttov profile structure
                               rtifc_print_profiles,   &! print rttov profile structure
                               rtifc_set_opts_sub,     &! set RTTOV options
                               rttov_options,          &! Options for RTTOV (type definition)
                               rttov_opts_def,         &! Default options for RTTOV
                               gas_unit_specconc,      &! specific concentration (kg/kg over wet air)
                               gas_unit_ppmv,          &! ppmv over wet air
                               gas_unit_ppmvdry         ! ppmv over dry air
#if defined(_RTTOV_ATLAS)
  use mo_rtifc_13,       only: rtifc_init_mw_atlas,    &!
                               rtifc_emis_atlas,       &! get emissivity from atlas
                               rtifc_emis_retrieve,    &! retrieve emissivity (by Karbou-method)
                               rtifc_init_brdf_atlas,  &!
                               rtifc_brdf_atlas
#endif

#elif (_RTTOV_VERSION == 12) 
  use mo_rtifc_12,       only: rtifc_vers,             &! RTTOV version number this interface is compiled for
                               rtifc_version,          &! version string
                               rtifc_init,             &! Initialize RTTOV, i.e. read coeffs
                               rtifc_coef_index,       &! Return index of coeffs for given satid/instr
                               rtifc_cleanup,          &! Cleanup RTTOV and mo_rtifc structures
                               rtifc_get_preslev,      &! Get RTTOV levels
                               rtifc_fill_input,       &! Fill rttov profile structure
                               rtifc_direct,           &! rttov_direct call
                               rtifc_k,                &! rttov_k call
                               rtifc_dealloc_profiles, &! deallocate rttov profile structure
                               rtifc_print_profiles,   &! print rttov profile structure
                               rtifc_set_opts_sub,     &! set RTTOV options
                               rttov_options,          &! Options for RTTOV (type definition)
                               rttov_opts_def,         &! Default options for RTTOV
                               gas_unit_specconc,      &! specific concentration (kg/kg over wet air)
                               gas_unit_ppmv,          &! ppmv over wet air
                               gas_unit_ppmvdry         ! ppmv over dry air
#if defined(_RTTOV_ATLAS)
  use mo_rtifc_12,       only: rtifc_init_mw_atlas,    &!
                               rtifc_emis_atlas,       &! get emissivity from atlas
                               rtifc_emis_retrieve,    &! retrieve emissivity (by Karbou-method)
                               rtifc_init_brdf_atlas,  &!
                               rtifc_brdf_atlas
#endif

#elif (_RTTOV_VERSION == 10)
  use mo_rtifc_10,       only: rtifc_vers,             &! RTTOV version this interface is compiled for
                               rtifc_version,          &! version string
                               rtifc_init,             &! Initialize RTTOV, i.e. read coeffs
                               rtifc_coef_index,       &! Return index of coeffs for given satid/instr
                               rtifc_cleanup,          &! Cleanup RTTOV and mo_rtifc structures
                               rtifc_get_preslev,      &! Get RTTOV levels
                               rtifc_fill_input,       &! Fill rttov profile structure
                               rtifc_direct,           &! rttov_direct call
                               rtifc_k,                &! rttov_k call
                               rtifc_dealloc_profiles, &! deallocate rttov profile structure
                               rtifc_print_profiles,   &! print rttov profile structure
                               rtifc_set_opts_sub,     &! set RTTOV options
                               rttov_options,          &! Options for RTTOV (type definition)
                               rttov_opts_def           ! Default options for RTTOV
#elif (_RTTOV_VERSION <= 0)
  use mo_rtifc_nort,     only: rtifc_vers,             &! RTTOV version this interface is compiled for
                               rtifc_version,          &! version string
                               rtifc_init,             &! Initialize RTTOV, i.e. read coeffs
                               rtifc_coef_index,       &! Return index of coeffs for given satid/instr
                               rtifc_cleanup,          &! Cleanup RTTOV and mo_rtifc structures
                               rtifc_get_preslev,      &! Get RTTOV levels
                               rtifc_fill_input,       &! Fill rttov profile structure
                               rtifc_direct,           &! rttov_direct call
                               rtifc_k,                &! rttov_k call
                               rtifc_dealloc_profiles, &! deallocate rttov profile structure
                               rtifc_print_profiles,   &! print rttov profile structure
                               rttov_options,          &! Options for RTTOV (type definition)
                               rttov_opts_def           ! Default options for RTTOV
#endif  


  implicit none

  !-------------
  ! Public stuff
  !-------------

  private

  ! RTTOV version
  public :: rtifc_vers           ! the rttov version this interface was compiled for

  ! subroutines
  public :: rtifc_check_config     ! Check RTTOV version and level number, set nlevs_top
  public :: rtifc_version          ! RTTOV and interface version string
  public :: rtifc_set_opts         ! set RTTOV options
  public :: rtifc_init             ! Initialise RTTOV modules, read coeffs
  public :: rtifc_coef_index       ! Returns index of coeffs for given satid/instr
  public :: rtifc_cleanup          ! frees memory allocated by rtifc_init
  public :: rtifc_fill_input       ! fills the profile-dependent part for RTTOV
  public :: rtifc_get_preslev      ! Get pressure levels
  public :: rtifc_direct           ! calls RTTOV direct routine
  public :: rtifc_k                ! calls RTTOV K routine
  public :: rtifc_dealloc_profiles ! deallocate the profiles
  public :: rtifc_print_profiles   ! print profiles
  public :: rtifc_errmsg           ! gives error message corresponding to exit status
#if defined(_RTTOV_ATLAS)
  ! Emissivity atlases
  public :: rtifc_init_mw_atlas
  public :: rtifc_init_brdf_atlas
  public :: rtifc_emis_atlas
  public :: rtifc_emis_retrieve
  public :: rtifc_brdf_atlas
#endif

  ! Reading/distribution of coeffs
  public :: read1pe

  ! error codes/messages
  public :: NO_ERROR             ! everything was ok.
  public :: WARN_RTTOV           ! warning

  ! options
  public :: rttov_options        ! type definition
  public :: rttov_opts_def       ! type(rttov_options) holding RTTOV default options

  ! output flags
  public :: OUT_ASB
  public :: OUT_CSB
  public :: OUT_ASR
  public :: OUT_CSR
  public :: OUT_VIS

  ! RTTOV levels above user levels
  public :: nlevs_top

  ! default profile values
  public :: default_wfetch     
  public :: default_fastem
  public :: default_watertype  
  public :: default_salinity   
  public :: default_o3_surf    
  public :: default_satazim    
  public :: default_sunzenangle
  public :: default_sunazangle 
  public :: default_ctp        
  public :: default_cfraction
#if (_RTTOV_VERSION >= 13)
  public :: default_ice_scheme
#else
  public :: default_idg
#endif
  public :: default_clw_scheme
  public :: default_gas_units  

  ! hard limits on profile variables
  public :: qmin_ifc
  public :: qmax_ifc
  public :: tmin_ifc
  public :: tmax_ifc

  ! RTTOV "constants"
  public :: min_od
#if (_RTTOV_VERSION >= 12)
  public :: gas_unit_specconc
  public :: gas_unit_ppmv
  public :: gas_unit_ppmvdry
#endif

  ! Regularization limits
  public :: chk_reg_lims
  public :: mask_lims_t
  public :: mask_lims_q 

  ! god (generalized optical depth) parameters
  public :: god_par_file
  public :: wr_god
  public :: out_path
  public :: chk_god
  public :: god_thresh

contains


  subroutine rtifc_set_opts(rttov_opts,         &!
                            init,               &!
                            addinterp,          &!
                            addrefrac,          &!
                            addclouds,          &!
                            addaerosl,          &!
                            addsolar,           &!
                            addpc,              &!
                            apply_reg_lims,     &!
                            verbose_reg_lims,   &!
                            crop_k_reg_lims,    &!
                            switchrad,          &!
                            conv_overc,         &!
                            fix_hgpl,           &!
                            fastem_version,     &!
                            ir_sea_emis_model,  &!
                            use_t2m_opdep,      &!
                            use_q2m,            &!
                            do_lambertian,      &!
                            cloud_overlap,      &!
                            do_checkinput,      &!
                            ozone_data,         &!
                            co2_data,           &!
                            n2o_data,           &!
                            co_data,            &!
                            ch4_data,           &!
                            so2_data,           &!
                            clw_data,           &!
                            dom_rayleigh,       &!
                            dom_nstreams,       &!
                            ir_scatt_model,     &!
                            vis_scatt_model,    &!
                            clip_gas_opdep      &!
                           )
    type(rttov_options), intent(inout), optional :: rttov_opts
    logical,             intent(in),    optional :: init
    logical,             intent(in),    optional :: addinterp
    logical,             intent(in),    optional :: addrefrac
    logical,             intent(in),    optional :: addclouds
    logical,             intent(in),    optional :: addaerosl
    logical,             intent(in),    optional :: addsolar
    logical,             intent(in),    optional :: addpc
    logical,             intent(in),    optional :: apply_reg_lims
    logical,             intent(in),    optional :: verbose_reg_lims
    logical,             intent(in),    optional :: crop_k_reg_lims
    logical,             intent(in),    optional :: switchrad
    logical,             intent(in),    optional :: conv_overc
    integer,             intent(in),    optional :: fix_hgpl
    integer,             intent(in),    optional :: fastem_version
    integer,             intent(in),    optional :: ir_sea_emis_model
    logical,             intent(in),    optional :: use_t2m_opdep
    logical,             intent(in),    optional :: use_q2m
    logical,             intent(in),    optional :: do_lambertian
    integer,             intent(in),    optional :: cloud_overlap
    logical,             intent(in),    optional :: do_checkinput
    logical,             intent(in),    optional :: ozone_data
    logical,             intent(in),    optional :: co2_data  
    logical,             intent(in),    optional :: n2o_data  
    logical,             intent(in),    optional :: co_data   
    logical,             intent(in),    optional :: ch4_data  
    logical,             intent(in),    optional :: so2_data  
    logical,             intent(in),    optional :: clw_data
    logical,             intent(in),    optional :: dom_rayleigh
    integer,             intent(in),    optional :: dom_nstreams
    integer,             intent(in),    optional :: ir_scatt_model
    integer,             intent(in),    optional :: vis_scatt_model
    logical,             intent(in),    optional :: clip_gas_opdep

#if (_RTTOV_VERSION == 13)
    call rtifc_set_opts_sub(rttov_opts        = rttov_opts,         &!
                            init              = init,               &!
                            addinterp         = addinterp,          &!
                            addrefrac         = addrefrac,          &!
                            addclouds         = addclouds,          &!
                            addaerosl         = addaerosl,          &!
                            addsolar          = addsolar,           &!
                            addpc             = addpc,              &!
                            apply_reg_lims    = apply_reg_lims,     &!
                            verbose_reg_lims  = verbose_reg_lims,   &!
                            crop_k_reg_lims   = crop_k_reg_lims,    &!
                            switchrad         = switchrad,          &!
                            conv_overc        = conv_overc,         &!
                            fix_hgpl          = fix_hgpl,           &!
                            fastem_version    = fastem_version,     &!
                            ir_sea_emis_model = ir_sea_emis_model,  &!
                            use_t2m_opdep     = use_t2m_opdep,      &!
                            use_q2m           = use_q2m,            &!
                            do_lambertian     = do_lambertian,      &!
                            cloud_overlap     = cloud_overlap,      &!
                            do_checkinput     = do_checkinput,      &!
                            ozone_data        = ozone_data,         &!
                            co2_data          = co2_data,           &!
                            n2o_data          = n2o_data,           &!
                            co_data           = co_data,            &!
                            ch4_data          = ch4_data,           &!
                            so2_data          = so2_data,           &!
                            clw_data          = clw_data,           &!
                            dom_rayleigh      = dom_rayleigh,       &!
                            dom_nstreams      = dom_nstreams,       &!
                            ir_scatt_model    = ir_scatt_model,     &!
                            vis_scatt_model   = vis_scatt_model,    &!
                            clip_gas_opdep    = clip_gas_opdep      &!
                           )
#elif (_RTTOV_VERSION == 12)
    call rtifc_set_opts_sub(rttov_opts        = rttov_opts,         &!
                            init              = init,               &!
                            addinterp         = addinterp,          &!
                            addrefrac         = addrefrac,          &!
                            addclouds         = addclouds,          &!
                            addaerosl         = addaerosl,          &!
                            addsolar          = addsolar,           &!
                            addpc             = addpc,              &!
                            apply_reg_lims    = apply_reg_lims,     &!
                            verbose_reg_lims  = verbose_reg_lims,   &!
                            crop_k_reg_lims   = crop_k_reg_lims,    &!
                            switchrad         = switchrad,          &!
                            conv_overc        = conv_overc,         &!
                            fix_hgpl          = fix_hgpl,           &!
                            fastem_version    = fastem_version,     &!
                            ir_sea_emis_model = ir_sea_emis_model,  &!
                            use_q2m           = use_q2m,            &!
                            do_lambertian     = do_lambertian,      &!
                            ozone_data        = ozone_data,         &!
                            co2_data          = co2_data,           &!
                            n2o_data          = n2o_data,           &!
                            co_data           = co_data,            &!
                            ch4_data          = ch4_data,           &!
                            so2_data          = so2_data,           &!
                            clw_data          = clw_data            &!
                           )
#elif (_RTTOV_VERSION == 10)
    call rtifc_set_opts_sub(rttov_opts        = rttov_opts,         &!
                            init              = init,               &!
                            addinterp         = addinterp,          &!
                            addrefrac         = addrefrac,          &!
                            addclouds         = addclouds,          &!
                            addaerosl         = addaerosl,          &!
                            addsolar          = addsolar,           &!
                            addpc             = addpc,              &!
                            apply_reg_lims    = apply_reg_lims,     &!
                            verbose_reg_lims  = verbose_reg_lims,   &!
                            switchrad         = switchrad,          &!
                            fastem_version    = fastem_version,     &!
                            use_q2m           = use_q2m             &!
                           )
#endif
    
  end subroutine rtifc_set_opts
                            

  function rtifc_errmsg(code, v) result(msg)
    character(len=120)            :: msg
    integer, intent(in)           :: code
    logical, intent(in), optional :: v

    logical :: v_

    msg = trim(errmsg(code))

    if (present(v)) then
      v_ = v
    else
      v_ = .false.
    end if
    if (v_) then
      msg = rtifc_version()//' '//trim(msg)
    end if

  end function rtifc_errmsg

end module mo_rtifc
