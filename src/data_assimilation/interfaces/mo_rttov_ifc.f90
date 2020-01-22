!
!+ Interface for the RTTOV library (version 9 and later)
!
MODULE mo_rttov_ifc
!
! Description:
!   This module contains functions and subroutines needed to work with RTTOV.
!   The following is provided:
!    - an initialization routine which initializes the RTTOV modules
!      and reads the instrument specific coefficients which are needed.
!      in addition it fills the default values for some parts
!      of the RTTOV profile structure such as Ozone, etc.
!    - a routine which provides a cleanup functionality for RTTOV
!      to allow a clean shutdown of the model.
!    - routines which fill the variable part of the rttov profile
!      including optional arguments which allow you
!      to change also the default initialization values
!    - routines which call rttov in forward, adjoint, tangent linear,
!      and gradient matrix mode
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
!  Merged some changes from Alexandre Lanciani into rttov_ifc.
!  Implementation of vectorized K-mode (Robin Faulwetter).
!  add namelist parameter for setting the rttov levels.
!  move module mo_rttov_ifc.f90 from directory basic to oo-model.
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
!  Unify mo_rttov_ifc with COSMO.
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
!=======================================================================

! Set Macros for the case, that this is part of ICON (and not DACE)
!GZ: this needs to be revised because ICON is supposed to use the library in any case
! #if defined(__ICON__) && defined(__USE_RTTOV) && defined(__DACE__)
! #define NO_RTTOV
! #endif
! This is needed to allow compiling ICON without RTTOV
#if defined(__ICON__) && !defined(__USE_RTTOV)
#define NO_RTTOV
#endif
#if defined(__ICON__) && !defined(NO_RTTOV)
#define _RTTOV_DO_DISTRIBCOEF
#define RTTOV12 
#define RTTOV_USE_OPENMP
#define RADSHARE
#ifdef __DACE__
#define RTTOV_IFC_USE_MPI_DACE
#else
! TODO use icon mpi routines?
#define RTTOV_IFC_USE_MPIF
#endif
#endif

!define RTTOV12 if neither RTTOV12 or RTTOV10 is defined
#if !defined(NO_RTTOV) && !defined(RTTOV10) && !defined(RTTOV12)
#define RTTOV12
#endif
!check if we shall distribute coefficients via MPI
#ifndef NOMPI
#if defined(_RTTOV_DO_DISTRIBCOEF)
!default is to use internal MPI interface (mo_mpi is not available
!everywhere
#if !defined(RTTOV_IFC_USE_MPI_DACE) && !defined(RTTOV_IFC_USE_MPIF)
!#define RTTOV_IFC_USE_MPIF
#define RTTOV_IFC_USE_MPI_DACE
#endif
#else
! make sure that all defines are disabled
#undef RTTOV_IFC_USE_MPI_DACE
#undef RTTOV_IFC_USE_MPIF
#endif
#else  /* NOMPI */
#undef RTTOV_IFC_USE_MPI_DACE
#undef RTTOV_IFC_USE_MPIF
#undef _RTTOV_DO_DISTRIBCOEF
#endif /* NOMPI */
! OpenMP
#ifdef _OPENMP
#if defined(RTTOV12)
#define RTTOV_USE_OPENMP
#endif
#endif

! TODO: think about an alternative to add_line
! #ifndef RADSHARE
! #define WR call add_line(trim(msg))
! #else
#define WR write(*,*) trim(msg)
! #endif

  !=============
  ! modules used
  !=============
!NOIEXPAND (rttov_clear_rad_var)
#if !defined(RADSHARE) && !defined(NO_RTTOV)
  use mo_rad,             only: t_radv                  ! derived type to store radiance obs.
#endif


  use iso_fortran_env,    only: stdout => output_unit, &!
                                iostat_end              !

#if defined(RTTOV_IFC_USE_MPIF) || defined(RTTOV_IFC_USE_MPI_DACE)
#ifdef  HAVE_MPI_MOD      /* prefer MPI module over mpif.h */
  use mpi
#endif
#endif

#ifndef NO_RTTOV
  !--------------
  ! rttov modules
  !--------------
  use rttov_const,        only: fastem_sp,             &! max. number of fastem surface parameters
                                errorstatus_success,   &! RTTOV errorstatus (0)
                                errorstatus_fatal,     &! RTTOV errorstatus (2)
                                gas_id_mixed,          &!
                                ngases_max,            &!
                                ncldtyp,               &!
                                qmin_rttov => qmin,    &! minimum allowed water vapour
                                qmax_rttov => qmax,    &! maximum allowed water vapour
                                tmin_rttov => tmin,    &! minimum allowed temperature
                                tmax_rttov => tmax      ! maximum allowed temperature

  use rttov_types,        only: rttov_coef,            &! Structure for RTTOV basic coefficients
                                rttov_coef_scatt_ir,   &! Structure for RTTOV IR scattering coefficients
                                rttov_optpar_ir         ! Structure for set of RTTOV IR scattering coefficients

  use parkind1,           only: jpim,                  &! default integer;
                                jpis,                  &! small integer
                                jprb                    ! default double precision (hopefully)

#if !defined(RADSHARE) && !defined(NO_RTTOV)
  use mo_rttov_ifc_tools, only: get_height,            &!
                                god_thresh              !
#endif

#if defined(RTTOV10)
  use rttov_const,        only: errorstatus_warning     ! RTTOV errorstatus (1)
#endif

#if defined(RTTOV12)
#ifndef RADSHARE
  use rttov_math_mod,       only: planck
  use mod_rttov_emis_atlas, only: rttov_emis_atlas_data ! Data type to hold atlas info
  use mod_mwatlas_m2,       only: telsem2_atlas_data    ! Data type to hold TELSEM2 atlas info
  use mod_cnrm_mw_atlas,    only: cnrm_mw_atlas_data    ! Data type to hold CNRM atlas info
  use mo_rttov_ifc_tools,   only: check_god_infl
  use rttov_god,            only: rttov_god2o,           &!
                                  rttov_d_god2o,         &!
                                  rttov_god_init          !
  use rttov_types,          only: rttov_god_par           !

#endif

  use rttov_const,        only: baran_ngauss,          &!
                                gas_unit_specconc,     &! specific concentration (kg/kg over wet air)
                                gas_unit_ppmv,         &! ppmv over wet air
                                gas_unit_ppmvdry,      &! ppmv over dry air
                                mair,                  &!
                                mh2o                    !

  use rttov_types,        only: rttov_phasefn_int,     &! interpolation for phase functions
                                rttov_coef_optpiclb,   &! interpolation factors for frequency
                                rttov_fast_coef,       &! fast coefficients
                                rttov_fast_coef_gas,   &! gas fast coefficients
                                rttov_nlte_coef,       &! internal NLTE coeffs
                                rttov_profile,         &! structure for atmospheric profiles
                                rttov_transmission,    &! Transmissions and optical depths (unitless)
                                rttov_radiance,        &! Radiance and corresponding brightness temperature
                                rttov_radiance2,       &! upwelling and downwelling radiances
                                rttov_emissivity,      &!
                                rttov_coef_pccomp2      ! Structure for RTTOV principal component coefficients

#else
  use rttov_types,        only: rttov_profile      =>  &! structure for atmospheric profiles
                                    profile_type,      &!
                                rttov_transmission =>  &! Transmissions and optical depths (unitless)
                                    transmission_type, &!
                                rttov_radiance     =>  &! Radiance and corresponding brightness temperature
                                    radiance_type       !

#endif

  use rttov_types,        only: rttov_options,         &! Structure for RTTOV options
                                rttov_coefs,           &! Structure for the RTTOV coefficients
                                rttov_coef_pccomp,     &! Structure for RTTOV principal component coefficients
                                rttov_coef_pccomp1,    &! Structure for RTTOV principal component coefficients
                                rttov_chanprof          !

#if !defined(RADSHARE) && !defined(NO_RTTOV)
  use rttov_types,        only: ipr_deb,               &!
                                pe_rt => pe
#endif

  use parkind1,           only: jplm                    ! default logical type

#ifdef RTTOV_IFC_USE_MPI_DACE
  use mo_mpi_dace,        only: p_bcast,               &!
                                dace                    ! MPI group info
#endif


#ifdef RTTOV_USE_OPENMP
#ifndef RADSHARE
  use mo_omp,             only: omp_get_max_threads     !
#endif
#endif

#endif  /* !NO_RTTOV */


  implicit none

 !------------------------------------------------------------------------------
 !================
 ! public entities
 !================
  private
 !-----------------
 ! module variables
 !-----------------
  public :: NO_ERROR             ! everything was ok.
  public :: ERR_ALLOC            ! allocation error
  public :: ERR_DIM              ! mismatch in dimensions
  public :: ERR_RTTOV_SETUP      ! error in rttov setup
  public :: ERR_CLOUD_AERO_MISSM ! aerosol and cloud classes change
                                 ! for different coef_scatt_ir structures
  public :: ERROR_RTTOV_CALL     ! unable to calculate bt for all profiles
  public :: ERROR_RTTOV_CALL_ANY ! unable to calculate bt for some profiles
  public :: WARN_RTTOV_DIR_ALL   ! warning for all profiles
                                 ! nevertheless bt were calculated
  public :: WARN_RTTOV_DIR_ANY   ! warning for some profiles
                                 ! nevertheless bt were calculated
  public :: ERR_RTTOV_MPI        ! mpi comm. error
                                 ! while reading/distributing the coeffs.
  public :: ERR_RTTOV_CL_AER     ! mismatch in rttov setup and profile
                                 ! initialization (cloud/aerosol treatment)
  public :: ERR_NO_RTTOV_LIB     ! mismatch in rttov setup and profile
                                 ! initialization (cloud/aerosol treatment)
  public :: ERR_INVALID_TSKIN    ! invalid t_skin
  public :: rttov_ifc_errMsg     ! holds the error messages
                                 ! for the errors above
  public :: RTTOV_IFC_VERSION    ! the rttov version this interface was
                                 ! compiled for
  public :: OUT_ASB
  public :: OUT_CSB
  public :: OUT_ASR
  public :: OUT_CSR
  public :: default_fastem_version
  public :: default_ir_emis_version
#ifdef RTTOV12
  public :: default_gas_units
  public :: gas_unit_specconc
  public :: gas_unit_ppmv
  public :: gas_unit_ppmvdry
#endif
#if defined(RTTOV12)
  public :: default_do_lambertian
#endif
#ifndef NO_RTTOV
  public :: default_use_q2m
  public :: default_salinity
  public :: fix_hgpl

  public :: app_reg_lims ! Apply regularization limits in RTTOV
  public :: chk_reg_lims ! Check regularization limits in RTTOV
  public :: mask_lims_t  ! mask for T check (app_reg_lims)
  public :: mask_lims_q  ! mask for Q check (app_reg_lims)
  public :: chk_god      ! Check influence of god smoothing
#endif
  !----------
  ! functions
  !----------
  public :: rttov_fill_input       ! fills the profile-dependent part for RTTOV
  public :: rttov_init             ! Initialise RTTOV modules, read coeffs
  public :: rttov_get_preslev      ! Get pressure levels
  public :: rttov_cleanup          ! frees memory allocated by rttov_init
  public :: rttov_direct_ifc       ! calls RTTOV direct routine
  public :: rttov_k_ifc            ! calls RTTOV K routine
  public :: rttov_dealloc_profiles ! deallocate the profiles
#if !defined(RADSHARE) && !defined(NO_RTTOV)
  public :: retrieve_emissivity    ! Emissivity dynamical retrieval
  public :: init_rttov_mw_atlas    ! Initialises emissivity atlases in RTTOV
  public :: get_emis_atlas         ! Fetchs emissivity from atlas
  public :: god_thresh   ! Threshold for chk_god
  public :: god_par_file
#endif
  public :: wr_god
  public :: out_path               ! output path
#ifndef NO_RTTOV
  public :: qmin_ifc
  public :: qmax_ifc
  public :: tmin_ifc
  public :: tmax_ifc
  public :: rttov_print_profiles
#ifdef RTTOV12 
#ifndef RADSHARE
  public :: mw_atlas              ! stores atlas information for use
#endif
#endif

#include "rttov_direct.interface"
#include "rttov_k.interface"
#ifdef RTTOV_USE_OPENMP
#include "rttov_parallel_direct.interface"
#include "rttov_parallel_k.interface"
#endif
#ifdef RTTOV10
#include "rttov_setup.interface"
#endif
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_errorhandling.interface"
#include "rttov_errorreport.interface"
#include "rttov_copy_prof.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_coeffname.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_init_transmission.interface"
#if defined(_RTTOV_DO_DISTRIBCOEF)
#include "rttov_nullify_coef.interface"
#include "rttov_nullify_coef_scatt_ir.interface"
#include "rttov_nullify_coef_pccomp.interface"
#include "rttov_nullify_optpar_ir.interface"
#endif
#include "rttov_read_coefs.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_coefs.interface"
#if defined(RTTOV12)
#include "rttov_convert_profile_units.interface"
#include "rttov_apply_reg_limits.interface"
#include "rttov_hdf_save.interface"
#include "rttov_setup_emis_atlas.interface"
#include "rttov_get_emis.interface"
#include "rttov_deallocate_emis_atlas.interface"
#elif defined(RTTOV10)
#include "rttov_checkinput.interface"
#endif

! mo_mpi does not yet support bcast for complex
! therefore we have to include mpif.h if using mo_mpi
#if defined(RTTOV_IFC_USE_MPIF) || defined(RTTOV_IFC_USE_MPI_DACE)
#ifndef HAVE_MPI_MOD      /* prefer MPI module over mpif.h */
INCLUDE "mpif.h"
#endif
#endif

#if defined (_FTRACE) && !defined (DISABLE_FTRACE_REGION)
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif

#endif

  !================================
  ! Module variables and data types
  !================================
  !--------------------
  ! data type precision
  !--------------------
  ! data type precision kinds 3dVar:
  integer      ,parameter::dp = selected_real_kind(13)! Double Precision
  integer      ,parameter::sp = selected_real_kind(6) ! Single Precision
  ! used data type precision kind
  integer      ,parameter::wp = dp                    ! Working precision
  !-------------
  ! MPI exchange
  !-------------
#ifndef NO_RTTOV
#ifdef _RTTOV_DO_DISTRIBCOEF
  !----------------------------------------
  ! interface p_bcast is defined in mo_mpi,
  ! extend it here
  !----------------------------------------
  interface p_bcast
    module procedure p_bcast_rttov_coefs
    module procedure p_bcast_rttov_coef
    module procedure p_bcast_rttov_coef_scatt_ir
    module procedure p_bcast_rttov_optpar_ir
#if defined(RTTOV12)
    module procedure p_bcast_rttov_coef_optpiclb
    module procedure p_bcast_rttov_phasefn_int
    module procedure p_bcast_rttov_fast_coef
    module procedure p_bcast_rttov_fast_coef_gas
    module procedure p_bcast_rttov_nlte_coef
#ifndef RADSHARE
    module procedure p_bcast_rttov_atlas
    module procedure p_bcast_rttov_telsem
    module procedure p_bcast_rttov_cnrm
    module procedure p_bcast_god_par
#endif
    module procedure p_bcast_rttov_coef_pccomp2
    module procedure p_bcast_rttov_sinteger_4d
#endif
    module procedure p_bcast_rttov_coef_pccomp
    module procedure p_bcast_rttov_coef_pccomp1
    module procedure p_bcast_rttov_complex_1d
    ! if we don't use mo_mpi we must provide some mpi functions
    ! by ourselves
#ifdef RTTOV_IFC_USE_MPIF
    module procedure p_bcast_rttov_integer_1d
    module procedure p_bcast_rttov_integer_2d
    module procedure p_bcast_rttov_integer_3d
    module procedure p_bcast_rttov_integer_4d
    module procedure p_bcast_rttov_real_1d
    module procedure p_bcast_rttov_real_2d
    module procedure p_bcast_rttov_real_3d
    module procedure p_bcast_rttov_real_4d
    module procedure p_bcast_rttov_char_1d
    module procedure p_bcast_rttov_bool
#endif
  end interface

  interface p_bcast_rttov_container
    module procedure p_bcast_rttov_cnt_coef
    module procedure p_bcast_rttov_cnt_coef_scatt_ir
    module procedure p_bcast_rttov_cnt_optpar_ir
#if defined(RTTOV12)
    module procedure p_bcast_rttov_cnt_coef_optpiclb
    module procedure p_bcast_rttov_cnt_phasefn_int
    module procedure p_bcast_rttov_cnt_fast_coef
    module procedure p_bcast_rttov_cnt_fast_coef_gas
    module procedure p_bcast_rttov_cnt_nlte_coef
#ifndef RADSHARE
    module procedure p_bcast_rttov_cnt_telsem
    module procedure p_bcast_rttov_cnt_cnrm
#endif
    module procedure p_bcast_rttov_cnt_coef_pccomp2
#endif
    module procedure p_bcast_rttov_cnt_coef_pccomp
    module procedure p_bcast_rttov_cnt_coef_pccomp1
  end interface
#endif

#endif
  !-------------------
  ! version identifier
  !-------------------
#if defined(RTTOV12)
  integer      ,parameter :: RTTOV_IFC_VERSION   = 12
#elif defined(RTTOV10)
  integer      ,parameter :: RTTOV_IFC_VERSION   = 10
#else
  integer      ,parameter :: RTTOV_IFC_VERSION   = -1
#endif
  !-----------------------------
  ! error numbers of this module
  !-----------------------------
  integer ,parameter :: NO_ERROR             =  0 ! everything was ok.
  integer ,parameter :: ERR_ALLOC            =  1 ! allocation error
  integer ,parameter :: ERR_DIM              =  2 ! mismatch in dimensions
  integer ,parameter :: ERR_RTTOV_SETUP      =  3 ! error in rttov setup
  integer ,parameter :: ERR_CLOUD_AERO_MISSM =  4 ! aerosol and cloud classes change
                                                  ! for different coef_scatt_ir structures
  integer ,parameter :: ERROR_RTTOV_CALL     =  5 ! unable to calculate bt for all profiles
  integer ,parameter :: ERROR_RTTOV_CALL_ANY =  6 ! unable to calculate bt for some profiles
  integer ,parameter :: WARN_RTTOV_DIR_ALL   =  7 ! warning for all profiles
                                                  ! nevertheless bt were calculated
  integer ,parameter :: WARN_RTTOV_DIR_ANY   =  8 ! warning for some profiles
                                                  ! nevertheless bt were calculated
  integer ,parameter :: ERR_RTTOV_MPI        =  9 ! mpi comm. error
                                                  ! while reading/distributing the coeffs.
  integer ,parameter :: ERR_RTTOV_CL_AER     = 10 ! mismatch in rttov setup and profile
                                                  ! initialization (cloud/aerosol treatment)
  integer ,parameter :: ERR_NO_RTTOV_LIB     = 11 ! mismatch in rttov setup and profile
                                                  ! initialization (cloud/aerosol treatment)
  integer ,parameter :: ERR_RTTOV_PREC       = 12 ! mismatch in real kind precision of
                                                  !  input variables and rttov library.
  integer ,parameter :: ERR_CLOUD_INCONSIS   = 14 ! missing cloud fraction for cloud calc.
  integer ,parameter :: ERR_GOD_FILE         = 15 ! missing cloud fraction for cloud calc.
  integer ,parameter :: ERR_WR_PROF          = 16 ! missing cloud fraction for cloud calc.
  integer ,parameter :: ERR_INVALID_TSKIN    = 17 ! invalid t_skin
  integer ,parameter :: ERR_INPUT            = 18 ! Invalid/unsupported input
  integer ,parameter :: ERR_NO_ATLAS         = 19 ! No atlas support in current configuration
  integer ,parameter :: ERR_ATLAS_INIT       = 20 ! Atlas was not initialized
  integer ,parameter :: ERR_INVALID_INSTR    = 21 ! Invalid instrument (e.g. not supported by atlas)

  integer ,parameter :: nerr                 = 21 ! number of different error messages
  character(len=256) :: rttov_ifc_errMsg (nerr)   ! holds the error messages
                                                  ! for the errors above
  logical, save      :: init_errMsg = .false.     ! switch: set to true if
                                                  ! rttov_ifc_errMsg initialized
#ifndef NO_RTTOV
  !----------------------------------
  ! RTTOV relating internal parameter
  !----------------------------------
  integer(jpim),parameter :: RTTOV_ERR_NONE  = 0 ! no error messages output
  integer(jpim),parameter :: RTTOV_ERR_FATAL = 1 ! FATAL errors only printed. these are errors which
                                                 ! mean that profile should be aborted (e.g. unphysical
                                                 ! profile input)
  integer(jpim),parameter :: RTTOV_ERR_WARN  = 2 ! WARNING errors only printed. Errors which can allow
                                                 ! the computation to continue but the results may be
                                                 ! suspect (e.g. profile outside basis profile limits)
  integer(jpim),parameter :: RTTOV_ERR_INFO  = 3 ! INFORMATION messages which inform the user about
                                                 ! the computation
  integer(jpim),parameter :: usedRttov_info_Level =  RTTOV_ERR_WARN ! current rttov information level
                                                 ! RTTOV_ERR_INFO  should be used; for operations
                                                 ! RTTOV_ERR_FATAL should be used
  integer(jpim),parameter :: USE_DEFAULT_ERR_UNIT = -1
  integer(jpim),parameter :: RTTOV_ALLOC   = 1   ! allocation switch in call of rttov_alloc_prof and
                                                 ! rttov_alloc_rad;
  integer(jpim),parameter :: RTTOV_DEALLOC = 0   ! deallocation switch in call of rttov_alloc_prof and
                                                 ! rttov_alloc_rad;
  ! other parameters:
  real(wp)     ,parameter :: awful_rttov = -9999.0d0

  real(jprb)   ,save      :: qmin_ifc   = qmin_rttov
  real(jprb)   ,save      :: qmax_ifc   = qmax_rttov
  real(jprb)   ,save      :: tmin_ifc   = tmin_rttov
  real(jprb)   ,save      :: tmax_ifc   = tmax_rttov


  !---------
  ! switches
  !---------
  type(rttov_options),save,target :: opts    ! common options passed to rttov calls
  logical(jplm),save,allocatable :: &
    instr_addclouds(:),             & ! cloud switch per instrument
    instr_addaerosl(:)                ! cloud switch per instrument
  logical,       pointer :: switchrad_p         => null()
  integer(jpim), pointer :: fastem_version_p    => null()
  integer(jpim), pointer :: ir_sea_emis_model_p => null()
  logical(jplm), pointer :: apply_reg_lims_p    => null()
  logical(jplm), pointer :: verbose_reg_lims_p  => null()
  logical(jplm), pointer :: addclouds_p         => null()
  logical(jplm), pointer :: addaerosl_p         => null()
  logical(jplm), pointer :: addinterp_p         => null()
  logical(jplm), pointer :: addsolar_p          => null()
  logical(jplm), pointer :: addrefrac_p         => null()
  logical(jplm), pointer :: addpc_p             => null()
  logical(jplm), pointer :: use_q2m_p           => null()
  logical(jplm), pointer :: do_lambertian_p     => null()
  logical(jplm), pointer :: conv_overc_p        => null()
  logical(jplm), pointer :: ozone_data_p        => null()
  logical(jplm), pointer :: co2_data_p          => null()
  logical(jplm), pointer :: n2o_data_p          => null()
  logical(jplm), pointer :: co_data_p           => null()
  logical(jplm), pointer :: ch4_data_p          => null()
  logical(jplm), pointer :: clw_data_p          => null()

  !-------------------------------
  ! check on regularization limits
  !-------------------------------
  logical                    :: app_reg_lims        = .false. ! Apply regularization limits in RTTOV
  integer                    :: chk_reg_lims        = 0       ! Check regularization limits in rttov_ifc
                                                              ! bit1 (1): print results (invalid profiles)
                                                              ! bit2 (2): set flag for use in calling prog
  logical, save, allocatable :: mask_lims_t(:)                ! mask for t exceeding limit
  logical, save, allocatable :: mask_lims_q(:)                ! mask for t exceeding limit
  !---------------------------------
  ! check influence of god smoothing
  !---------------------------------
  integer                    :: chk_god             = 0       ! Check influence of god smoothing
                                                              ! bit1 (1): print results (invalid profiles)
                                                              ! bit2 (2): set flag for use in calling prog

#if defined(RTTOV12)
  integer, parameter         :: nprof_aux = 2
#elif defined(RTTOV10)
  integer, parameter         :: nprof_aux = 1
#endif
  type(rttov_profile),pointer:: profiles_aux(:) => NULL()     ! auxiliary profiles


  !----------------------------
  ! general dimension variables
  !----------------------------
  integer(jpim)          :: nprof              ! no.of profiles processed in one rttov call
  integer(jpim)          :: nmaxchansperinst   ! maximum number of channels of all instrs. wanted;
  integer(jpim)          :: nlevs              ! number of model levels of external grid;
  integer(jpim)          :: nlevs_top          ! additional levels added to top in case of RTTOV10;
  !-------------------------------------------
  ! general variables for RTTOV_... interfaces
  !-------------------------------------------
  integer(jpim)                             :: alloc_status(20)
  !------------------------------------
  ! variables for RTTOV_setup interface
  !------------------------------------
  type(rttov_coefs),target ,allocatable,save:: coefs(:)        ! container for coeffs
  integer,                  allocatable,save:: nchans_instr(:) ! number of channels for each instrument
  !-------------------------------------
  ! variables for RTTOV_direct interface
  !-------------------------------------
  type(rttov_profile),     pointer :: profiles    (:) => NULL()   ! atm. data
  type(rttov_transmission),save    :: transmission            ! transmittances,layer optical depths
  type(rttov_radiance)             :: radiance                ! radiances, brightness temperatures
#if defined(RTTOV12)
  type(rttov_radiance2)            :: radiance2               ! upselling and downwelling radiances(RTTOV12)
#endif
  real(kind=jprb),         pointer :: height_aux(:) => NULL() ! height of weighting function
  !--------------------------------
  ! variables for RTTOV_k interface
  !--------------------------------
  type(rttov_profile),     pointer :: profiles_k(:) => NULL() ! atm. data
  type(rttov_transmission),save    :: transmission_k          ! transmittances,layer optical depths
  type(rttov_radiance)             :: radiance_k              ! radiances, brightness temperatures
  !-------------------------------------------------------------
  ! define default values for rttov profiles_type initialization
  !-------------------------------------------------------------
  real(jprb)      ,parameter:: default_o3_surf           = 0.031438_jprb! o3 surface
  real(jprb)      ,parameter:: default_wfetch            = 100000._jprb ! wind fetch
  real(jprb)      ,parameter:: default_fastem(fastem_sp) = &            ! fastem coefficients relevant for land/ice
                             (/3.0_jprb,5.0_jprb,15.0_jprb,0.1_jprb,0.3_jprb/)
  ! TODO: consider switch to FASTEM6
  integer(jpim)   ,save     :: default_fastem_version    = 5            ! only for RTTOV10/12
#if defined(RTTOV12)
  integer(jpim)   ,save     :: default_ir_emis_version   = 2            ! only for RTTOV12
#else
  integer(jpim)   ,save     :: default_ir_emis_version   = 1            ! only for RTTOV10
#endif
  real(jprb)      ,parameter:: default_ctp               = 500.0_jprb   ! cloud top pressure
  real(jprb)      ,parameter:: default_cfraction         =   0.0_jprb   ! cloud fraction
  real(jprb)      ,parameter:: default_satazim           =   0.0_jprb   ! satellite azimuth angle
  real(jprb)      ,parameter:: default_sunzenangle       =   0.0_jprb   ! solar zenith angle
  real(jprb)      ,parameter:: default_sunazangle        =   0.0_jprb   ! solar azimuth angle
  integer(jpim)   ,parameter:: default_watertype         =   1_jpim     ! water type (fresh 0/ocean 1)
  integer(jpim)   ,parameter:: default_idg               =   4_jpim     ! Scheme for IWC to eff
  integer(jpim)   ,parameter:: default_ice_scheme        =   1_jpim     ! shape of the ice crystals
  logical         ,parameter:: default_addsolar          =  .false.     ! incl. reflected sol. rad.
  logical         ,parameter:: default_addrefrac         =  .true.      ! incl. refr. in path calc.
#if defined(RTTOV12)
  integer         ,save     :: default_gas_units         =  gas_unit_specconc
  logical         ,save     :: default_do_lambertian     = .false.
#endif
  logical         ,save     :: default_use_q2m           = .false.
  real(jprb)      ,save     :: default_salinity          = 0._jprb     ! salinity
  integer(jpim)   ,save     :: fix_hgpl                  = 0

#else /* defined(NO_RTTOV) */

  integer         ,save     :: default_fastem_version  = -1      ! otherwise
  integer         ,save     :: default_ir_emis_version = -1      ! otherwise
! real(jprb)      ,save     :: qmin_ifc                = -1._jprb
! real(jprb)      ,save     :: qmax_ifc                = -1._jprb
! real(jprb)      ,save     :: tmin_ifc                = -1._jprb
! real(jprb)      ,save     :: tmax_ifc                = -1._jprb

#endif
  !-------------------------------------------------------------
  ! Flag variables for output arrays
  !-------------------------------------------------------------
  integer         ,parameter:: OUT_ASB=0                         ! all sky brightness temp.
  integer         ,parameter:: OUT_CSB=1                         ! clear sky brightness temp.
  integer         ,parameter:: OUT_ASR=2                         ! all sky radiances
  integer         ,parameter:: OUT_CSR=3                         ! clear sky radiances

  !--------------
  ! MPI Variables
  !--------------
  integer                   :: mpi_my_proc_id         ! mpi id of this processor

  !--------------------------
  ! generalized optical depth
  !--------------------------
  character(len=300) ,save  :: god_par_file = ''
  logical            ,save  :: wr_god       = .false.

  character(len=300) ,save  :: out_path     = ''


  !--------------------------
  ! ME emissivity atlases
  !--------------------------
#if defined(RTTOV12) && ! defined(RADSHARE)
  type(rttov_emis_atlas_data), save :: mw_atlas(4)
#endif

#ifdef RADSHARE
  integer :: ipr_deb, pe_rt
#endif


!======================================================================
contains
!======================================================================

  function ini_rttov_ifc_errHandle()
  integer :: ini_rttov_ifc_errHandle
  !----------------------------------------------------------
  ! initialization routine for the rttov error message vector
  !----------------------------------------------------------
    integer :: irun1

    ini_rttov_ifc_errHandle = NO_ERROR
    if (.not.init_errMsg) then
       init_errMsg = .true.
       do irun1=1, size(rttov_ifc_errMsg)
          rttov_ifc_errMsg(irun1) = ''
       enddo
       rttov_ifc_errMsg (1)  = 'an allocation error occurred'
       rttov_ifc_errMsg (2)  = 'a mismatch in the dimensions of some arrays occurred'
       rttov_ifc_errMsg (3)  = 'an error in rttov setup occurred'
       rttov_ifc_errMsg (4)  = 'this error occurs if the aerosol and cloud classes change for different ' // &
                               'rttovids coef_scatt_ir structures rewrite code for single profiles ' // &
                               'allocation due to aerosol and cloud treatment'
       rttov_ifc_errMsg (5)  = 'rttov was not able to calculate the brightness temperatures for all ' // &
                               'of the input profiles'
       rttov_ifc_errMsg (6)  = 'rttov was not able to calculate the brightness temperatures for some ' // &
                               'of the input profiles'
       rttov_ifc_errMsg (7)  = 'a warning occurred for all of the input profiles; nevertheless the ' // &
                               'brightness temperatures (increments) were calculated'
       rttov_ifc_errMsg (8)  = 'a warning occurred for some of the input profiles; nevertheless the ' // &
                               'brightness temperatures (increments) were calculated'
       rttov_ifc_errMsg (9)  = 'an mpi communication error occurred while reading/distributing the coeffs.!'
       rttov_ifc_errMsg (10) = 'mismatch in call of rttov setup routine and profile initialization ' // &
                               '(cloud and/or aerosol treatment)!'
       rttov_ifc_errMsg (11) = 'the rttov library is not provided.'
       rttov_ifc_errMsg (12) = 'mismatch in real kind precision of input variables and rttov library.'
       rttov_ifc_errMsg (13) = 'more than one cloud type in one level'
       rttov_ifc_errMsg (14) = 'cloud cover and cloud water content arrays inconsistent'
       rttov_ifc_errMsg (15) = 'failed to read god_par_file.'
       rttov_ifc_errMsg (16) = 'failed to write hdf5 profile file.'
       rttov_ifc_errMsg (17) = 'some invalid t_surf/T_G/ts_fg.'
       rttov_ifc_errMsg (18) = 'invalid/unsupported input.'
       rttov_ifc_errMsg (19) = 'no atlas support in current configuration'
       rttov_ifc_errMsg (20) = 'atlas was not initialized'
       rttov_ifc_errMsg (21) = 'invalid instrument (e.g. not supported by atlas)'
    endif

#ifdef NO_RTTOV
    ini_rttov_ifc_errHandle = ERR_NO_RTTOV_LIB
#else
    if (jprb /= wp) ini_rttov_ifc_errHandle = ERR_RTTOV_PREC
#endif

  end function ini_rttov_ifc_errHandle
  !======================================================================

  function rttov_init(instruments,channels,nchansPerInst,my_proc_id,n_proc,  &
                      io_proc_id,mpi_comm_type,appRegLim,readAeros,readCloud,&
                      pathCoefs)
  integer,intent(in) :: instruments(:,:)       ! info on processed instruments:
                                               ! (1,1:nInstr): rttov platform IDs
                                               ! (2,1:nInstr): rttov satellite IDs
                                               ! (3,1:nInstr): rttov instrument IDs
                                               ! nInstr: number of all instruments
  integer,intent(in) :: channels(:,:)          ! channels to extract:
                                               ! (1:nchansPerInst,1:nInstr);
                                               ! numbers in channels are the indices
                                               ! of the channels in the coeff. file
  integer,intent(in) :: nchansPerInst(:)       ! number of channels for each instrument
                                               ! as defined in channels
  integer,intent(in) :: my_proc_id             ! ID of local processors
  integer,intent(in) :: n_proc                 ! no. of processors in mpi communication domain
  integer,intent(in) :: io_proc_id             ! ID of IO processor
  integer,intent(in) :: mpi_comm_type          ! mpi communicator type
  logical,intent(in),optional:: appRegLim      ! reset the profiles variables
                                               ! to the regression limits if outside
  logical,intent(in),optional:: readAeros(:)   ! read also the aerosol (ir scattering) coeffs.
  logical,intent(in),optional:: readCloud(:)   ! read also the cloud  (ir scattering) coeffs
  character(*),intent(in),optional:: pathCoefs ! path to coefficient files
  integer                         :: rttov_init! function result
  !--------------------------------------------------------------------------
  ! initialization routine which initializes the rttov modules
  ! and reads the instrument specific coefficients which are needed (wanted).
  ! in addition it fills the default values for some parts of the rttov
  ! profile structure such as Ozone, etc.
  !--------------------------------------------------------------------------
#ifndef NO_RTTOV
    integer(jpim),allocatable:: setup_errorstatus(:)
    integer(jpim)            :: nInstr
    integer(jpim)            :: run1
    integer(jpim)            :: error_status
#endif

    rttov_init = ini_rttov_ifc_errHandle()
    if (rttov_init /= NO_ERROR) return
#ifdef NO_RTTOV
    rttov_init = ERR_NO_RTTOV_LIB
    return
#else

FTRACE_BEGIN('rttov_init')

    ! initialize some module variables
    nprof            = 0
    nlevs            = 0
    nmaxchansperinst = int(maxval(nchansPerInst(:)),jpim)
    nInstr           = size(instruments, 2)
    alloc_status(:)  = 0_jpim

#if defined(RTTOV12)
!    opts%config%verbose = .true.
#ifndef RADSHARE
    opts%config%crop_k_reg_limits = .false. ! Option by DWD (RF):
                                  ! If the humidity input to rttov12 exceeds the regression limits, we get
                                  ! humi_k=0. above level_hum_dum, which is consistent. The humi_k values
                                  ! above level_hum_dum are not relevant for the assimilation. But they
                                  ! are written into the *RTOVP* files. In order to have realistic profiles in the
                                  ! *RTOVP* files, we prevent RTTOV12 from setting the results to zero, where
                                  ! the reg_limits were exceeded.
    conv_overc_p        => opts%config%conv_overc
    opts%config%fix_hgpl = fix_hgpl
#endif
    apply_reg_lims_p    => opts%config%apply_reg_limits
    verbose_reg_lims_p  => opts%config%verbose
    switchrad_p         => opts%rt_all%switchrad
    fastem_version_p    => opts%rt_mw%fastem_version
    ir_sea_emis_model_p => opts%rt_ir%ir_sea_emis_model
    addclouds_p         => opts%rt_ir%addclouds
    addaerosl_p         => opts%rt_ir%addaerosl
    addinterp_p         => opts%interpolation%addinterp
    addsolar_p          => opts%rt_ir%addsolar
    addrefrac_p         => opts%rt_all%addrefrac
    addpc_p             => opts%rt_ir%pc%addpc
    use_q2m_p           => opts%rt_all%use_q2m
    do_lambertian_p     => opts%rt_all%do_lambertian
    ozone_data_p        => opts%rt_ir%ozone_data
    co2_data_p          => opts%rt_ir%co2_data
    n2o_data_p          => opts%rt_ir%n2o_data
    co_data_p           => opts%rt_ir%co_data
    ch4_data_p          => opts%rt_ir%ch4_data
    clw_data_p          => opts%rt_mw%clw_data

#else
    apply_reg_lims_p    => opts%apply_reg_limits
    verbose_reg_lims_p  => opts%verbose_checkinput_warnings
    switchrad_p         => opts%switchrad
    fastem_version_p    => opts%fastem_version
    addclouds_p         => opts%addclouds
    addaerosl_p         => opts%addaerosl
    addinterp_p         => opts%addinterp
    addsolar_p          => opts%addsolar
    addrefrac_p         => opts%addrefrac
    addpc_p             => opts% addpc
    use_q2m_p           => opts% use_q2m
!   conv_overc_p        => opts%conv_overc
    ozone_data_p        => opts%ozone_data
    co2_data_p          => opts%co2_data
    n2o_data_p          => opts%n2o_data
    co_data_p           => opts%co_data
    ch4_data_p          => opts%ch4_data
    clw_data_p          => opts%clw_data
#endif

    mpi_my_proc_id = my_proc_id

    ! initialization of RTTOV_setup input variables
#if defined(RTTOV12)
    call rttov_errorhandling(USE_DEFAULT_ERR_UNIT)
#else
    call rttov_errorhandling(USE_DEFAULT_ERR_UNIT,usedRttov_info_Level)
#endif

    ! allocate, read and initialise coefficients
    allocate(setup_errorstatus(nInstr),stat=alloc_status(1))
    allocate(coefs            (nInstr),stat=alloc_status(2))
    allocate(instr_addclouds  (nInstr),stat=alloc_status(3))
    allocate(instr_addaerosl  (nInstr),stat=alloc_status(4))
    allocate(nchans_instr     (nInstr),stat=alloc_status(3))

    if(any(alloc_status /= 0)) then
       rttov_init = ERR_ALLOC
       return
    end if

    switchrad_p         = .true.
    fastem_version_p    = default_fastem_version
#if defined(RTTOV12)
    ir_sea_emis_model_p = default_ir_emis_version
    do_lambertian_p     = default_do_lambertian
#endif
    use_q2m_p           = default_use_q2m

    if (present(readCloud)) then
      instr_addclouds = readCloud(1:nInstr)
    else
      instr_addclouds = .False.
    endif

    if (present(readAeros)) then
      instr_addaerosl = readAeros(1:nInstr)
    else
      instr_addaerosl = .False.
    endif
    nchans_instr(1:nInstr) = nChansPerInst(1:nInstr)


    ! call RTTOV
    if (present(appRegLim)) app_reg_lims = appRegLim
    apply_reg_lims_p   = app_reg_lims
#if defined(RTTOV12)
    if (.not.apply_reg_lims_p) opts%config%do_checkinput = .false.
#endif
    verbose_reg_lims_p = .false.


    do run1 = 1, nInstr
      addclouds_p = instr_addclouds(run1)
      addaerosl_p = instr_addaerosl(run1)
#if defined(_RTTOV_DO_DISTRIBCOEF)
      if (io_proc_id == my_proc_id) then
#endif

#if defined(RTTOV12)
        call rttov_read_coefs(                                  &
             setup_errorstatus(run1),                           &! -> return code
             coefs(run1),                                       &! -> coefficients
             opts,                                              &! <- options
             channels   = channels(1:nchansPerInst(run1),run1), &! <- channels to load
             instrument = instruments(:,run1),                  &!
             path       = pathCoefs)

#else
        ! rttov_setup() requires that the coefficient files reside within
        ! the current directory. This is very ugly and we want to load
        ! the files from some arbitrary directory. Therefore, load them
        ! manually
        setup_errorstatus(run1) = rttov_read_coef_path(&
                 coefs(run1)                         ,&! -> coefs
                 instruments(:,run1),                 &! <- instr triplet.
                 channels(1:nchansPerInst(run1),run1),&! <- channels
                 pathCoefs)                            ! <- path
#endif

#if defined(_RTTOV_DO_DISTRIBCOEF)
      else
         call rttov_nullify_coef(coefs(run1) % coef)
         call rttov_nullify_coef_pccomp(coefs(run1) % coef_pccomp)
         call rttov_nullify_coef_scatt_ir(coefs(run1) % coef_scatt_ir)
         call rttov_nullify_optpar_ir(coefs(run1) % optp)
      endif

      ! distribute error status to all PEs
      if (n_proc > 1 ) &
           call p_bcast(setup_errorstatus(run1:run1), io_proc_id, mpi_comm_type)
#endif

      if (setup_errorstatus(run1) /= errorstatus_success) cycle

#if defined(_RTTOV_DO_DISTRIBCOEF)
      ! distribute loaded coefficients to all PEs
      if (n_proc > 1 ) &
           call p_bcast(coefs(run1),io_proc_id,mpi_comm_type)
#endif

#if ! defined(RTTOV12)
      call rttov_init_coefs(setup_errorstatus(run1),&! ->  return code
                            opts                   ,&! <-  options
                            coefs(run1)             )! <-> coefficients
#endif

   enddo

    if(any(setup_errorstatus(:) /= ERRORSTATUS_SUCCESS)) then
       print *, 'rttov_init: setup_errorstatus =', setup_errorstatus
       rttov_init = ERR_RTTOV_SETUP
       return
    endif

    if(my_proc_id == io_proc_id .and. usedRttov_info_Level >= RTTOV_ERR_INFO) then
       print *, 'successfully initialized RTTOV'
    endif

    if(associated(profiles)) then
#if defined(RTTOV12)
      error_status = dealloc_rttov_arrays(profs=profiles, rads=radiance, rads2=radiance2, transm=transmission)
#else
      error_status = dealloc_rttov_arrays(profs=profiles, rads=radiance, transm=transmission)
#endif
    end if
    if(associated(profiles_k)) &
      error_status = dealloc_rttov_arrays(profs=profiles_k, rads=radiance_k, transm=transmission_k)

    if (nInstr > 0) then
      nlevs = coefs(1)%coef%nlevels
      allocate(mask_lims_t(nlevs),mask_lims_q(nlevs))
      mask_lims_t = .true.
      mask_lims_q = .true.
    end if

#ifndef RADSHARE
#if defined(RTTOV12) 
    if (god_par_file /= '') then
      if (io_proc_id == my_proc_id) call read_god_par(setup_errorstatus(1))
#if defined(_RTTOV_DO_DISTRIBCOEF)
#ifndef NOMPI
      if (n_proc > 1 ) call p_bcast(setup_errorstatus(1:1), io_proc_id, mpi_comm_type)
#endif
#endif
      if (setup_errorstatus(1) /= 0) then
        rttov_init = setup_errorstatus(1)
        return
      end if
#if defined(_RTTOV_DO_DISTRIBCOEF)
      call bcast_god_par(my_proc_id,io_proc_id,mpi_comm_type)
#endif
    end if
#endif
#endif

    deallocate(setup_errorstatus)

FTRACE_END('rttov_init')


#endif

  end function rttov_init
  !======================================================================

#ifndef RADSHARE
#if defined(RTTOV12)
  subroutine read_god_par(stat)
    integer, intent(inout) :: stat

    type(rttov_coef),    pointer :: coef => null()
    integer                      :: iu, istat, i, icoef, i_entry, ichan, n_used
    logical                      :: ldum
    ! File content
    character(len=10000)         :: line
    character(len=100)           :: tag
    integer                      :: version
    ! god parameters
    integer                      :: instr, chan, layer, mode
    real(kind=jprb)              :: opdep, par1
    type(rttov_god_par), pointer :: gp

#define ERR stat = ERR_GOD_FILE ; close(iu) ; return
    stat = 0
    if (god_par_file /= '') then
      inquire(file=trim(god_par_file), exist=ldum)
      if (.not.ldum) then
        stat = ERR_GOD_FILE ; return
      end if
      iu = get_unit_number()
      open(iu, file=trim(god_par_file))
      version = 1
      i_entry = 0
      n_used  = 0
      read_god_loop: do
        read(iu, '(A)', iostat=istat) line
        !print*,'line',trim(line)
        if (istat /= 0) then
          if (istat == iostat_end) then
            write(stdout,*)
            write(stdout,*) 'Read god-parameters from "'//trim(god_par_file)//'":'
            write(stdout,*) '  #Entries          : ',i_entry
            write(stdout,*) '  #Useful entries   : ',n_used
            close(iu)
            exit
          else
            ERR
          end if
        end if
        line = adjustl(line)

        if (line == ''                            ) cycle read_god_loop ! Empty line
        if (line(1:1) == '#' .or. line(1:1) == '!') cycle read_god_loop ! Comment line

        if (i_entry == 0) then
          read(line,*,iostat=istat) tag, i
          if (istat == 0) then
            if (tag == 'version' .or. tag == 'Version') then
              version = i
              cycle read_god_loop
            end if
          end if
        end if

        istat = -1
        select case(version)
        case(1)
          read(line,*,iostat=istat) instr, chan, layer, mode, opdep, par1
        case default
          write(0,*) 'Invalid version ',version,' of god_par_file "'//trim(god_par_file)//'"'
          ERR
        end select
        if (istat == 0) then
          i_entry = i_entry + 1
          ldum = .false.
          do icoef = 1, size(coefs)
            coef => coefs(icoef)%coef
            if (coef%id_inst == instr) then
              if (layer < 1 .or. layer > coef%nlayers) then
                write(0,*) 'Invalid layer in god_par_file entry ',i_entry,': '//trim(line)
                ERR
              end if
              if (.not.associated(coef%god)) then
                allocate(coef%god(coef%nlayers, coef%fmv_chn))
                coef%god%ntr = 0
              end if
              do ichan = 1, coef%fmv_chn
                if (coef%ff_ori_chn(ichan) == chan) then
                  gp => coef%god(layer,ichan)
                  gp%ntr = gp%ntr + 1
                  if (gp%ntr > size(gp%tr)) then
                    write(0,*) 'Too many entries for instr=',instr,' chan=',chan,' layer=',layer,&
                         &' in god_par_file.'
                    ERR
                  end if
                  gp%tr(gp%ntr)%mode  = mode
                  gp%tr(gp%ntr)%opdep = opdep
                  gp%tr(gp%ntr)%par1  = par1
                  ldum = .true.
                end if
              end do
            end if
          end do
          if (ldum) then
            n_used = n_used + 1
          else
            print*,'not used entry: '//trim(line)
          end if
        else
          write(0,*) 'Invalid entry "',trim(line)//'" in god_par_file.'
          ERR
        end if
      end do read_god_loop
      close(iu)

      do icoef = 1, size(coefs)
        coef => coefs(icoef)%coef
        if (associated(coef%god)) then
          do ichan = 1, coef%fmv_chn
            do layer = 1, coef%nlayers
              gp => coef%god(layer, ichan)
              if (gp%ntr > 0) then
                call rttov_god_init(p=gp)
                if (wr_god) call write_god(gp, coef%id_inst, coef%ff_ori_chn(ichan), layer)
              end if
            end do
          end do
        end if
      end do

    end if
#undef ERR
  end subroutine read_god_par

#if defined(_RTTOV_DO_DISTRIBCOEF)

  subroutine bcast_god_par(my_proc_id,io_proc_id,mpi_comm_type)
    integer,          intent(in) :: my_proc_id             ! ID of processor
    integer,          intent(in) :: io_proc_id             ! ID of IO processor
    integer,          intent(in) :: mpi_comm_type          ! mpi communicator type
#ifndef NOMPI
    type(rttov_coef), pointer    :: coef => null()
    integer(kind=jpim)           :: icoef, ichan, layer
    logical                      :: l_io, l_assoc
    l_io = (io_proc_id == my_proc_id)
    do icoef = 1, size(coefs)
      coef => coefs(icoef)%coef
      l_assoc = associated(coef%god)
      call p_bcast(l_assoc,io_proc_id,mpi_comm_type)
      if (l_assoc) then
        if (.not.l_io) allocate(coef%god(coef%nlayers, coef%fmv_chn))
        do ichan = 1, coef%fmv_chn
          do layer = 1, coef%nlayers
            call p_bcast(coef%god(layer,ichan),io_proc_id,mpi_comm_type)
          end do
        end do
      end if
    end do
#endif
  end subroutine bcast_god_par

  subroutine p_bcast_god_par(buffer,source,comm)
  type(rttov_god_par), intent(inout)   :: buffer
  integer,             intent(in)      :: source
  integer,    optional,intent(in)      :: comm
  !------------------------------------------------------------------
  ! Broadcast an rttov_coef container across all available processors
  !------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#endif
  end subroutine p_bcast_god_par

#endif /* _RTTOV_DO_DISTRIBCOEF */

  subroutine write_god(god, instr, chan, layer)
    type(rttov_god_par), intent(in) :: god
    integer,             intent(in) :: instr
    integer,             intent(in) :: chan
    integer,             intent(in) :: layer

    integer, parameter  :: n = 1000
    integer             :: i, iu
    real(kind=wp)       :: x, dx, od
    character(len=100)  :: fname

    if (god%ntr > 0) then
      od = god%tr(1)%opdep
      dx = 6*od / (n-1)
      iu = get_unit_number()
      write(fname, '("god_instr",I3.3,"_chan",I5.5,"_layer",I3.3,".dat")') instr, chan, layer
      if (out_path /= '') fname = trim(out_path)//'/'//trim(fname)
      open(iu, file=trim(fname))
      print*,trim(fname),god%ntr
      do i = 1, n
        x = -2*od + (i-1) * dx
        write(iu, *) x, rttov_god2o(x, god)
      end do
      write(iu, *) '&'
      do i = 1, n
        x = -2*od + (i-1) * dx
        write(iu, *) x, rttov_d_god2o(god, god=x)
      end do

      write(iu, *) '&'
      do i = 1, n
        x = -2*od + (i-1) * dx
        write(iu, *) x, rttov_d_god2o(god, o=rttov_god2o(x, god))
      end do

      close(iu)
    end if
  end subroutine write_god

  function get_unit_number() result(iu)
    integer :: iu
    logical :: opened
    do iu = 10, 100000
      inquire(iu, opened=opened)
      if (.not.opened) return
    end do
  end function get_unit_number

#endif /* RTTOV12 */
#endif /* RADSHARE */
  !======================================================================


  function rttov_cleanup()
  integer :: rttov_cleanup ! function result
  !-----------------------------------------------------------------------
  ! the routine frees memory allocated by rttov_init (...)
  !-----------------------------------------------------------------------
#ifndef NO_RTTOV
    integer(jpim) :: irun1
    integer(jpim) :: errorstatus
    integer(jpim) :: nrttov_ids
#endif
    integer :: k

    rttov_cleanup = ini_rttov_ifc_errHandle()
    if (rttov_cleanup /= NO_ERROR) return

#ifdef NO_RTTOV
    rttov_cleanup = ERR_NO_RTTOV_LIB
    return
#else

    ! deallocation of rttov permanent arrays
    errorstatus = ERR_ALLOC
    if (associated(profiles)) then
#if defined(RTTOV12)
       errorstatus = dealloc_rttov_arrays(profs=profiles,rads=radiance,rads2=radiance2,transm=transmission)
#else
       errorstatus = dealloc_rttov_arrays(profs=profiles,rads=radiance,transm=transmission)
#endif
       if(errorstatus /= NO_ERROR) return
    endif
    if (associated(profiles_k)) then
       errorstatus = dealloc_rttov_arrays(profs=profiles_k,rads=radiance_k,transm=transmission_k)
       if(errorstatus /= NO_ERROR) return
    endif
    errorstatus = 0

    ! deallocate coefficients
    if ((allocated(coefs))) then
       nRttov_ids = size(coefs)
       do irun1=1,nRttov_ids
          call rttov_dealloc_coefs(errorstatus, coefs(irun1))
          if(errorstatus /= ERRORSTATUS_SUCCESS) then
             rttov_cleanup = ERR_ALLOC
             return
          endif
       enddo
       deallocate(coefs)
    endif

    if (allocated(mask_lims_t)) deallocate(mask_lims_t)
    if (allocated(mask_lims_q)) deallocate(mask_lims_q)
#endif

   ! Cleanup mw emis atlases
#if defined(RTTOV12) && ! defined(RADSHARE)
    do k = 1, size(mw_atlas)
      if (.not. mw_atlas(k)% init) cycle
      call rttov_deallocate_emis_atlas(mw_atlas(k))
    end do
#endif
  end function rttov_cleanup

  !======================================================================

  function rttov_get_preslev(instrIdx, preslev)
  integer,  intent(in)  :: instrIdx
  real(wp), intent(out) :: preslev(:)
  integer               :: rttov_get_preslev

#ifndef NO_RTTOV
    integer :: nlevs_user

    nlevs_user = size(preslev)

    rttov_get_preslev = 0
    if (nlevs_user > nlevs .or. nlevs_user < nlevs - 1) then
      rttov_get_preslev = ERR_DIM
      return
    end if
    if (.not.allocated(coefs)) then
      rttov_get_preslev = ERR_RTTOV_SETUP
      return
    end if
    if (size(coefs) < 1) then
      rttov_get_preslev = ERR_RTTOV_SETUP
      return
    end if

    nlevs_top = nlevs - nlevs_user
    preslev(:) = coefs(instrIdx)%coef%ref_prfl_p(1+nlevs_top:)
#else
    rttov_get_preslev = ERR_NO_RTTOV_LIB
#endif

  end function rttov_get_preslev
  !======================================================================

#ifndef NO_RTTOV

  subroutine rttov_prof_pc(prof,header,my_proc_id)
  type(rttov_profile),intent(in) :: prof      ! prof structure to be printed
  character(len=*)  ,intent(in) :: header    ! header string specifies the content of prof.
  integer           ,intent(in) :: my_proc_id! id of processor
  !-----------------------------------------------------------------------
  ! this subroutine prints the content of a rttov_profile structure
  ! c.f., module rttov_types in the rttov library for more information
  ! on this structure;
  !-----------------------------------------------------------------------
    logical,parameter        :: addaerosl1=.false.
    logical,parameter        :: addclouds1=.false.
    integer(jpim)            :: irun1,irun2
    Real(Kind=jprb),pointer  :: o3(:),co2(:),n2o(:),co(:),ch4(:),clw(:)
    Real(Kind=jprb),target   :: dummy(size(prof%p))

    dummy = -1.0_jprb

    o3  => dummy
    co2 => dummy
    n2o => dummy
    co  => dummy
    ch4 => dummy
    clw => dummy

    if(associated(prof% o3 )) o3  => prof% o3
    if(associated(prof% co2)) co2 => prof% co2
    if(associated(prof% n2o)) n2o => prof% n2o
    if(associated(prof% co )) co  => prof% co
    if(associated(prof% ch4)) o3  => prof% ch4
    if(associated(prof% clw)) o3  => prof% clw

    print *,my_proc_id,': content of rttov rttov_profile structure - ', &
                       trim(adjustl(header)),':'
    print *,my_proc_id,': ===========================================',&
                       repeat('=',len_trim(adjustl (header))),'='
    print *, my_proc_id, ': nlevels    : ', prof% nlevels
    print *, my_proc_id, ': idg        : ', prof% idg
#if defined(RTTOV12)
    print *, my_proc_id, ': ice_scheme : ', prof% ice_scheme
#else
    print *, my_proc_id, ': ish        : ', prof% ish
#endif
    print *, my_proc_id, ': zenangle   : ', prof% zenangle
    print *, my_proc_id, ': azangle    : ', prof% azangle
    print *, my_proc_id, ': sunzenangle: ', prof% sunzenangle
    print *, my_proc_id, ': sunazangle : ', prof% sunazangle
    print *, my_proc_id, ': elevation  : ', prof% elevation
    print *, my_proc_id, ': latitude   : ', prof% latitude
    print *, my_proc_id, ': ctp        : ', prof% ctp
    print *, my_proc_id, ': cfraction  : ', prof% cfraction

    do irun1=1,prof% nlevels
       print *, my_proc_id, ': ',                          &
            irun1, prof%p(irun1), prof%t(irun1),           &
            prof%q(irun1), o3(irun1)
       print *, my_proc_id, ': ',                          &
            irun1, co2(irun1), n2o(irun1),       &
            co(irun1), ch4(irun1), clw(irun1)
    enddo

    print *, my_proc_id, ': 2m/10m parameter:'
    print *, my_proc_id, ': ============='
    print *, my_proc_id, ': t2m  : ', prof% s2m% t
    print *, my_proc_id, ': q2m  : ', prof% s2m% q
    print *, my_proc_id, ': o2m  : ', prof% s2m% o
    print *, my_proc_id, ': p2m  : ', prof% s2m% p
    print *, my_proc_id, ': u10m : ', prof% s2m% u
    print *, my_proc_id, ': v10m : ', prof% s2m% v
    print *, my_proc_id, ': wfetc: ', prof% s2m% wfetc
    print *, my_proc_id, ': skin parameter:'
    print *, my_proc_id, ': ============='
    print *, my_proc_id, ': surftype : ', prof% skin% surftype
    print *, my_proc_id, ': watertype: ', prof% skin% watertype
    print *, my_proc_id, ': t        : ', prof% skin% t
    print *, my_proc_id, ': fastem   : ', prof% skin% fastem(:)
    if(addaerosl1) then
       print *, my_proc_id, ': aerosol parameter:'
       print *, my_proc_id, ': =================='
       do irun1=1,prof% nlevels
          do irun2 = 1, size(prof% aerosols,1)
             print *, my_proc_id, ': ', irun1, irun2, ': ', prof% aerosols(irun2,irun1)
          enddo
       enddo
    endif
    if( addclouds1 ) then
       print *, my_proc_id, ': cloud parameter:'
       print *, my_proc_id, ': ================'
       do irun1=1,prof% nlevels
          do irun2=1,size(prof% cloud,1)
             print *, my_proc_id, ': ',  irun1, irun2, ': ', prof% cloud(irun2,irun1)
          enddo
       enddo
    endif

  end subroutine rttov_prof_pc
  !======================================================================

  subroutine rttov_rad_pc(rad,header,my_proc_id)
  type (rttov_radiance),intent(in) :: rad
  character (len=*)   ,intent(in) :: header
  integer             ,intent(in) :: my_proc_id
  !-----------------------------------------------------------------------
  ! this subroutine prints the content of a rttov_radiance structure
  ! c.f., module rttov_types in the rttov library for more information
  ! on this structure;
  !-----------------------------------------------------------------------
    integer(jpim) :: irun1,irun2
    integer(jpim) :: n_levs, nchans

    nchans  = size(rad% clear)
    n_levs  = size(rad% overcast,1)

    print *, my_proc_id, ': content of rttov rttov_radiance structure - ', &
                         trim(adjustl(header)),':'
    print *, my_proc_id, ': ===========================================', &
                         repeat('=', len_trim(adjustl(header))),'='
    print *, my_proc_id, ': nlevs  : ', n_levs
    print *, my_proc_id, ': nchans : ', nchans

    do irun1=1,nchans
       print *, my_proc_id, ': ',irun1,rad%clear(irun1),rad%cloudy(irun1),&
                                 rad%total(irun1),rad%bt(irun1)
#if defined(RTTOV12)
       print *, my_proc_id, ': ',irun1,rad%bt_clear(irun1),rad%refl_clear(irun1)
#else
       print *, my_proc_id, ': ',irun1,rad%bt_clear(irun1),rad%upclear(irun1),&
                                 rad%dnclear(irun1),rad%reflclear(irun1)
#endif
       do irun2=1,n_levs
          print *, my_proc_id, ': ',irun1,irun2,rad%overcast(irun2,irun1)
       enddo
    enddo

  end subroutine rttov_rad_pc

#endif
  !======================================================================

  function rttov_fill_input (press,temp,humi,t2m,q2m,psurf,hsurf,u10m,v10m,stemp,    &
                             stype,lat,lon,satzenith,sunZenith,                      &
#if !defined(RADSHARE) && !defined(NO_RTTOV)
                             rad,                                                    &
#endif
                             fastem,satazim,   &
                             sunazangle,o3,co2,n2o,co,ch4,clw,aerosols,cloud,cfrac,  &
                             idg,ice_scheme,watertype,o3_surf,wfetc,ctp,cfraction,addsolar, &
                             addrefrac,addinterp,init,ivect,rttov9_compat,istart,    &
                             iend, pe, wr_profs, wr_profs_fmt)
  ! --- Either these arguments ...
  real(wp),     intent(in), optional, target :: press     (:,:) ! press. grid (nlevs,nprofs) [hPa]
  real(wp),     intent(in), optional, target :: temp      (:,:) ! temperature profiles (nlevs,nprofs)[K]
  real(wp),     intent(in), optional, target :: humi      (:,:) ! water vapour profile (nlevs,nprofs) [ppmv]
  real(wp),     intent(in), optional, target :: t2m         (:) ! 2m temperature (nprofs) [K]
  real(wp),     intent(in), optional, target :: q2m         (:) ! 2m humidity (nprofs) [ppmv]
  real(wp),     intent(in), optional, target :: psurf       (:) ! surface pressure (nprofs) [hPa]
  real(wp),     intent(in), optional, target :: hsurf       (:) ! surface height (nprofs) [km]
  real(wp),     intent(in), optional, target :: u10m        (:) ! 10m U wind component (nprofs) [m/s]
  real(wp),     intent(in), optional, target :: v10m        (:) ! 10m V wind component (nprofs) [m/s]
  real(wp),     intent(in), optional, target :: stemp       (:) ! radiative skin temperature (nprofs) [K]
  integer ,     intent(in), optional, target :: stype       (:) ! surface type (nprofs):0,1,2=land,sea,sea-ice
  real(wp),     intent(in), optional, target :: lat         (:) ! latitude of ground point (nprofs) [deg]
  real(wp),     intent(in), optional, target :: lon         (:) ! longitude (nprofs) [deg]
  real(wp),     intent(in), optional, target :: satzenith   (:) ! local satellite zenith angle (nprofs) [deg]
  real(wp),     intent(in), optional, target :: sunZenith   (:) ! local solar zenith angle (nprofs) [deg]
  real(wp),     intent(in), optional, target :: ctp         (:) ! cloud top pressure (of a black body cloud)
                                                                !   (nprofs) [hPa] default: 500
  real(wp),     intent(in), optional, target :: cfraction   (:) ! cloud fraction (of a black body cloud)
                                                                !   (nprofs) range (0-1) default: 0
#if !defined(RADSHARE) && !defined(NO_RTTOV)
  ! --- or this must be present
  type(t_radv), intent(in), optional, target :: rad      ! derived type to store all the above info
#endif
  ! --------------------------------------

  real(wp),     intent(in), optional, target :: satazim     (:) ! local satellite azimuth (nprofs) [deg]
                                                                !   range: 0-360; east=90
  real(wp),     intent(in), optional, target :: sunazangle  (:) ! local solar azimuth angle (nprofs) [deg]
                                                                ! range: 0-360; east=90
  real(wp),     intent(in), optional         :: fastem      (:) ! fastem-2/3  land/sea-ice surf. param.(fastem_sp)
                                                                !   default: (/3.0, 5.0, 15.0, 0.1, 0.3/)
  real(wp),     intent(in), optional         :: o3        (:,:) ! ozone (nlevs,nprofs) [ppmv]
  real(wp),     intent(in), optional         :: co2       (:,:) ! carbon dioxide (nlevs,nprofs) [ppmv]
  real(wp),     intent(in), optional         :: n2o       (:,:) ! nitrous oxide (nlevs,nprofs) [ppmv]
  real(wp),     intent(in), optional         :: co        (:,:) ! carbon monoxide (nlevs,nprofs) [ppmv]
  real(wp),     intent(in), optional         :: ch4       (:,:) ! methane (nlevs,nprofs) [ppmv]
  real(wp),     intent(in), optional         :: clw       (:,:) ! cld liquid H2O (MW only) (nlevs,nprofs) [kg/kg]
  real(wp),     intent(in), optional         :: aerosols(:,:,:) ! aerosols (nAerosPar,nlevs,nprofs) [cm^-3]
  real(wp),     intent(in), optional, target :: cloud   (:,:,:) ! cloud water/ice (IR only)
                                                                ! (nCloudPar,nlevs,nprofs) [g m^-3]
  real(wp),     intent(in), optional, target :: cfrac     (:,:) ! cloud fractional cover (IR only)
                                                                ! (nCloudPar,nlevs,nprofs) range (0-1)
  integer ,     intent(in), optional         :: idg         (:) ! Scheme for IWC to eff. diameter (nprofs)
                                                                !   1: Scheme by Ou and Liou,
                                                                !      1995, Atmos. Res., 35, 127-138.
                                                                !   2: Scheme by Wyser et al.
                                                                !      (see McFarquhar et al. (2003))
                                                                !   3: Scheme by Boudala et al.,
                                                                !      2002 , Int. J. Climatol., 22, 1267-1284
                                                                !   4: Scheme by McFarquhar et al. (2003)
                                                                !   default: 4
  integer ,     intent(in), optional         :: ice_scheme  (:) ! shape of the ice crystals (nprofs)
                                                                !   1: hexagonal; 2: ice aggregates; default: 1
  integer ,     intent(in), optional         :: watertype   (:) ! water type of ground point (nprofs)
                                                                !   0: fresh water 1: ocean water  default: 1
  real(wp),     intent(in), optional         :: o3_surf     (:) ! surface ozone concentration (nprofs) [ppmv]
  real(wp),     intent(in), optional         :: wfetc       (:) ! wind fetch (nprofs) [m] default: 100000
  logical ,     intent(in), optional         :: addsolar        ! incl. reflected sol. rad. default: .false.
  logical ,     intent(in), optional         :: addrefrac       ! incl. refr. in path calc. default: .true.
  logical ,     intent(in), optional         :: addinterp       ! switch for interpolator default: .true.
  logical ,     intent(in), optional         :: init            ! initialization of profile structures
                                                                !   after alloc. components default: .true.
  integer,      intent(in), optional         :: ivect           ! Vectorization
                                                                ! ivect =1 : vectorize profiles
                                                                ! ivect/=1 : vectorize levels
  logical ,     intent(in), optional         :: rttov9_compat   ! Add 0.005 hPa level to RTTOV10 default: .true..
  integer ,     intent(in), optional         :: istart   ! first profile (in rad) to be used
  integer ,     intent(in), optional         :: iend     ! last  profile (in rad) to be used
  integer ,     intent(in), optional         :: pe
  integer ,     intent(in), optional         :: wr_profs(:)
  character(len=*),intent(in), optional      :: wr_profs_fmt


  integer                                    :: rttov_fill_input! function result
  !-----------------------------------------------------------------------
  ! this routine fills the variable part for a call of rttov for a
  ! rttov profile structure array (for processing more than one profile at once).
  ! in addition it allows you to change the default initialization values
  ! (set by rttov_init) via several optional arguments
  !-----------------------------------------------------------------------
#ifndef NO_RTTOV
    logical           :: init1
    integer(jpim)     :: errorstatus
    integer(jpim)     :: irun1
    integer(jpim)     :: ivect_loc
    integer(jpim)     :: i, j
    integer(jpim)     :: n_levs
    integer(jpim)     :: n_levs_top
    integer(jpim)     :: n_lev_layers
    integer           :: is, ie
    real(wp), pointer :: press_const_p (:) => null()
    real(wp), pointer :: press_p     (:,:) => null()
    real(wp), pointer :: temp_p      (:,:) => null()
    real(wp), pointer :: humi_p      (:,:) => null()
    real(wp), pointer :: t2m_p         (:) => null()
    real(wp), pointer :: q2m_p         (:) => null()
    real(wp), pointer :: psurf_p       (:) => null()
    real(wp), pointer :: hsurf_p       (:) => null()
    real(wp), pointer :: u10m_p        (:) => null()
    real(wp), pointer :: v10m_p        (:) => null()
    real(wp), pointer :: stemp_p       (:) => null()
    real(wp), pointer :: ctp_p         (:) => null() ! simple blackbody cloud top
    real(wp), pointer :: clf_p         (:) => null() ! simple blackbody cloud top
    real(wp), pointer :: cld_p     (:,:,:) => null() ! complex cloud water content
    real(wp), pointer :: cfr_p       (:,:) => null() ! complex cloud fraction
    integer , pointer :: stype_p       (:) => null()
    real(wp), pointer :: lat_p  (:)        => null()
    real(wp), pointer :: lon_p  (:)        => null()
    real(wp), pointer :: satzenith_p   (:) => null()
    real(wp), pointer :: sunZenith_p   (:) => null()
    real(wp), pointer :: satazim_p     (:) => null()
    real(wp), pointer :: sunazangle_p  (:) => null()
    integer , pointer :: id_p          (:) => null()
    logical           :: l_input
    ! writing of profiles
    integer           :: stat
    character(len=300):: fname             =  ""
#endif
    integer :: pe_loc

    rttov_fill_input = ini_rttov_ifc_errHandle()
    if (rttov_fill_input /= NO_ERROR) return

#ifdef NO_RTTOV
    rttov_fill_input = ERR_NO_RTTOV_LIB
    return
#else

    if (present(pe)) then
      pe_loc = pe
    else
      pe_loc = -1
    end if
    pe_rt = pe_loc


#define DIM_ERROR(text) rttov_fill_input = ERR_DIM ; print*,text ; return

FTRACE_BEGIN('rttov_fill_input')
      satazim_p    => null()
      sunazangle_p => null()
      id_p         => null()

      l_input = .false.

#ifndef RADSHARE
    if (present(rad)) then
      l_input = .true.
      if (present(istart)) then
        is = istart
      elseif (associated(rad%t_fg)) then
        is = lbound(rad%t_fg, 2)
      else
        DIM_ERROR('rad%t_fg lbound')
      end if
      if (present(iend)) then
        ie = iend
      elseif (associated(rad%t_fg)) then
        ie = ubound(rad%t_fg, 2)
      else
        DIM_ERROR('rad%t_fg ubound')
      end if
      if (size(rad%p, 2) == 1) then
        press_const_p=> rad% p(:,1)
        press_p      => NULL()
      else
        press_p      => rad% p
        press_const_p=> NULL()
      end if
      if (.not.associated(rad%t_fg  )) then ; DIM_ERROR('rad%t_fg') ; endif
      temp_p       => rad% t_fg(:,is:ie)
      if (.not.associated(rad%q_fg  )) then ; DIM_ERROR('rad%q_fg') ; endif
      humi_p       => rad% q_fg(:,is:ie)
      if (.not.associated(rad%t2m   )) then ; DIM_ERROR('rad%t2m') ; endif
      t2m_p        => rad% t2m(is:ie)
      if (.not.associated(rad%q2m   )) then ; DIM_ERROR('rad%q2m') ; endif
      q2m_p        => rad% q2m(is:ie)
      if (.not.associated(rad%ps_fg )) then ; DIM_ERROR('rad%ps_fg') ; endif
      psurf_p      => rad% ps_fg(is:ie)
      if (.not.associated(rad%shgt  )) then ; DIM_ERROR('rad%shgt') ; endif
      hsurf_p      => rad% shgt(1,is:ie)
      if (.not.associated(rad%u10_fg)) then ; DIM_ERROR('rad%u10_fg') ; endif
      u10m_p       => rad% u10_fg(is:ie)
      if (.not.associated(rad%v10_fg)) then ; DIM_ERROR('rad%v10_fg') ; endif
      v10m_p       => rad% v10_fg(is:ie)
      if (.not.associated(rad%ts_fg )) then ; DIM_ERROR('rad%ts_fg') ; endif
      stemp_p      => rad% ts_fg(is:ie)
      if (.not.associated(rad%stype )) then ; DIM_ERROR('rad%stype') ; endif
      stype_p      => rad% stype(1,is:ie)
      if (.not.associated(rad%dlat  )) then ; DIM_ERROR('rad%dlat') ; endif
      lat_p        => rad% dlat(is:ie)
      if (.not.associated(rad%dlon  )) then ; DIM_ERROR('rad%dlon') ; endif
      lon_p        => rad% dlon(is:ie)
      if (.not.associated(rad%stzen )) then ; DIM_ERROR('rad%stzen') ; endif
      satzenith_p  => rad% stzen(is:ie)

      if (associated(rad% cld_top) .and. associated(rad% cld_frc)) then
        ctp_p      => rad% cld_top(is:ie)
        clf_p      => rad% cld_frc(1,is:ie)
      elseif (associated(rad% cld_fg) .and. associated(rad% cld_frc)) then
        cld_p      => rad% cld_fg (:,:,is:ie)
        cfr_p      => rad% cld_frc(  :,is:ie)
      end if
      if (associated(rad% sunzen))  sunZenith_p  => rad% sunzen(is:ie)
      if (associated(rad% stazi))   satazim_p    => rad% stazi(is:ie)
      if (associated(rad% sunazi))  sunazangle_p => rad% sunazi(is:ie)
      if (associated(rad% obsnum))  id_p         => rad% obsnum(is:ie)
    end if
#endif
    if (.not.l_input) then
      if ( present(press)      .and. present(temp)      .and. &
           present(humi)       .and. present(stemp)     .and. &
           present(t2m)        .and. present(q2m)       .and. &
           present(v10m)       .and. present(u10m)      .and. &
           present(psurf)      .and. present(hsurf)     .and. &
           present(lat)        .and. present(satzenith) .and. &
           present(stype))    then
        l_input = .true.
        press_p      => press
        press_const_p=> NULL()
        temp_p       => temp
        humi_p       => humi
        t2m_p        => t2m
        q2m_p        => q2m
        psurf_p      => psurf
        hsurf_p      => hsurf
        u10m_p       => u10m
        v10m_p       => v10m
        stemp_p      => stemp
        stype_p      => stype
        lat_p        => lat
        if (present(lon)) lon_p => lon
        satzenith_p  => satzenith

        if (present(ctp) .and. present(cfraction)) then
          ctp_p      => ctp
          clf_p      => cfraction
        elseif (present(cloud) .and. present(cfrac)) then
          cld_p      => cloud
          cfr_p      => cfrac(:,:)
        end if
        if (present(satazim))    satazim_p    => satazim
        if (present(sunazangle)) sunazangle_p => sunazangle
        if (present(sunzenith))  sunZenith_p  => sunzenith
      end if
    end if
    if (.not.l_input) then
      DIM_ERROR('no input')
    end if


    ! switch for the interpolator
    addinterp_p = .true.
    if(present(addinterp)) addinterp_p = addinterp

    n_levs     = size(temp_p, 1)
    nprof      = size(temp_p, 2)

    n_levs_top   = 0
    n_lev_layers = n_levs - 1
    if (.not.addinterp_p) then
      if ( .not. present(rttov9_compat)) then
        n_levs_top = 1
      else
        if (rttov9_compat) then
          n_levs_top = 1
        endif
      endif
      n_lev_layers = n_levs + n_levs_top - 1
    end if


    errorstatus = 0

    ! Initialize structures after alloc. components
    init1 = .true.
    if(present(init)) init1 = init

    ! vectorization
    ivect_loc = 0
    if (present(ivect)) ivect_loc = ivect

    if (associated(press_p)) then
      if (size(press_p,    1)   /= n_levs) then
        DIM_ERROR('press_p')
      end if
    elseif (associated(press_const_p)) then
      if (size(press_const_p,1) /= n_levs)  then
        DIM_ERROR('press_const_p')
      end if
    else
      DIM_ERROR('no pressure')
    end if

    if   (size(humi_p,    1) /= n_levs) then
      DIM_ERROR('humi_p n_levs')
    end if
    if (( size(temp_p,    2) /= nprof) .or. &
         (size(humi_p,    2) /= nprof) .or. &
         (size(t2m_p       ) /= nprof) .or. &
         (size(q2m_p       ) /= nprof) .or. &
         (size(psurf_p     ) /= nprof) .or. &
         (size(u10m_p      ) /= nprof) .or. &
         (size(v10m_p      ) /= nprof) .or. &
         (size(stemp_p     ) /= nprof) .or. &
         (size(lat_p       ) /= nprof) .or. &
         (size(satzenith_p ) /= nprof) .or. &
         (size(stype_p     ) /= nprof)) then
      DIM_ERROR('nprof')
    end if

    if (present(idg))        then
       if (size(idg)         /= nprof) then
         write(0,*) size(idg), nprof
         DIM_ERROR('idg')
       end if
    endif
    if (present(ice_scheme))        then
       if (size(ice_scheme)  /= nprof) then
         DIM_ERROR('ice_scheme')
       end if
    endif
    if (present(watertype))  then
       if (size(watertype)   /= nprof) then
         DIM_ERROR('watertype')
       end if
    endif
    if (associated(satazim_p))    then
       if (size(satazim_p)     /= nprof) then
         DIM_ERROR('satazim')
       end if
    endif
    if (associated(sunzenith_p)) then
       if (size (sunzenith_p) /= nprof) then
         DIM_ERROR('sunzenith')
       end if
    endif
    if (associated(sunazangle_p)) then
       if (size (sunazangle_p) /= nprof) then
         DIM_ERROR('sunazangle')
       end if
    endif
    if (present(o3_surf))    then
       if (size(o3_surf)     /= nprof) then
         DIM_ERROR('o3_surf')
       end if
    endif
    if (present(wfetc))      then
       if (size (wfetc)      /= nprof) then
         DIM_ERROR('wfetc')
       end if
    endif
    if (present(ctp))        then
       if (size(ctp)         /= nprof) then
         DIM_ERROR('ctp')
       end if
    endif
    if (present(cfraction))  then
       if (size(cfraction)   /= nprof) then
         DIM_ERROR('cfraction')
       end if
    endif
    if (present(fastem))     then
       if (size(fastem)      /= fastem_sp) then
         DIM_ERROR('fastem')
       end if
    endif
    if (present(o3)) then
       if ((size(o3,1) /=n_levs).and.(size(o3 ,2)/=nprof)) then
         DIM_ERROR('o3')
       end if
    endif
    if (present(co2)) then
       if ((size(co2,1)/=n_levs).and.(size(co2,2)/=nprof)) then
         DIM_ERROR('co2')
       end if
    endif
    if (present(n2o)) then
       if ((size(n2o,1)/=n_levs).and.(size(n2o,2)/=nprof)) then
         DIM_ERROR('n2o')
       end if
    endif
    if (present(co))  then
       if ((size(co,1) /=n_levs).and.(size(co ,2)/=nprof)) then
         DIM_ERROR('co')
       end if
    endif
    if (present(ch4)) then
       if ((size(ch4,1)/=n_levs).and.(size(ch4,2)/=nprof)) then
         DIM_ERROR('ch4')
       end if
    endif
    if (present(clw)) then
       if ((size(clw,1)/=n_levs).and.(size(clw,2)/=nprof)) then
         DIM_ERROR('clw')
       end if
    endif

    if (present(aerosols)) then
       if ((size(aerosols,1)>max(1_jpim,maxval(coefs(:) % coef_scatt_ir % fmv_aer_comp))) .or.    &
           (size(aerosols,2)/=n_lev_layers).or.(size(aerosols,3)/=nprof)) then
         DIM_ERROR('aerosols')
       end if
    endif
    if (present(cloud)) then
       if ((size(cloud,1)>max (1_jpim,maxval(coefs(:) % coef_scatt_ir % fmv_wcl_comp+1_jpim))) .or.    &
           (size(cloud,2)/=n_lev_layers).or.(size(cloud,3)/=nprof)) then
         DIM_ERROR('cloud')
       end if
    endif
!#if defined(RTTOV12)
    if (present(cfrac)) then
       if ((size(cfrac,1)/=n_lev_layers).or.(size(cfrac,2)/=nprof)) then
         DIM_ERROR('cfrac')
       end if
    endif
! #else
!     if (present(cfrac)) then
!        if ((size(cfrac,1)>ncldtyp) .or.    &
!            (size(cfrac,2)/=n_lev_layers).or.(size(cfrac,3)/=nprof)) then
!          DIM_ERROR('cfrac')
!        end if
!     endif
! #endif
    rttov_fill_input = 0

#undef DIM_ERROR

    if (present(cloud) .neqv. present(cfrac)) then
       rttov_fill_input = ERR_CLOUD_INCONSIS
       return
    endif

    if (present(addsolar)) then
      addsolar_p = addsolar
    else
      addsolar_p = default_addsolar
    endif
    if (present(addrefrac)) then
      addrefrac_p = addrefrac
    else
      addrefrac_p = default_addrefrac
    endif
    ozone_data_p = present(o3)
    co2_data_p   = present(co2)
    n2o_data_p   = present(n2o)
    co_data_p    = present(co)
    ch4_data_p   = present(ch4)
    clw_data_p   = present(clw)
    addaerosl_p  = present(aerosols)
    addclouds_p  = present(cloud)



    !-----------------------------------------------
    ! allocate the direct structures if not yet done
    !-----------------------------------------------
    errorstatus = realloc_rttov_arrays( nprof,                   &
                                        n_levs + n_levs_top,     &
                                        nprof * nmaxchansperinst,&
                                        profs    = profiles      )

    if (errorstatus /= NO_ERROR) then
      rttov_fill_input    = int(errorstatus)
      nlevs     = 0
      nlevs_top = 0
      return
    else
      nlevs     = n_levs
      nlevs_top = n_levs_top
    endif

    !-------------------------
    ! fill RTTOV profile input
    !-------------------------
    profiles(1:nprof)% nlevels  = n_levs + n_levs_top
    profiles(1:nprof)% nlayers   = n_levs + n_levs_top - 1

#if defined(RTTOV12)
    profiles(1:nprof)% gas_units  = default_gas_units
    profiles(1:nprof)% mmr_cldaer = .false. ! g/cm^3 instead of kg/kg
#endif


    do irun1=1,nprof
      profiles(irun1)%p = huge(1._wp)
    end do

    if (ivect_loc == 1) then ! Vectorize profiles
      if (.not.addinterp_p .and. n_levs_top==1) then
        do irun1=1,nprof
          ! profile variables:
          profiles(irun1)% p(1) = 0.005_jprb
          profiles(irun1)% t(1) = min(tmax_ifc,max(tmin_ifc,temp_p(1,irun1)))
          profiles(irun1)% q(1) = min(qmax_ifc,max(qmin_ifc,humi_p(1,irun1)))
        enddo
      end if

!cdir novector
      do i=1,n_levs
        if (associated(press_p)) then
!cdir vector
          do irun1=1,nprof
            ! profile variables:
            profiles(irun1)% p(i+nlevs_top) = press_p(i,irun1)
          end do
        else
!cdir vector
          do irun1=1,nprof
            profiles(irun1)% p(i+nlevs_top) = press_const_p(i)
          end do
        end if
!cdir vector
        do irun1=1,nprof
          ! profile variables:
          profiles(irun1)% t(i+nlevs_top) = min(tmax_ifc,max(tmin_ifc,temp_p (i,irun1)))
          profiles(irun1)% q(i+nlevs_top) = min(qmax_ifc,max(qmin_ifc,humi_p(i,irun1)))
        enddo
      enddo

!cdir vector
       do irun1=1,nprof
          if (associated(id_p)) then
            write(profiles(irun1)% id,*) id_p(irun1)
          end if
          ! 2 meter air variables:
          profiles(irun1)% s2m% t          = min(tmax_ifc,max(tmin_ifc,t2m_p(irun1)))
          profiles(irun1)% s2m% q          = min(qmax_ifc,max(qmin_ifc,q2m_p(irun1)))
          profiles(irun1)% s2m% p          = psurf_p     (irun1)
          profiles(irun1)% s2m% u          = u10m_p      (irun1)
          profiles(irun1)% s2m% v          = v10m_p      (irun1)

          ! skin variables
          profiles(irun1)% skin% t         = stemp_p     (irun1)
          profiles(irun1)% skin% surftype  = int(stype (irun1), jpim)
          profiles(irun1)% skin% salinity  = default_salinity

          ! cloud top pressure and fraction
          if (associated(ctp_p)) then
            profiles(irun1)% ctp           = ctp_p       (irun1)
          else
            profiles(irun1)% ctp           = default_ctp
          end if
          if (associated(clf_p)) then
            profiles(irun1)% cfraction     = clf_p       (irun1)
          else
            profiles(irun1)% cfraction     = default_cfraction
          end if

          ! other variables from interactive inputs
          profiles(irun1)% elevation       = hsurf_p     (irun1)
          profiles(irun1)% zenangle        = abs(satzenith_p (irun1))
          profiles(irun1)% latitude        = lat_p(irun1)
          if (associated(lon_p)) then
            profiles(irun1)% longitude     = lon_p(irun1)
          end if

          if (present(idg)) then
             profiles(irun1)% idg = int(idg(irun1),jpim)
          else
             profiles(irun1)% idg = default_idg
          endif

#if defined(RTTOV12)
          if (present(ice_scheme)) then
             profiles(irun1)% ice_scheme = int(ice_scheme(irun1),jpim)
          else
             profiles(irun1)% ice_scheme = default_ice_scheme
          endif
#else
          if (present(ice_scheme)) then
             profiles(irun1)% ish        = int(ice_scheme(irun1),jpim)
          else
             profiles(irun1)% ish        = default_ice_scheme
          endif
#endif

          if (present(watertype)) then
             profiles(irun1)% skin% watertype = int(watertype(irun1),jpim)
          else
             profiles(irun1)% skin% watertype = default_watertype
          endif

          if (present(fastem)) then
             profiles(irun1)% skin% fastem = fastem
          else
             profiles(irun1)% skin% fastem = default_fastem
          endif

          if (associated(satazim_p)) then
             profiles(irun1)% azangle = satazim_p(irun1)
          else
             profiles(irun1)% azangle = default_satazim
          endif

          if (associated(sunazangle_p)) then
             profiles(irun1)% sunazangle = sunazangle_p(irun1)
          else
             profiles(irun1)% sunazangle = default_sunazangle
          endif

          if (associated(sunzenith_p)) then
             profiles(irun1)% sunzenangle = sunzenith_p(irun1)
          else
             profiles(irun1)% sunzenangle = default_sunzenangle
          endif

          if (present(o3_surf)) then
             profiles(irun1)% s2m% o = o3_surf(irun1)
          else
             profiles(irun1)% s2m% o = default_o3_surf
          endif

          if (present(wfetc)) then
             profiles(irun1)% s2m% wfetc = wfetc(irun1)
          else
             profiles(irun1)% s2m% wfetc = default_wfetch
          endif

       enddo


       if (present(o3)) then
         do irun1 = 1, nprof
            profiles(irun1)% o3(1) = o3(1,irun1)
         enddo
!cdir novector
         do i=1,n_levs
!cdir vector
           do irun1 = 1, nprof
              profiles(irun1)% o3(i+nlevs_top) = o3(i,irun1)
           enddo
         enddo
       endif

       if (present(co2)) then
         do irun1 = 1, nprof
            profiles(irun1)% co2(1) = co2(1,irun1)
         enddo
!cdir novector
         do i=1,n_levs
!cdir vector
           do irun1 = 1, nprof
              profiles(irun1)% co2(i + nlevs_top) = co2(i,irun1)
           enddo
         enddo
       endif

       if (present(n2o)) then
         do irun1 = 1, nprof
            profiles(irun1)% n2o(1) = n2o(1,irun1)
         enddo
!cdir novector
         do i=1,n_levs
!cdir vector
           do irun1 = 1, nprof
              profiles(irun1)% n2o(i + nlevs_top) = n2o(i,irun1)
           enddo
         enddo
       endif

       if (present(co)) then
         do irun1 = 1, nprof
            profiles(irun1)% co(1) = co(1,irun1)
         enddo
!cdir novector
         do i=1,n_levs
!cdir vector
           do irun1 = 1, nprof
              profiles(irun1)% co(i + nlevs_top) = co(i,irun1)
           enddo
         enddo
       endif

       if (present(ch4)) then
         do irun1 = 1, nprof
            profiles(irun1)% ch4(1) = ch4(1,irun1)
         enddo
!cdir novector
         do i=1,n_levs
!cdir vector
           do irun1 = 1, nprof
              profiles(irun1)% ch4(i + nlevs_top) = ch4(i,irun1)
           enddo
         enddo
       endif

       if (present(clw)) then
         do irun1 = 1, nprof
            profiles(irun1)% clw(1) = clw(1,irun1)
         enddo
!cdir novector
         do i=1,n_levs
!cdir vector
           do irun1 = 1, nprof
              profiles(irun1)% clw(i + nlevs_top) = clw(i,irun1)
           enddo
         enddo
       endif

       if (present(aerosols)) then
!cdir novector
          do i=1,n_lev_layers
!cdir novector
             do j=1,size(aerosols,1)
!cdir vector
                do irun1=1,nprof
                   profiles(irun1)% aerosols(j,i) = aerosols(j,i,irun1)
                enddo
             enddo
          enddo
       endif

       if (associated(cld_p)) then
!cdir novector
          do i=1,n_lev_layers
!cdir novector
             do j=1,size(cld_p,1)
!cdir vector
                do irun1=1,nprof
                   profiles(irun1)% cloud(j,i) = cld_p(j,i,irun1)
                enddo
             enddo
          enddo
       endif

       if (associated(cfr_p)) then
#if defined(RTTOV12)
!cdir novector
          do i=1,n_lev_layers
             do irun1=1,nprof
                profiles(irun1)% cfrac(i) = cfr_p(i,irun1)
             enddo
          enddo
#else
!cdir novector
          do j=1,size(profiles(irun1)% cfrac,1)
!cdir novector
             do i=1,n_lev_layers
                do irun1=1,nprof
                   profiles(irun1)% cfrac(j,i) = cfr_p(i,irun1)
                enddo
             enddo
          enddo
#endif
       endif
    else ! Vectorize levels
       do irun1=1,nprof
          if (associated(id_p)) then
            write(profiles(irun1)% id,*) id_p(irun1)
          end if
          !  profile variables:
          if (.not.addinterp_p .and. n_levs_top==1) then
             profiles(irun1)% p(1)            = 0.005_jprb
             profiles(irun1)% t(1)            = min(tmax_ifc,max(tmin_ifc,temp_p(1,irun1)))
             profiles(irun1)% q(1)            = min(qmax_ifc,max(qmin_ifc,humi_p(1,irun1)))
          end if
          if (associated(press_p)) then
            profiles(irun1)% p(1+nlevs_top:) = press_p(1:n_levs,irun1)
          else
            profiles(irun1)% p(1+nlevs_top:) = press_const_p(1:n_levs)
          end if
          profiles(irun1)% t(1+nlevs_top:) = min(tmax_ifc,max(tmin_ifc,temp_p(1:n_levs,irun1)))
          do i=1,n_levs
            profiles(irun1)% q(i+nlevs_top) = min(qmax_ifc,max(qmin_ifc,humi_p(i,irun1)))
          enddo

          !  2 meter air variables:
          profiles(irun1)% s2m% t          = min(tmax_ifc,max(tmin_ifc,t2m_p(irun1)))
          profiles(irun1)% s2m% q          = min(qmax_ifc,max(qmin_ifc,q2m_p(irun1)))
          profiles(irun1)% s2m% p          = psurf_p     (irun1)
          profiles(irun1)% s2m% u          = u10m_p      (irun1)
          profiles(irun1)% s2m% v          = v10m_p      (irun1)

          !  skin variables
          profiles(irun1)% skin% t         = stemp_p     (irun1)
          profiles(irun1)% skin% surftype  = int(stype_p (irun1),jpim)
          profiles(irun1)% skin% salinity  = default_salinity

          ! cloud top pressure and fraction
          if (associated(ctp_p)) then
            profiles(irun1)% ctp           = ctp_p       (irun1)
          else
            profiles(irun1)% ctp           = default_ctp
          end if
          if (associated(clf_p)) then
            profiles(irun1)% cfraction     = clf_p       (irun1)
          else
            profiles(irun1)% cfraction     = default_cfraction
          end if

          !  other variables from interactive inputs
          profiles(irun1)% elevation       = hsurf_p     (irun1)
          profiles(irun1)% zenangle        = abs(satzenith_p (irun1))
          profiles(irun1)% latitude        = lat_p(irun1)
          if (associated(lon_p)) then
            profiles(irun1)% longitude     = lon_p(irun1)
          end if

          if (present(idg)) then
            profiles(irun1)% idg           = int(idg(irun1),jpim)
          else
            profiles(irun1)% idg           = default_idg
          end if
#if defined(RTTOV12)
          if (present(ice_scheme)) then
            profiles(irun1)%ice_scheme     = int(ice_scheme(irun1),jpim)
          else
            profiles(irun1)% ice_scheme    = default_ice_scheme
          end if
#else
          if (present(ice_scheme)) then
            profiles(irun1)%ish            = int(ice_scheme(irun1),jpim)
          else
            profiles(irun1)% ish           = default_ice_scheme
          end if
#endif

          if (present(watertype)) then
            profiles(irun1)%skin%watertype   = int(watertype (irun1),jpim)
          else
            profiles(irun1)% skin% watertype = default_watertype
          end if
          if (present(fastem)) then
            profiles(irun1)%skin%fastem   = fastem
          else
            profiles(irun1)% skin% fastem    = default_fastem
          end if
          if (associated(satazim_p)) then
            profiles(irun1)%azangle     = satazim_p      (irun1)
          else
            profiles(irun1)% azangle         = default_satazim
          end if
          if (associated(sunazangle_p)) then
            profiles(irun1)%sunazangle  = sunazangle_p   (irun1)
          else
            profiles(irun1)% sunazangle      = default_sunazangle
          end if
          if (associated(sunZenith_p)) then
            profiles(irun1)%sunzenangle = sunZenith_p    (irun1)
          else
            profiles(irun1)%sunzenangle = default_sunzenangle
          end if
          if (present(o3_surf)) then
            profiles(irun1)%s2m% o        = o3_surf       (irun1)
          else
            profiles(irun1)% s2m% o          = default_o3_surf
          end if
          if (present(wfetc)) then
            profiles(irun1)%s2m% wfetc    = wfetc         (irun1)
          else
            profiles(irun1)% s2m% wfetc      = default_wfetch
          end if
          if (associated(cld_p))        profiles(irun1)% cloud (1:size(cld_p,1),1:n_lev_layers )=&
               cld_p(1:size(cld_p,1),1:n_lev_layers,irun1)
#if defined(RTTOV12)
          if (associated(cfr_p))        profiles(irun1)% cfrac (1:n_lev_layers )=cfr_p(1:n_lev_layers,irun1)
#else
          if (associated(cfr_p)) then
             do j = 1, size(profiles(irun1)% cfrac, 1)
                profiles(irun1)% cfrac (j,1:n_lev_layers ) = cfr_p(1:n_lev_layers,irun1)
             end do
          end if
#endif
          if (present(aerosols))   profiles(irun1)% aerosols(1:size(aerosols,1),1:n_lev_layers )=&
               aerosols(1:size(aerosols,1),1:n_lev_layers,irun1)

          if (present(o3)) then
             profiles(irun1)% o3(1)               = o3  (1       ,irun1)
             profiles(irun1)% o3(1+nlevs_top:)    = o3  (1:n_levs,irun1)
          endif
          if (present(co2)) then
             profiles(irun1)% co2(1)               = co2 (1       ,irun1)
             profiles(irun1)% co2(1+nlevs_top:)    = co2 (1:n_levs,irun1)
          endif
          if (present(n2o)) then
             profiles(irun1)% n2o(1)               = n2o  (1       ,irun1)
             profiles(irun1)% n2o(1+nlevs_top:)    = n2o  (1:n_levs,irun1)
          endif
          if (present(co)) then
             profiles(irun1)% co(1)                = co  (1       ,irun1)
             profiles(irun1)% co(1+nlevs_top:)     = co  (1:n_levs,irun1)
          endif
          if (present(ch4)) then
             profiles(irun1)% ch4(1)               = ch4  (1       ,irun1)
             profiles(irun1)% ch4(1+nlevs_top:)    = ch4  (1:n_levs,irun1)
          endif
          if (present(clw)) then
             profiles(irun1)% clw(1)               = clw  (1       ,irun1)
             profiles(irun1)% clw(1+nlevs_top:)    = clw  (1:n_levs,irun1)
          endif

       enddo

    endif

    ! Missing t_surf/T_G in input to var3d/mec results in skin%t == 0. This will crash
    ! in RTTOV if the do_checkinput option (apply_reg_lims option in var3d/mec) is not set.
    ! Thus, we check here for invalid skin temp. in order to get a useful error message.
    if (any(profiles(1:nprof)% skin% t < tmin_ifc)) then
      rttov_fill_input = ERR_INVALID_TSKIN
      return
    end if


#if defined(RTTOV12)
    if (present(wr_profs)) then
      if (size(wr_profs) > 0) then
        do irun1 = 1, nprof
          if (any(id_p(irun1) == wr_profs(:))) then
            write(fname,trim(wr_profs_fmt)) id_p(irun1)
            call rttov_write_profile(profiles(irun1:irun1), opts, trim(fname), stat)
            if (stat /= 0) then
              rttov_fill_input = ERR_WR_PROF
              return
            end if
          end if
        end do
      end if
    end if
#endif

FTRACE_END('rttov_fill_input')

#endif

  end function rttov_fill_input
  !======================================================================


#if defined(RTTOV12)
  subroutine rttov_write_profile(prof, opts, name, stat)
    use hdf5
    use rttov_hdf_mod
    type(rttov_profile), intent(in)  :: prof(:)
    type(rttov_options), intent(in)  :: opts
    character(len=*),    intent(in)  :: name
    integer,             intent(out) :: stat

    call open_hdf(.true., stat)
    if (stat /= 0) return
    call rttov_hdf_save(stat, trim(name), '/PROFILES', create=.true., profiles=prof(:))
    if (stat /= 0) return
    call rttov_hdf_save(stat, trim(name), '/OPTIONS', create=.false., options=opts)
    if (stat /= 0) return
    call close_hdf(stat)
    if (stat /= 0) return

  end subroutine rttov_write_profile
#endif

  function rttov_direct_ifc (instrIdx,lprofs,chans,emissiv,t_b,satazim,satzenith,t_b_clear,rad, &
                             radclear,radupclear,raddnclear,refdnclear,radovercast,radtotal,    &
                             transm,transmtotal,opdep,height,istore,errorstatus,reg_lim,rflag,  &
                             dealloc,iprint,rad_out_flg,pe,l_pio, l_opdep)
  integer, intent(in)          :: instrIdx           ! RTTOV instrument index,
                                                     ! as in the initialization routine
  integer, intent(in)          :: lprofs         (:) ! list of profile indices (nchans*nprof)
  integer, intent(in)          :: chans          (:) ! list of channel indices (nchans*nprof)
                                                     ! if the DWD optimized version of
                                                     ! RTTOV9.3 is used, the results are only
                                                     ! correct if profiles are the outer loop
                                                     ! and channels the inner loop,
                                                     ! e.g.:  profile channel
                                                     !        1       1
                                                     !        1       2
                                                     !        1       3
                                                     !        2       1
                                                     !        2       2
  real(wp),intent(inout)       :: emissiv      (:,:) ! emissivities (nchans,nprof),
                                                     ! if < 0.01, they will be calculated,
                                                     ! else they will be used
  real(wp),intent(out)         :: t_b          (:,:) ! brightness temp. (nchans,nprof) [K]
  real(wp),intent(in) ,optional:: satazim        (:) ! local satellite azimuth angle (nprofs)
                                                     ! [deg] range: 0-360; east=90

  real(wp),intent(in) ,optional:: satzenith      (:) ! local satellite zenith angle (nprofs)
                                                     ! [deg] range: 0-360; east=90;
  real(wp),intent(out),optional:: t_b_clear    (:,:) ! calc. clear sky b.t. (nchans,nprof) [K]
  real(wp),intent(out),optional:: rad          (:,:) ! calculated total radiance
                                                     ! (nchans,nprof) [mW/cm^-1/ster/m^2]
  real(wp),intent(out),optional:: radclear     (:,:) ! calculated clear sky radiances
                                                     ! (nchans,nprof) [mW/cm^-1/ster/m^2]
  real(wp),intent(out),optional:: radupclear   (:,:) ! calculated clear sky upwelling radiances
                                                     ! includes surface emission term but not downwelling reflected
                                                     ! (nchans,nprof) [mW/cm^-1/ster/m^2]
  real(wp),intent(out),optional:: raddnclear   (:,:) ! calculated clear sky downwelling radiances at surface
                                                     ! (nchans,nprof) [mW/cm^-1/ster/m^2]
  real(wp),intent(out),optional:: refdnclear   (:,:) ! reflected clear sky downwelling radiances at TOA
                                                     ! (nchans,nprof) [mW/cm^-1/ster/m^2]
  real(wp),intent(out),optional:: radovercast(:,:,:) ! overcast radiances (nlevs,nchans,nprof)
                                                     ! [mW/cm^-1/ster/m^2]
  real(wp),intent(out),optional:: radtotal     (:,:) ! total radiances (nchans,nprof)
                                                     ! [mW/cm^-1/ster/m^2]
  real(wp),intent(out),optional:: transm     (:,:,:) ! transmission (nlevs,nchans,nprof)
  real(wp),intent(out),optional:: transmtotal  (:,:) ! total surface to TOA transmission (nchans,nprof)

  real(wp),intent(out),optional:: opdep      (:,:,:) ! original optical depth (nlevs,nchans,nprof)

  real(wp),intent(out),optional:: height       (:,:) ! height estimated on basis of weighting function
  integer, intent(in) ,optional:: istore       (:,:) ! put result i into array position (istore(i,1),istore(i,2))
  integer, intent(out),optional:: errorstatus    (:) ! errorstatus (nprof)
  integer, intent(out),optional:: reg_lim    (:,:,:) ! result of apply_reg_limits (nlevs,nvars,nprof)
  integer, intent(out),optional:: rflag        (:,:) ! results of other test like e.g. the chk_god test
  integer, intent(in) ,optional:: iprint             ! profile to printed in RTTOV (debug)
  logical, intent(in) ,optional:: dealloc
  integer, intent(in) ,optional:: rad_out_flg        ! bit field (see OUT_* parameters)
  integer, intent(in) ,optional:: pe
  logical, intent(in) ,optional:: l_pio
  logical, intent(in) ,optional:: l_opdep

  integer :: rttov_direct_ifc                        ! function result
  integer :: ipr
  integer :: rad_out
  !-----------------------------------------------------------------------
  ! this function renews the instrument dependent part
  ! of a call of rttov_direct in the case that there are more
  ! than one instrument looking at the same ground point
  ! from the same satellite (or at least if their measurements
  ! are interpolated to the same ground point by e.g., aapp)
  !-----------------------------------------------------------------------
#ifndef NO_RTTOV
   integer(jpim)              :: irun1
   integer(jpim)              :: nchans
   integer(jpim)              :: nchansprofs
   integer(jpim)              :: nprof_store
   integer(jpim)              :: nprof_calc
   integer(jpim)              :: iprof_s
   integer(jpim)              :: iprof_e
   integer(jpim)              :: sind
   integer(jpim)              :: eind
   integer(jpim)              :: i, k, minchan, maxchan, minprof, maxprof
   integer(jpim)              :: stat
   integer(jpim)              :: lprofs_aux    (size(lprofs))
   integer(jpim), allocatable :: index_prof    (:)
   logical(jplm), allocatable :: prof_used     (:)
   logical                    :: l_dealloc
   logical(jplm)              :: calcemis      (size(chans))
#if defined(RTTOV12)
   type(rttov_emissivity)     :: emissivity    (size(chans))
#else
   real(jprb)                 :: emissivity_out(size(chans))
#if defined(RTTOV10)
   real(jprb)                 :: emissivity    (size(chans))
#endif
#endif
   integer(jpim)              :: rttov_errorstatus(1)
   type(rttov_chanprof)       :: chanprof(size(chans))
   integer                    :: lims_flag(nlevs,2)
#endif

#ifndef NO_RTTOV
   integer :: pe_loc
   logical :: lpio
   logical :: lopdep, l_transm

   type(rttov_profile), pointer :: p => null()
   character(len=1000) :: msg  = ''
   character(len=1000) :: msg_ = ''
#endif

   rttov_direct_ifc = ini_rttov_ifc_errHandle()
   if (rttov_direct_ifc /= NO_ERROR) return
#ifdef NO_RTTOV
   rttov_direct_ifc      = ERR_NO_RTTOV_LIB
   return
#else

   lpio = .false.
   if (present(l_pio)) lpio=l_pio

    lopdep = .false.
    if (present(l_opdep) .and. present(opdep)) lopdep=l_opdep
    l_transm = .false.
    if (present(transm)) then
      l_transm = all(shape(transm) > 0)
    end if

FTRACE_BEGIN('rttov_direct_ifc')

    if (present(iprint)) then
       ipr = iprint
    else
       ipr = 0
    end if
    if (ipr > 0) then
      ipr_deb = ipr
      call print_profile(ipr)
    end if

    if (present(pe)) then
      pe_loc = pe
    else
      pe_loc = -1
    end if
    pe_rt = pe_loc

    if (present(rad_out_flg)) then
      rad_out = rad_out_flg
    else
      rad_out = 0
      rad_out = ibset(rad_out,OUT_ASB)
      rad_out = ibset(rad_out,OUT_CSB)
      rad_out = ibset(rad_out,OUT_ASR)
      rad_out = ibset(rad_out,OUT_CSR)
    end if

    l_dealloc = .true. ; if (present(dealloc)) l_dealloc = dealloc

    nchansprofs    = size(chans)
    if (present(istore)) then
       nchans      = maxval(istore(1:nchansprofs,1))
       nprof_store = maxval(istore(1:nchansprofs,2))
    else
       nchans      = nchansprofs/nprof
       nprof_store = nprof
    end if


    ! RTTOV(10) crashes, if the number of profiles supplied to it is bigger than
    ! the number of forward-calculations. This can happen if e.g.
    ! chans = (/ 1, 2, 1, 2 /) and lprofs = (/ 1, 1, 7, 7 /)  (7 profiles and 4
    ! forward calculations. Therefore the profiles array has to be packed in the
    ! RTTOV call: pack(profiles(iprof_s:iprof_e), mask=prof_used)
    ! In the following an appropriate lprofs-array "lprofs_aux" is derived.
    ! 1. Determine the used profiles
    iprof_s     = minval(lprofs(:))
    iprof_e     = maxval(lprofs(:))
    allocate(prof_used(iprof_s:iprof_e))
    prof_used(:) = .false.
    do i = 1, nchansprofs
      prof_used(lprofs(i)) = .true.
    end do
    nprof_calc  = count(prof_used(:))
    if (nchansprofs < nprof_calc) then
      write(0,*) 'Less RTTOV-calculations than profiles! RTTOV will crash'
      write(0,*) 'nchansprofs',nchansprofs
      write(0,*) 'nprof_calc',nprof_calc
      write(0,*) 'iprof_s, iprof_e', iprof_s, iprof_e
      write(0,*) 'lprofs', lprofs(:)
      stop
    end if
    ! 2. Determine an array to translate lprofs into lprofs_aux
    allocate(index_prof(iprof_s:iprof_e))
    k = 1
    do i = iprof_s, iprof_e
      if (prof_used(i)) then
        index_prof(i) = k
        k = k + 1
        if (i==ipr_deb) ipr_deb = index_prof(i)
      end if
    end do
    ! 3. Translate lprofs into lprofs_aux
    do i = 1, nchansprofs
      lprofs_aux(i) = index_prof(lprofs(i))
    end do
    deallocate(index_prof)

    if (.not.present(istore)) then
      sind = 1
      eind = nchans
      do irun1=1,nprof
#if defined(RTTOV12)
        emissivity(sind:eind)%emis_in  = emissiv(1:nchans,irun1)
        emissivity(sind:eind)%emis_out = 0._jprb
#elif defined(RTTOV10)
        emissivity(sind:eind)         = emissiv(1:nchans,irun1)
#endif
        sind = eind + 1
        eind = eind + nchans
      enddo
    else
!cdir nodep
      do i = 1, nchansprofs
#if defined(RTTOV12)
        emissivity(i)%emis_in  = emissiv(istore(i,1),istore(i,2))
        emissivity(i)%emis_out = 0._jprb
#elif defined(RTTOV10)
        emissivity(i)         = emissiv(istore(i,1),istore(i,2))
#endif
      end do
    end if

#define DIM_ERROR(text) rttov_direct_ifc = ERR_DIM ; print*,text ; return

    if (.not.present(istore)) then
      ! check dimensions of input and output arrays:
      if (size(lprofs ) /= nchansprofs) then
        DIM_ERROR('lprofs')
      endif
      if (size(emissiv) /= nchansprofs) then
        DIM_ERROR('emissiv')
      endif
      if ((size(t_b,1) /= nchans)  .or. &
           (size(t_b,2) /= nprof_store )) then
        DIM_ERROR('t_b')
      endif
      if (btest(rad_out,OUT_CSB)) then
        if (present(t_b_clear)) then
          if ((size(t_b_clear,1) /= nchans) .or. &
               (size(t_b_clear,2) /= nprof_store )) then
            DIM_ERROR('t_b_clear')
          endif
        endif
      end if
      if (btest(rad_out,OUT_ASR)) then
        if (present(rad   )) then
          if ((size(rad,1) /= nchans ) .or. &
               (size(rad,2) /= nprof_store  )) then
            DIM_ERROR('rad')
          endif
        endif
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(radclear)) then
          if ((size(radclear,1) /= nchans) .or. &
               (size(radclear,2) /= nprof_store )) then
            DIM_ERROR('radclear')
          endif
        endif
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(radupclear)) then
          if ((size(radupclear,1) /= nchans) .or. &
               (size(radupclear,2) /= nprof_store )) then
            DIM_ERROR('radupclear')
          endif
        endif
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(raddnclear)) then
          if ((size(raddnclear,1) /= nchans) .or. &
               (size(raddnclear,2) /= nprof_store )) then
            DIM_ERROR('raddnclear')
          endif
        endif
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(refdnclear)) then
          if ((size(refdnclear,1) /= nchans) .or. &
               (size(refdnclear,2) /= nprof_store )) then
            DIM_ERROR('refdnclear')
          endif
        endif
      end if
      if (present(radovercast)) then
        if ((size(radovercast,1) /= nlevs )  .or. &
             (size(radovercast,2) /= nchans)  .or. &
             (size(radovercast,3) /= nprof_store )) then
          DIM_ERROR('radovercast')
        endif
      endif
      if (present(radtotal)) then
        if ((size(radtotal,1) /= nchans)  .or. &
             (size(radtotal,2) /= nprof_store )) then
          DIM_ERROR('radtotal')
        endif
      endif
      if (present(transmtotal)) then
        if ((size(transmtotal,1) /= nchans)  .or. &
             (size(transmtotal,2) /= nprof_store )) then
          DIM_ERROR('transmtotal')
        endif
      endif
      if (present(height)) then
        if ((size(height,1) /= nchans)  .or. &
             (size(height,2) /= nprof_store )) then
          DIM_ERROR('height')
        endif
      endif
    else
      minchan = minval(istore(1:nchansprofs,1))
      maxchan = maxval(istore(1:nchansprofs,1))
      minprof = minval(istore(1:nchansprofs,2))
      maxprof = maxval(istore(1:nchansprofs,2))
      if ((minchan < lbound(t_b, 1)) .or. (maxchan > ubound(t_b, 1))) then
        DIM_ERROR('t_b(1)')
      endif
      if ((minprof < lbound(t_b, 2)) .or. (maxprof > ubound(t_b, 2))) then
        DIM_ERROR('t_b(2)')
      endif
      if (btest(rad_out,OUT_CSB)) then
        if (present(t_b_clear)) then
          if ((minchan < lbound(t_b_clear, 1)) .or. (maxchan > ubound(t_b_clear, 1))) then
            DIM_ERROR('t_b_clear')
          endif
        end if
      end if
      if (btest(rad_out,OUT_ASR)) then
        if (present(rad)) then
          if ((minchan < lbound(rad, 1)) .or. (maxchan > ubound(rad, 1))) then
            DIM_ERROR('rad')
          endif
        end if
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(radclear)) then
          if ((minchan < lbound(radclear, 1)) .or. (maxchan > ubound(radclear, 1))) then
            DIM_ERROR('radclear')
          endif
        end if
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(radupclear)) then
          if ((minchan < lbound(radupclear, 1)) .or. (maxchan > ubound(radupclear, 1))) then
            DIM_ERROR('radupclear')
          endif
        end if
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(raddnclear)) then
          if ((minchan < lbound(raddnclear, 1)) .or. (maxchan > ubound(raddnclear, 1))) then
            DIM_ERROR('raddnclear')
          endif
        end if
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(refdnclear)) then
          if ((minchan < lbound(refdnclear, 1)) .or. (maxchan > ubound(refdnclear, 1))) then
            DIM_ERROR('refdnclear')
          endif
        end if
      end if
      if (present(radovercast)) then
        if ((minchan < lbound(radovercast, 1)) .or. (maxchan > ubound(radovercast, 1))) then
          DIM_ERROR('radovercast')
        endif
      end if
      if (present(radtotal)) then
        if ((minchan < lbound(radtotal, 1)) .or. (maxchan > ubound(radtotal, 1))) then
          DIM_ERROR('radtotal')
        endif
      end if
      if (present(transmtotal)) then
        if ((minchan < lbound(transmtotal, 1)) .or. (maxchan > ubound(transmtotal, 1))) then
          DIM_ERROR('transmtotal')
        endif
      end if
      if (present(height)) then
        if ((minchan < lbound(height, 1)) .or. (maxchan > ubound(height, 1))) then
          DIM_ERROR('height')
        endif
      end if
    endif

    if (present(satazim)) then
       if (size(satazim(:)) < nprof) then
         DIM_ERROR('satazim')
       endif
       do irun1=1,nprof
          profiles(irun1)% azangle  = satazim(irun1)
       enddo
    endif
    if (present(satzenith)) then
       if (size(satzenith(:)) < nprof) then
         DIM_ERROR('satzenith')
       endif
       do irun1=1,nprof
          profiles(irun1)% zenangle = satzenith(irun1)
       enddo
    endif
    if (chk_reg_lims /= 0) then
      if (present(reg_lim)) then
        if (size(reg_lim,1) < nlevs       .or. &
             size(reg_lim,2) < 2           .or. &
             size(reg_lim,3) < nprof_store) then
          DIM_ERROR('reg_lim')
        endif
      else
        DIM_ERROR('reg_lim not present')
      endif
    end if

#undef DIM_ERROR

#if defined(RTTOV12)
#ifndef RADSHARE
    transmission%l_opdep = lopdep !TEST
#endif
    stat = realloc_rttov_arrays(                             &
                                  n_profiles = nchansprofs,      &
                                  n_levels   = nlevs + nlevs_top,&
                                  n_channels = nchansprofs,      &
                                  rads       = radiance,         &
                                  rads2      = radiance2,        &
                                  transm     = transmission)
#else

    stat = realloc_rttov_arrays(                             &
                                  n_profiles = nchansprofs,      &
                                  n_levels   = nlevs + nlevs_top,&
                                  n_channels = nchansprofs,      &
                                  rads       = radiance,         &
                                  transm     = transmission)
#endif
    if (stat /= NO_ERROR) then
      rttov_direct_ifc = int(stat)
      errorstatus(:)   = int(stat)
      return
    endif

    ! clear variable output part of rttov structures
    call rttov_clear_rad_var(nchansprofs,radiance)
#if defined(RTTOV12)
    call rttov_clear_rad2_var(nchansprofs,radiance2)
    call rttov_init_transmission(transmission)
#else
    call rttov_init_transmission(nchansprofs, transmission)
#endif

    addclouds_p = instr_addclouds(instrIdx)
    addaerosl_p = instr_addaerosl(instrIdx)

    ! fill rttov simulated obs. specific input
#if defined(RTTOV12)
    calcemis(1:nchansprofs) = emissivity    (1:nchansprofs)%emis_in < 0.01_jprb
#elif defined(RTTOV10)
    calcemis(1:nchansprofs) = emissivity    (1:nchansprofs) < 0.01_jprb
#endif

    rttov_errorstatus(:) = errorstatus_success

    !---------------------
    ! fill in the chanprof
    !---------------------
    do irun1=1,nchansprofs
      chanprof(irun1)% chan = chans(irun1)
      chanprof(irun1)% prof = lprofs_aux(irun1)
    enddo
    !-----------
    ! call rttov
    !-----------
#if defined(RTTOV12)
#ifdef RTTOV_USE_OPENMP
    call rttov_parallel_direct (                            &
           rttov_errorstatus(1)                            ,& ! --> error flag
           chanprof(1:nchansprofs)                         ,& ! <-- channels and profiles to calculate
           opts                                            ,& ! <-- options
           pack(profiles(iprof_s:iprof_e), mask=prof_used) ,& ! <--  profile array
           coefs(instrIdx)                                 ,& ! <--  coefs array
           transmission                                    ,& ! --> array of transmittances
           radiance                                        ,& ! <-> computed radiance array
           radiance2                                       ,& ! <--> computed secondary radiance array
           calcemis=calcemis(1:nchansprofs)                ,& ! <-- flag for internal emissivity calc
           emissivity=emissivity(1:nchansprofs)            ,& ! <-- input emissivities per channel
           nthreads=omp_get_max_threads()                   )
#else
    call rttov_direct (                                     &
           rttov_errorstatus(1)                            ,& ! --> error flag
           chanprof(1:nchansprofs)                         ,& ! <-- channels and profiles to calculate
           opts                                            ,& ! <-- options
           pack(profiles(iprof_s:iprof_e), mask=prof_used) ,& ! <--  profile array
           coefs(instrIdx)                                 ,& ! <--  coefs array
           transmission                                    ,& ! --> array of transmittances
           radiance                                        ,& ! <-> computed radiance array
           radiance2                                       ,& ! <--> computed secondary radiance array
           calcemis=calcemis(1:nchansprofs)                ,& ! <-- flag for internal emissivity calc
           emissivity=emissivity(1:nchansprofs)            )  ! <-- input emissivities per channel
#endif
#elif defined(RTTOV10)
    call rttov_direct (                                     &
           rttov_errorstatus(1)                            ,& ! --> error flag
           chanprof(1:nchansprofs)                         ,& ! <-- channels and profiles to calculate
           opts                                            ,& ! <-- options
           pack(profiles(iprof_s:iprof_e), mask=prof_used) ,& ! <--  profile array
           coefs(instrIdx)                                 ,& ! <--  coefs array
           calcemis(1:nchansprofs)                         ,& ! <-- flag for internal emissivity calc
           emissivity(1:nchansprofs)                       ,& ! <-- input emissivities per channel
           emissivity_out(1:nchansprofs)                   ,& ! --> emissivities used by rttov
           transmission                                    ,& ! --> array of transmittances
           radiance                                        )  ! <-> computed radiance array
#endif
    !----------------------------------------------
    ! check for completely awful rttov error status
    !----------------------------------------------
    if (rttov_errorstatus(1) == ERRORSTATUS_FATAL) then
       print *, 'after rttov direct -> fatal error'
       rttov_direct_ifc = ERROR_RTTOV_CALL
       return
    endif

    if (present(height)) then
#ifndef RADSHARE
      stat = realloc_rttov_arrays(                           &
                                  n_profiles = nchansprofs,      &
                                  n_levels   = nlevs + nlevs_top,&
                                  n_channels = nchansprofs,      &
                                  height     = height_aux)
      if (stat /= NO_ERROR) then
        rttov_direct_ifc = int(stat)
        errorstatus(:)   = int(stat)
        return
      endif

      call get_height(transmission, profiles(iprof_s:iprof_e), chanprof(1:nchansprofs), &
           height_aux, stat, pe=pe)
      if (stat /= NO_ERROR) then
        rttov_direct_ifc = int(stat)
        errorstatus(:)   = int(stat)
        return
      endif
#else
      rttov_direct_ifc = ERR_INPUT
      return
#endif
    end if

    ! fill result into output arrays
    if (.not.present(istore)) then
      sind = 1
      eind = nchans
      do irun1=1,nprof
#if defined(RTTOV12)
        emissiv(1:nchans,irun1) = dble(emissivity(sind:eind)%emis_out)
#else
        emissiv(1:nchans,irun1) = dble(emissivity_out(sind:eind))
#endif
        t_b(1:nchans,irun1) = dble(radiance% bt(sind:eind))
        if (btest(rad_out,OUT_CSB)) then
          if (present(t_b_clear)) t_b_clear(1:nchans,irun1)=dble(radiance%bt_clear(sind:eind))
        end if
        if (btest(rad_out,OUT_ASR)) then
          if (present(rad      )) rad      (1:nchans,irun1)=dble(radiance%total   (sind:eind))
        end if
        if (btest(rad_out,OUT_CSR)) then
          if (present(radclear )) radclear (1:nchans,irun1)=dble(radiance%clear   (sind:eind))
        end if
        if (btest(rad_out,OUT_CSR)) then
#if defined(RTTOV12)
          if (present(radupclear )) radupclear (1:nchans,irun1)=dble(radiance2%upclear   (sind:eind))
#else
          if (present(radupclear )) radupclear (1:nchans,irun1)=dble(radiance%upclear   (sind:eind))
#endif
        end if
        if (btest(rad_out,OUT_CSR)) then
#if defined(RTTOV12)
          if (present(raddnclear )) raddnclear (1:nchans,irun1)=dble(radiance2%dnclear   (sind:eind))
          if (present(refdnclear )) refdnclear (1:nchans,irun1)=dble(radiance2%refldnclear   (sind:eind))
#else
          if (present(raddnclear )) raddnclear (1:nchans,irun1)=dble(radiance%dnclear   (sind:eind))
          if (present(refdnclear )) refdnclear (1:nchans,irun1)=dble(radiance%reflclear   (sind:eind))
#endif
        end if
        if (present(radtotal))  radtotal (1:nchans,irun1)=dble(radiance%clear   (sind:eind)) !clear -> total ???
        if (present(radovercast)) &
             radovercast(:,1:nchans,irun1) = dble(radiance% overcast(:,sind:eind))
        if (present(height))    height   (1:nchans,irun1)=height_aux(sind:eind)
        if (l_transm) &
             transm(:,1:nchans,irun1) = dble(transmission%tau_levels(1+nlevs_top:,sind:eind))
        if (present(transmtotal)) transmtotal(1:nchans,irun1) = dble(transmission%tau_total(sind:eind))
#if defined(RTTOV12) && !defined(RADSHARE)
        if (lopdep) &
             opdep (:,1:nchans,irun1) = dble(transmission%opdep_ref (1+nlevs_top:,sind:eind))
#endif
        sind = eind + 1
        eind = eind + nchans
      enddo
    else
!cdir nodep
      do i = 1, nchansprofs
#if defined(RTTOV12)
        emissiv(istore(i,1),istore(i,2)) = dble(emissivity(i)%emis_out)
#else
        emissiv(istore(i,1),istore(i,2)) = dble(emissivity_out(i))
#endif
        t_b(istore(i,1),istore(i,2)) = dble(radiance% bt(i))
        if (btest(rad_out,OUT_CSB)) then
          if (present(t_b_clear  )) t_b_clear  (istore(i,1),istore(i,2))   = dble(radiance% bt_clear(i))
        end if
        if (btest(rad_out,OUT_ASR)) then
          if (present(rad        )) rad        (istore(i,1),istore(i,2))   = dble(radiance% total   (i))
        end if
        if (btest(rad_out,OUT_CSR)) then
          if (present(radclear   )) radclear   (istore(i,1),istore(i,2))   = dble(radiance% clear   (i))
        end if
        if (btest(rad_out,OUT_CSR)) then
#if defined(RTTOV12)
          if (present(radupclear   )) radupclear   (istore(i,1),istore(i,2))   = dble(radiance2% upclear  (i))
#else
          if (present(radupclear   )) radupclear   (istore(i,1),istore(i,2))   = dble(radiance% upclear   (i))
#endif
#if defined(RTTOV12)
          if (present(raddnclear   )) raddnclear   (istore(i,1),istore(i,2))   = dble(radiance2% dnclear  (i))
          if (present(refdnclear   )) refdnclear   (istore(i,1),istore(i,2))   = dble(radiance2% refldnclear  (i))
#else
          if (present(raddnclear   )) raddnclear   (istore(i,1),istore(i,2))   = dble(radiance% dnclear   (i))
          if (present(refdnclear   )) refdnclear   (istore(i,1),istore(i,2))   = dble(radiance% reflclear   (i))
#endif
        end if
        if (present(radtotal   )) radtotal    (istore(i,1),istore(i,2))   = dble(radiance% clear   (i))
        if (present(radovercast)) radovercast (:,istore(i,1),istore(i,2)) = dble(radiance% overcast(:,i))
        if (present(height     )) height      (istore(i,1),istore(i,2))   = height_aux(i)
        if (present(transm     )) transm(:,istore(i,1),istore(i,2)) = dble(transmission%tau_levels(1+nlevs_top:,i))
        if (present(transmtotal)) transmtotal(istore(i,1),istore(i,2)) = dble(transmission%tau_total(i))
#if defined(RTTOV12) && !defined(RADSHARE)
        if (lopdep) &
             opdep(:,istore(i,1),istore(i,2)) = dble(transmission%opdep_ref(1+nlevs_top:,i))
#endif
      end do
    end if

    if (chk_reg_lims /= 0 .and. nprof > 0) then
      stat = realloc_rttov_arrays(nprof_aux,profiles(1)%nlevels,0,profs=profiles_aux)
      if (stat /= NO_ERROR) then
        rttov_direct_ifc = int(stat)
        errorstatus(:)   = int(stat)
        return
      endif
      do irun1 = iprof_s, iprof_e
        if (prof_used(irun1)) then
          call check_prof(irun1, coefs(instrIdx), lims_flag, lpio)
          if (chk_reg_lims > 1 .and. present(reg_lim)) then
            if (.not.present(istore)) then
              do i = 1, nprof
                if (lprofs((i-1)*nchans+1) ==irun1) then
                  reg_lim(:,:,i) = lims_flag(:,:)
                  exit
                end if
              end do
            else
              do i = 1, nchansprofs
                if (lprofs(i)==irun1) then
                  reg_lim(:,:,istore(i,2)) = lims_flag(:,:)
                  exit
                end if
              end do
            end if
          end if
        end if
      end do
      stat = dealloc_rttov_arrays(profs=profiles_aux)
      if (stat /= NO_ERROR) then
        rttov_direct_ifc = int(stat)
        errorstatus(:)   = int(stat)
        return
      endif
    end if

#if defined(RTTOV12)
    if (chk_god /= 0 .and. nprof > 0 .and. present(rflag)) then
#ifndef RADSHARE
      do i = 1, nchansprofs
        stat = 0
        if (associated(coefs(instrIdx)%coef%god)) then
          if (lprofs(i) == ipr) then
            print*,'check_god_ifc prof',profiles(ipr)%id,chans(i),coefs(instrIdx)%coef%ff_ori_chn(chans(i))
            print*,'check_god_ifc ntr',chans(i),any(coefs(instrIdx)%coef%god(:,chans(i))%ntr > 0)
          end if
          if (any(coefs(instrIdx)%coef%god(:,chans(i))%ntr > 0)) then
            call check_god_infl(transmission%tau_levels(:,i),         &
                                coefs(instrIdx)%coef%god(:,chans(i)), &
                                stat, msg_)
            if (stat /= 0) then
              if (btest(chk_god, 0)) then
                p => profiles(lprofs(i))
                write(msg,'(2x,A,I5.5," god-smoothing too large in ",A," (",f7.3,"N",f8.3,"E)")') &
                     'profile '//trim(p%id)//' chan ',coefs(instrIdx)%coef%ff_ori_chn(chans(i)),&
                     trim(msg_),p%latitude,p%longitude
                WR
              end if
            end if
          end if
        end if
        if ((btest(chk_god,1) .or. btest(chk_god,2)) .and. all(shape(rflag) > 0)) then
          if (present(istore)) then
            if (lprofs(i) == ipr) then
              print*,'check_god_ifc stat',stat,istore(i,:)
            end if
            rflag(istore(i,1), istore(i,2)) = stat
          else
            ! TODO !!!
          end if
        end if
      end do
#else
      rttov_direct_ifc = ERR_INPUT
      return
#endif
    end if
#endif

#ifdef RTTOV10
    if (any(rttov_errorstatus == errorstatus_warning)) &
       rttov_direct_ifc = WARN_RTTOV_DIR_ANY
    if (all(rttov_errorstatus == errorstatus_warning)) &
       rttov_direct_ifc = WARN_RTTOV_DIR_ALL
#endif

    if (present(errorstatus)) errorstatus(1:nprof) = rttov_errorstatus(1)

FTRACE_END('rttov_direct_ifc')

    if (l_dealloc) then
#if defined(RTTOV12)
      stat = dealloc_rttov_arrays( rads       = radiance,   &
                                   rads2      = radiance2,  &
                                   transm     = transmission)
#else
      stat = dealloc_rttov_arrays( rads       = radiance,   &
                                   transm     = transmission)
#endif
      if (stat /= NO_ERROR) errorstatus = stat
    end if
#endif
   end function rttov_direct_ifc
  !======================================================================


  function rttov_k_ifc (instrIdx,lprofs,chans,emissiv,emissiv_k,temp_k,     &
                        humi_k,t2m_k,q2m_k,stemp_k,t_b,t_b_clear,rad,       &
                        radclear,radtotal,radovercast,transm,opdep,psurf_k,u10m_k,v10m_k,&
                        o3_surf_k,wfetc_k,ctp_k,cfraction_k,clw_k,o3_k,     &
                        co2_k,n2o_k,co_k,ch4_k,istore,reg_lim,rflag,dealloc, &
                        iprint,rad_out_flg,conv_overc,pe,l_pio,l_opdep)
  integer ,intent(in)          :: instrIdx           ! rttov instrument index
                                                     !   as for in the initialization routine
  integer ,intent(in)          :: lprofs         (:) ! list of profile indices(nchans*nprof)
  integer ,intent(in)          :: chans          (:) ! list of channel indices(nchans*nprof)
  real(wp),intent(inout)       :: emissiv      (:,:) ! emissivities - calculated if < 0.01,
                                                     !   else they will be used (nchans*nprof)
  real(wp),intent(inout)       :: emissiv_k    (:,:) ! k-matrix of surf. emiss. (nchans*nprof)
  real(wp),intent(out)         :: temp_k     (:,:,:) ! grad. of temp. (nlevs,nchans,nprof) [K/K]
  real(wp),intent(out)         :: humi_k     (:,:,:) ! grad. of w.v. (nlevs,nchans,nprof) [K/ppmv]
  real(wp),intent(out)         :: t2m_k        (:,:) ! grad. of 2m temp. (nchans,nprof) [K/K]
  real(wp),intent(out)         :: q2m_k        (:,:) ! grad. of 2m hum. (nchans,nprof) [K/ppmv]
  real(wp),intent(out)         :: stemp_k      (:,:) ! grad. of surf. skin t. (nchans,nprof) [K/K]
  real(wp),intent(out)         :: t_b          (:,:) ! calculated b.t. (nchans,nprof) [K]
  real(wp),intent(out),optional:: t_b_clear    (:,:) ! calc. clear sky b.t. (nchans,nprof) [K]
  real(wp),intent(out),optional:: rad          (:,:) ! calculated total radiance
                                                     ! (nchans,nprof) [mW/cm^-1/ster/m^2]
  real(wp),intent(out),optional:: radclear     (:,:) ! calculated clear sky radiances
                                                     ! (nchans,nprof) [mW/cm^-1/ster/m^2]
  real(wp),intent(out),optional:: radtotal     (:,:) ! calculated total radiances
                                                     !   (nchans,nprof) [mW/cm^-1/ster/m^2]
  real(wp),intent(out),optional:: radovercast(:,:,:) ! calculated overcast radiances
                                                     !   (nlevs,nchans,nprof) [mW/cm^-1/ster/m^2]
  real(wp),intent(out),optional:: transm     (:,:,:) ! transmission (nlevs,nchans,nprof)

  real(wp),intent(out),optional:: opdep      (:,:,:) ! transmission (nlevs,nchans,nprof)

  real(wp),intent(out),optional:: psurf_k      (:,:) ! grad. of surf. press. (nchans,nprof) [K/hPa]
  real(wp),intent(out),optional:: u10m_k       (:,:) ! grad. of 10m U wind (nchans,nprof) [K/(m/s)]
  real(wp),intent(out),optional:: v10m_k       (:,:) ! grad. of 10m V wind (nchans,nprof) [K/(m/s)]
  real(wp),intent(out),optional:: o3_surf_k    (:,:) ! k  of surf. o3 conc. (nchans,nprof) [K/ppmv]
  real(wp),intent(out),optional:: wfetc_k      (:,:) ! grad. of wind fetch (nchans,nprof) [K/m]
  real(wp),intent(out),optional:: ctp_k        (:,:) ! grad. w.resp. cloud top pressure   [K/hPa]
  real(wp),intent(out),optional:: cfraction_k  (:,:) ! grad. w.resp. cloud fraction       [K]
  real(wp),intent(out),optional:: clw_k      (:,:,:) ! grad. of cloud liquid water (microwave only)
                                                     !   (nlevs,nchans,nprof) [K/(kg/kg)]
  real(wp),intent(out),optional:: o3_k       (:,:,:) ! grad. of o3  (nlevs,nchans,nprof) [K/ppmv]
  real(wp),intent(out),optional:: co2_k      (:,:,:) ! grad. of co2 (nlevs,nchans,nprof) [K/ppmv]
  real(wp),intent(out),optional:: n2o_k      (:,:,:) ! grad. of n2o (nlevs,nchans,nprof) [K/ppmv]
  real(wp),intent(out),optional:: co_k       (:,:,:) ! grad. of co  (nlevs,nchans,nprof) [K/ppmv]
  real(wp),intent(out),optional:: ch4_k      (:,:,:) ! grad. of ch4 (nlevs,nchans,nprof) [K/ppmv]
  integer, intent(in) ,optional:: istore       (:,:) ! put result i into array position (istore(i,1),istore(i,2))
  integer, intent(out),optional:: reg_lim    (:,:,:) ! result of apply_reg_limits (nlevs,nvars,nprof)
  integer, intent(out),optional:: rflag        (:,:) ! results of other test like e.g. the chk_god test
  logical, intent(in) ,optional:: dealloc
  integer, intent(in) ,optional:: rad_out_flg        ! bit field (see OUT_* parameters)
  integer, intent(in) ,optional:: iprint
  logical, intent(in) ,optional:: conv_overc
  integer, intent(in) ,optional:: pe
  logical, intent(in) ,optional:: l_pio
  logical, intent(in) ,optional:: l_opdep


  integer                      :: rttov_k_ifc        ! function returns an error code or
                                                     !   NO_ERROR (=0) in case of success
  !-----------------------------------------------------------------------
  ! this function renews the instrument dependent part
  ! of a call of rttov_k in the case that there are more
  ! than one instrument looking at the same ground point
  ! from the same satellite (or at least if their measurements
  ! are interpolated to the same ground point by e.g., aapp)
  !----------------------------------------------------------------------
  integer                      :: rad_out
  integer                      :: pe_loc
  integer                      :: ipr
#ifndef NO_RTTOV
    integer(jpim)              :: irun1
    integer(jpim)              :: irun2
    integer(jpim)              :: i, k
    integer(jpim)              :: stat
    integer(jpim)              :: nchans
    integer(jpim)              :: nchansprofs
    integer(jpim)              :: nprof_store
    integer(jpim)              :: nprof_calc
    integer(jpim)              :: iprof_s
    integer(jpim)              :: iprof_e
    integer(jpim)              :: sind
    integer(jpim)              :: eind
    integer(jpim)              :: count1
    integer(jpim)              :: nusedchans
    integer(jpim)              :: lprofs_aux      (size(lprofs))
    integer(jpim), allocatable :: index_prof      (:)
    logical(jplm), allocatable :: prof_used       (:)
    logical                    :: l_dealloc
    logical(jplm)              :: calcemis        (size(chans ))
#if defined(RTTOV12)
   type(rttov_emissivity)      :: emissivity      (size(chans ))
   type(rttov_emissivity)      :: emissivity_k    (size(chans ))
#else
    real(jprb)                 :: emissivity_out  (size(chans ))
    real(jprb)                 :: emissivity_out_k(size(chans ))
#if defined(RTTOV10)
    real(jprb)                 :: emissivity(size(chans))
    real(jprb)                 :: emissivity_k(size(chans))
#endif
#endif
    integer(jpim)              :: rttov_errorstatus(1)
    type(rttov_chanprof)       :: chanprof(size(chans))
    integer(jpim)              :: n_lev_layers ! This variable is used to switch between
                                               ! RTTOV9 levels and RTTOV10 layers
                                               ! Currently, we just cur ot the surface level
                                               ! as in the italian RTTOV10 modification
    integer                    :: lims_flag(nlevs,2)
#endif

#ifndef NO_RTTOV
    logical :: lpio
    logical :: lopdep, l_transm

    type(rttov_profile), pointer :: p => null()
    character(len=1000) :: msg  = ''
    character(len=1000) :: msg_ = ''
#endif

#ifdef NO_RTTOV
    rttov_k_ifc = ERR_NO_RTTOV_LIB
    return
#else
    rttov_k_ifc = ini_rttov_ifc_errHandle()
    if (rttov_k_ifc /= NO_ERROR) return

    lpio = .false.
    if (present(l_pio)) lpio=l_pio

    lopdep = .false.
    if (present(l_opdep) .and. present(opdep)) lopdep=l_opdep
    l_transm = .false.
    if (present(transm)) then
      l_transm = all(shape(transm) > 0)
    end if

    if (present(pe)) then
      pe_loc = pe
    else
      pe_loc = -1
    end if
    pe_rt = pe_loc

    if (present(rad_out_flg)) then
      rad_out = rad_out_flg
    else
      rad_out = ibset(0,OUT_ASB)
    end if

#ifdef NO_RTTOV
    rttov_k_ifc           = ERR_NO_RTTOV_LIB
    return
#else

FTRACE_BEGIN('rttov_k_ifc')

    n_lev_layers    = profiles(1)% nlayers

#ifndef NO_RTTOV
    ipr=-1
    if (present(iprint)) then
      ipr = iprint
      ipr_deb = ipr
    end if
    if (ipr > 0) call print_profile(ipr)
#endif

    l_dealloc = .true. ; if (present(dealloc)) l_dealloc = dealloc

    nchansprofs  = size(chans)
    if (present(istore)) then
      nchans      = maxval(istore(1:nchansprofs,1))
      nprof_store = maxval(istore(1:nchansprofs,2))
    else
      nchans      = nchansprofs / nprof
      nprof_store = nprof
    end if

    ! RTTOV(10) crashes, if the number of profiles supplied to it is bigger than
    ! the number of forward-calculations. This can happen if e.g.
    ! chans = (/ 1, 2, 1, 2 /) and lprofs = (/ 1, 1, 7, 7 /)  (7 profiles and 4
    ! forward calculations. Therefore the profiles array has to be packed in the
    ! RTTOV call: pack(profiles(iprof_s:iprof_e), mask=prof_used)
    ! In the following an appropriate lprofs-array "lprofs_aux" is derived.
    ! 1. Determine the used profiles
    iprof_s     = minval(lprofs(:))
    iprof_e     = maxval(lprofs(:))
    allocate(prof_used(iprof_s:iprof_e))
    prof_used(:) = .false.
    do i = 1, nchansprofs
      prof_used(lprofs(i)) = .true.
    end do
    nprof_calc  = count(prof_used(:))
    if (nchansprofs < nprof_calc) then
      write(0,*) 'Less RTTOV-calculations than profiles! RTTOV will crash'
      write(0,*) 'nchansprofs',nchansprofs
      write(0,*) 'nprof_calc',nprof_calc
      write(0,*) 'iprof_s, iprof_e', iprof_s, iprof_e
      write(0,*) 'lprofs', lprofs(:)
      stop
    end if
    ! 2. Determine an array to translate lprofs into lprofs_aux
    allocate(index_prof(iprof_s:iprof_e))
    k = 1
    do i = iprof_s, iprof_e
      if (prof_used(i)) then
        index_prof(i) = k
        k = k + 1
        if (i==ipr_deb) ipr_deb = index_prof(i)
      end if
    end do
    ! 3. Translate lprofs into lprofs_aux
    do i = 1, nchansprofs
      lprofs_aux(i) = index_prof(lprofs(i))
    end do
    deallocate(index_prof)

    stat  = NO_ERROR

    if (.not.present(istore)) then
      sind = 1
      eind = nchans
      do irun1=1,nprof
#if defined(RTTOV12)
        emissivity  (sind:eind)%emis_in  = real (emissiv  (1:nchans,irun1),jprb)
        emissivity_k(sind:eind)%emis_in  = 0.0_jprb
        emissivity  (sind:eind)%emis_out = 0.0_jprb
        emissivity_k(sind:eind)%emis_out = 0.0_jprb
#elif defined(RTTOV10)
        emissivity(sind:eind)       = real (emissiv  (1:nchans,irun1),jprb)
        emissivity_k(sind:eind)     = 0.0_jprb
        emissivity_out_k(sind:eind) = 0.0_jprb
#endif
        sind = eind + 1
        eind = eind + nchans
      enddo
    else
!cdir nodep
      do i = 1, nchansprofs
#if defined(RTTOV12)
        emissivity  (i)%emis_in  = real (emissiv  (istore(i,1),istore(i,2)),jprb)
        emissivity_k(i)%emis_in  = 0.0_jprb
        emissivity  (i)%emis_out = 0.0_jprb
        emissivity_k(i)%emis_out = 0.0_jprb
#elif defined(RTTOV10)
        emissivity(i)       = real (emissiv  (istore(i,1),istore(i,2)),jprb)
        emissivity_k(i)     = 0.0_jprb
        emissivity_out_k(i) = 0.0_jprb
#endif
      end do
    end if

#define DIM_ERROR(text) rttov_k_ifc = ERR_DIM ; write(0,*) 'dim_error: '//text ; return
    ! check dimensions of input and output arrays
    if ((size(temp_k,1) /= nlevs) .or. &
        (size(humi_k,1) /= nlevs)) then
      DIM_ERROR('temp_k/humi_k levels')
    endif

    if ((size(temp_k ,2) < nchans) .or. &
        (size(humi_k ,2) < nchans) .or. &
        (size(t2m_k  ,1) < nchans) .or. &
        (size(q2m_k  ,1) < nchans) .or. &
        (size(t_b    ,1) < nchans) .or. &
        (size(stemp_k,1) < nchans)) then
      DIM_ERROR('*_k channels')
    endif

    if ((size(temp_k ,3) <  nprof_store) .or. &
        (size(humi_k ,3) <  nprof_store) .or. &
        (size(t2m_k  ,2) <  nprof_store) .or. &
        (size(q2m_k  ,2) <  nprof_store) .or. &
        (size(t_b    ,2) <  nprof_store) .or. &
        (size(stemp_k,2) <  nprof_store)) then
      DIM_ERROR('*_k,t_b profiles')
    endif

    if ((size(lprofs   ) <  nchansprofs) .or. &
        (size(emissiv  ) <  nchansprofs) .or. &
        (size(emissiv_k) <  nchansprofs)) then
      DIM_ERROR('lprofs,emissiv* nchansprofs')
    endif

    if (present(radtotal)) then
       if ((size(radtotal,1) <  nchans) .or. &
           (size(radtotal,2) <  nprof_store )) then
         DIM_ERROR('radtotal')
       endif
    endif

    if (present(radovercast)) then
       if ((size(radovercast,1) /= n_lev_layers ) .or. &
           (size(radovercast,2) <  nchans) .or. &
           (size(radovercast,3) <  nprof_store )) then
         write(0,*) 'dim_error:',size(radovercast,1),n_lev_layers,&
              size(radovercast,2),nchans,&
              size(radovercast,3),nprof_store
         DIM_ERROR('radovercast')
       endif
    endif

    if (present(psurf_k)) then
       if ((size(psurf_k,1) <  nchans) .or. &
           (size(psurf_k,2) <  nprof_store )) then
         DIM_ERROR('psurf_k')
       endif
    endif

    if (present(u10m_k)) then
       if ((size(u10m_k,1) <  nchans) .or. &
           (size(u10m_k,2) <  nprof_store )) then
         DIM_ERROR('u10m_k')
       endif
    endif

    if (present(v10m_k)) then
       if ((size(v10m_k,1) <  nchans) .or. &
           (size(v10m_k,2) <  nprof_store )) then
         DIM_ERROR('v10m_k')
       endif
    endif

    if (present(o3_surf_k)) then
       if ((size(o3_surf_k,1) <  nchans) .or. &
           (size(o3_surf_k,2) <  nprof_store )) then
         DIM_ERROR('o3_surf_k')
       endif
    endif

    if (present(wfetc_k)) then
       if ((size(wfetc_k,1) <  nchans) .or. &
           (size(wfetc_k,2) <  nprof_store )) then
         DIM_ERROR('wfetc_k')
       endif
    endif

    if (present(ctp_k)) then
       if ((size(ctp_k,1) <  nchans) .or. &
           (size(ctp_k,2) <  nprof_store )) then
         DIM_ERROR('ctp_k')
       endif
    endif

    if (present(cfraction_k)) then
       if ((size(cfraction_k,1) <  nchans) .or. &
           (size(cfraction_k,2) <  nprof_store )) then
         DIM_ERROR('cfraction_k')
       endif
    endif

    if (present(clw_k)) then
       if ((size(clw_k,1) /= nlevs ) .or. &
           (size(clw_k,2) <  nchans) .or. &
           (size(clw_k,3) <  nprof_store )) then
         DIM_ERROR('clw_k')
       endif
    endif

    if (present(o3_k)) then
       if ((size(o3_k,1) /= nlevs ) .or. &
           (size(o3_k,2) <  nchans) .or. &
           (size(o3_k,3) <  nprof_store )) then
         DIM_ERROR('o3_k')
       endif
    endif

    if (present(co2_k)) then
       if ((size(co2_k,1) /= nlevs ) .or. &
           (size(co2_k,2) <  nchans) .or. &
           (size(co2_k,3) <  nprof_store )) then
         DIM_ERROR('co2_k')
       endif
    endif

    if (present(n2o_k)) then
       if ((size(n2o_k,1) /= nlevs ) .or. &
           (size(n2o_k,2) <  nchans) .or. &
           (size(n2o_k,3) <  nprof_store )) then
         DIM_ERROR('n2o_k')
       endif
    endif

    if (present(co_k)) then
       if ((size(co_k,1) /= nlevs ) .or. &
           (size(co_k,2) <  nchans) .or. &
           (size(co_k,3) <  nprof_store )) then
         DIM_ERROR('co_k')
       endif
    endif

    if (present(ch4_k)) then
       if ((size(ch4_k,1) /= nlevs ) .or. &
           (size(ch4_k,2) <  nchans) .or. &
           (size(ch4_k,3) <  nprof_store )) then
         DIM_ERROR('ch4_k')
       endif
    endif

    if (chk_reg_lims /= 0) then
      if (present(reg_lim)) then
        if (size(reg_lim,1) < nlevs       .or. &
            size(reg_lim,2) < 2           .or. &
            size(reg_lim,3) < nprof_store) then
          DIM_ERROR('reg_lim')
        endif
      else
        DIM_ERROR('reg_lim not present')
      endif
    end if

#undef DIM_ERROR

    addclouds_p = instr_addclouds(instrIdx)
    addaerosl_p = instr_addaerosl(instrIdx)

#if defined(RTTOV12) && !defined(RADSHARE)
    transmission%l_opdep = lopdep !TEST
    transmission_k%l_opdep = .false. !TEST
#endif

    ! allocate the direct structures if not done yet
    ! ATTENTION: k - matrix has requires as much profiles per
    !            instrument as channels exist
    stat = realloc_rttov_arrays(n_profiles = nchansprofs ,      &
                                    n_levels   = nlevs + nlevs_top, &
                                    n_channels = nchansprofs,       &
                                    profs      = profiles_k,        &
                                    rads       = radiance_k,        &
                                    transm     = transmission_k)
    stat = realloc_rttov_arrays(n_profiles = nchansprofs,       &
                                    n_levels   = nlevs + nlevs_top, &
                                    n_channels = nchansprofs,       &
                                    rads       = radiance,          &
                                    transm     = transmission)

    if (stat /= NO_ERROR) then
       rttov_k_ifc = int(stat)
       return
    endif

    ! clear variable output part of rttov structures
    !emissiv_k(:,:) = 0.0_wp
    call rttov_clear_rad_var(nchansprofs,radiance)
    call rttov_clear_prof   (profiles_k(1:nchansprofs))
#if defined(RTTOV12)
    call rttov_init_transmission(transmission)
    call rttov_init_transmission(transmission_k)
#else
    call rttov_init_transmission(nchansprofs, transmission)
    call rttov_init_transmission(nchansprofs, transmission_k)
#endif

    addclouds_p = instr_addclouds(instrIdx)
    addaerosl_p = instr_addaerosl(instrIdx)

#if defined(RTTOV12)
    profiles_k(1:nchansprofs)% gas_units = default_gas_units
#endif

#if defined(RTTOV12)
    calcemis(1:nchansprofs) = emissivity(1:nchansprofs)%emis_in < 0.01_jprb
#elif defined(RTTOV10)
    calcemis(1:nchansprofs) = emissivity(1:nchansprofs) < 0.01_jprb
#endif

    rttov_errorstatus(:) = errorstatus_success

    !---------------------------------------------------
    ! TODO: the following may be optimized,
    !       we only have to clear 1:nchansprofs elements
    ! set the requested perturbations
    !---------------------------------------------------
    call rttov_clear_rad_var(nchansprofs,radiance_k)

    if (switchrad_p) then
      radiance_k% bt(1:nchansprofs)    = 1.0_jprb
      radiance_k% total(1:nchansprofs) = 0.0_jprb
    else
      radiance_k% bt(1:nchansprofs)    = 0.0_jprb
      radiance_k% total(1:nchansprofs) = 1.0_jprb
    endif

    if (present(conv_overc)) then
#if defined(RTTOV12)
      conv_overc_p = conv_overc
#endif
    end if

    !---------------------
    ! fill in the chanprof
    !---------------------
    do irun1=1,nchansprofs
      chanprof(irun1)% chan = chans(irun1)
      chanprof(irun1)% prof = lprofs_aux(irun1)
    enddo

#if defined(RTTOV12)
#ifdef RTTOV_USE_OPENMP
    call rttov_parallel_k(                                  &
           rttov_errorstatus(1)                            ,& !  --> error flag
           chanprof(1:nchansprofs)                         ,& ! <-- channels and profiles to calculate
           opts                                            ,& ! <-- options
           pack(profiles(iprof_s:iprof_e), mask=prof_used) ,& ! <--  profile array
           profiles_k(1:nchansprofs)                       ,& ! <--> profile increments;
           coefs(instrIdx)                                 ,& ! <--  coefs array
           transmission                                    ,& ! <--> transmittances
           transmission_k                                  ,& ! <--> K matrix of transmittances
           radiance                                        ,& ! <--> direct model output radiances
           radiance_k                                      ,& ! <--> K matrix of radiances
           calcemis(1:nchansprofs)                         ,& ! <--  flag for internal emissivity calc
           emissivity(1:nchansprofs)                       ,& ! <-- input emissivities per channel
           emissivity_k(1:nchansprofs)                     ,& ! <--> k matrix on input surface emissivity
           nthreads=omp_get_max_threads()                   )
#else
    call rttov_k(                                           &
           rttov_errorstatus(1)                            ,& !  --> error flag
           chanprof(1:nchansprofs)                         ,& ! <-- channels and profiles to calculate
           opts                                            ,& ! <-- options
           pack(profiles(iprof_s:iprof_e), mask=prof_used) ,& ! <--  profile array
           profiles_k(1:nchansprofs)                       ,& ! <--> profile increments;
           coefs(instrIdx)                                 ,& ! <--  coefs array
           transmission                                    ,& ! <--> transmittances
           transmission_k                                  ,& ! <--> K matrix of transmittances
           radiance                                        ,& ! <--> direct model output radiances
           radiance_k                                      ,& ! <--> K matrix of radiances
           calcemis(1:nchansprofs)                         ,& ! <--  flag for internal emissivity calc
           emissivity(1:nchansprofs)                       ,& ! <-- input emissivities per channel
           emissivity_k(1:nchansprofs)                     )  ! <--> k matrix on input surface emissivity
#endif
#elif defined(RTTOV10)
    call rttov_k(                                           &
           rttov_errorstatus(1)                            ,& !  --> error flag
           chanprof(1:nchansprofs)                         ,& ! <-- channels and profiles to calculate
           opts                                            ,& ! <-- options
           pack(profiles(iprof_s:iprof_e), mask=prof_used) ,& ! <--  profile array
           profiles_k(1:nchansprofs)                       ,& ! <--> profile increments;
           coefs(instrIdx)                                 ,& ! <--  coefs array
           calcemis(1:nchansprofs)                         ,& ! <--  flag for internal emissivity calc
           emissivity(1:nchansprofs)                       ,& ! <-- input emissivities per channel
           emissivity_k(1:nchansprofs)                     ,& ! <--> k matrix on input surface emissivity
           emissivity_out(1:nchansprofs)                   ,& ! -->  emissivities used in calculations
           emissivity_out_k(1:nchansprofs)                 ,& ! <--> k matrix on output surface emissivity
           transmission                                    ,& ! <--> transmittances
           transmission_k                                  ,& ! <--> K matrix of transmittances
           radiance                                        ,& ! <--> direct model output radiances
           radiance_k                                      )  ! <--> K matrix of radiances
#endif

    !----------------------------------------------
    ! check for completely awful rttov error status
    !----------------------------------------------
    if (rttov_errorstatus(1) == errorstatus_fatal) then
       print *, 'after rttov_k -> fatal error'
       rttov_k_ifc = ERROR_RTTOV_CALL
       return
    endif

    ! fill result into output arrays
    FTRACE_BEGIN('rttov_k_ifc:store')
    if (.not.present(istore)) then
      sind   = 0
      eind   = 0
      count1 = 0
      do irun1=1,nprof
        nusedchans = count (lprofs == irun1)
        sind       = eind + 1
        eind       = eind + nusedchans
        do irun2=1,count(lprofs==irun1)
          count1 = count1 + 1
          temp_k(:,irun2,irun1)=dble(profiles_k(count1)%t  (1+nlevs_top:))
          humi_k(:,irun2,irun1)=dble(profiles_k(count1)%q  (1+nlevs_top:))
          if (present(clw_k)) clw_k(:,irun2,irun1)=dble(profiles_k(count1)%clw(1+nlevs_top:))
          if (present( o3_k))  o3_k(:,irun2,irun1)=dble(profiles_k(count1)%o3 (1+nlevs_top:))
          if (present(co2_k)) co2_k(:,irun2,irun1)=dble(profiles_k(count1)%co2(1+nlevs_top:))
          if (present(n2o_k)) n2o_k(:,irun2,irun1)=dble(profiles_k(count1)%n2o(1+nlevs_top:))
          if (present( co_k))  co_k(:,irun2,irun1)=dble(profiles_k(count1)%co (1+nlevs_top:))
          if (present(ch4_k)) ch4_k(:,irun2,irun1)=dble(profiles_k(count1)%ch4(1+nlevs_top:))
          if(nlevs_top==1) then
            !--------------------------
            ! rttov9 compatibility mode
            !--------------------------
            temp_k(1,irun2,irun1) = temp_k(1,irun2,irun1) + dble(profiles_k(count1)%t(1))
            humi_k(1,irun2,irun1) = humi_k(1,irun2,irun1) + dble(profiles_k(count1)%q(1))
            if (present(clw_k)) clw_k(1,irun2,irun1) = clw_k(1,irun2,irun1) + dble(profiles_k(count1)%clw(1))
            if (present( o3_k))  o3_k(1,irun2,irun1) =  o3_k(1,irun2,irun1) + dble(profiles_k(count1)%o3 (1))
            if (present(co2_k)) co2_k(1,irun2,irun1) = co2_k(1,irun2,irun1) + dble(profiles_k(count1)%co2(1))
            if (present(n2o_k)) n2o_k(1,irun2,irun1) = n2o_k(1,irun2,irun1) + dble(profiles_k(count1)%n2o(1))
            if (present( co_k))  co_k(1,irun2,irun1) =  co_k(1,irun2,irun1) + dble(profiles_k(count1)%co (1))
            if (present(ch4_k)) ch4_k(1,irun2,irun1) = ch4_k(1,irun2,irun1) + dble(profiles_k(count1)%ch4(1))
          endif
        enddo
        t2m_k    (1:nusedchans,irun1) = dble(profiles_k(sind:eind)% s2m % t)
        q2m_k    (1:nusedchans,irun1) = dble(profiles_k(sind:eind)% s2m % q)
        stemp_k  (1:nusedchans,irun1) = dble(profiles_k(sind:eind)% skin% t)
        t_b      (1:nusedchans,irun1) = dble(radiance%bt(sind:eind)        )
#if defined(RTTOV12)
        emissiv  (1:nusedchans,irun1) = dble(emissivity  (sind:eind)%emis_out)
        emissiv_k(1:nusedchans,irun1) = dble(emissivity_k(sind:eind)%emis_out)
#else
        emissiv  (1:nusedchans,irun1) = dble(emissivity_out(sind:eind))
        emissiv_k(1:nusedchans,irun1) = dble(emissivity_out_k(sind:eind))
#endif
        if (btest(rad_out,OUT_CSB)) then
          if (present(t_b_clear)) t_b_clear(1:nusedchans,irun1)=dble(radiance%bt_clear(sind:eind))
        end if
        if (btest(rad_out,OUT_ASR)) then
          if (present(rad      )) rad      (1:nusedchans,irun1)=dble(radiance%total   (sind:eind))
        end if
        if (btest(rad_out,OUT_CSR)) then
          if (present(radclear )) radclear (1:nusedchans,irun1)=dble(radiance%clear   (sind:eind))
        end if
        if (present(   radtotal)) radtotal (  1:nchans,irun1)=dble(radiance%clear   (  sind:eind) )
        if (present(radovercast)) then
          radovercast (:,1:nusedchans,irun1)=dble(radiance%overcast(:,sind:eind) )
        endif
        if (l_transm) &
             transm(:,1:nchans,irun1) = dble(transmission%tau_levels(1+nlevs_top:,sind:eind))
#if defined(RTTOV12) && !defined(RADSHARE)
        if (lopdep) &
             opdep(:,1:nchans,irun1) = dble(transmission%opdep_ref(1+nlevs_top:,sind:eind))
#endif
        if (present(    psurf_k))     psurf_k (  1:nusedchans,irun1)=dble(profiles_k(sind:eind)%s2m%p    )
        if (present(     u10m_k))      u10m_k (  1:nusedchans,irun1)=dble(profiles_k(sind:eind)%s2m%u    )
        if (present(     v10m_k))      v10m_k (  1:nusedchans,irun1)=dble(profiles_k(sind:eind)%s2m%v    )
        if (present(  o3_surf_k))   o3_surf_k (  1:nusedchans,irun1)=dble(profiles_k(sind:eind)%s2m%o    )
        if (present(    wfetc_k))     wfetc_k (  1:nusedchans,irun1)=dble(profiles_k(sind:eind)%s2m%wfetc)
        if (present(      ctp_k))       ctp_k (  1:nusedchans,irun1)=dble(profiles_k(sind:eind)%ctp      )
        if (present(cfraction_k)) cfraction_k (  1:nusedchans,irun1)=dble(profiles_k(sind:eind)%cfraction)
      enddo
    else
      do i = 1, nchansprofs
        ! Forward calc
        t_b(istore(i,1),istore(i,2)) = dble(radiance% bt(i))
        if (btest(rad_out,OUT_CSB)) then
          if (present(t_b_clear  )) t_b_clear  (istore(i,1),istore(i,2)) = dble(radiance% bt_clear(i))
        end if
        if (btest(rad_out,OUT_ASR)) then
          if (present(rad        )) rad        (istore(i,1),istore(i,2)) = dble(radiance% total   (i))
        end if
        if (btest(rad_out,OUT_CSR)) then
          if (present(radclear   )) radclear   (istore(i,1),istore(i,2)) = dble(radiance% clear   (i))
        end if
        if (present(radtotal   )) radtotal   (istore(i,1),istore(i,2))   = dble(radiance% clear     (i))
        if (present(radovercast)) radovercast(:,istore(i,1),istore(i,2)) = dble(radiance% overcast(:,i))
        if (l_transm) &
             transm(:,istore(i,1),istore(i,2)) = dble(transmission%tau_levels(1+nlevs_top:,i))
#if defined(RTTOV12)  && !defined(RADSHARE)
        if (lopdep) &
             opdep(:,istore(i,1),istore(i,2)) = dble(transmission%opdep_ref(1+nlevs_top:,i))
#endif
        ! K (profile variables)
        temp_k(:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%t (1+nlevs_top:))
        humi_k(:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%q (1+nlevs_top:))
        if (present(clw_k )) clw_k (:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%clw(1+nlevs_top:))
        if (present(o3_k  )) o3_k  (:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%o3 (1+nlevs_top:))
        if (present(co2_k )) co2_k (:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%co2(1+nlevs_top:))
        if (present(n2o_k )) n2o_k (:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%n2o(1+nlevs_top:))
        if (present(co_k  )) co_k  (:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%co (1+nlevs_top:))
        if (present(ch4_k )) ch4_k (:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%ch4(1+nlevs_top:))
        if(nlevs_top==1) then
          !--------------------------
          ! rttov9 compatibility mode
          !--------------------------
          temp_k(1,istore(i,1),istore(i,2)) = temp_k(1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%t(1))
          humi_k(1,istore(i,1),istore(i,2)) = humi_k(1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%q(1))
          if (present(clw_k )) clw_k (1,istore(i,1),istore(i,2)) = clw_k (1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%clw(1))
          if (present(o3_k  )) o3_k  (1,istore(i,1),istore(i,2)) = o3_k  (1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%o3 (1))
          if (present(co2_k )) co2_k (1,istore(i,1),istore(i,2)) = co2_k (1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%co2(1))
          if (present(n2o_k )) n2o_k (1,istore(i,1),istore(i,2)) = n2o_k (1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%n2o(1))
          if (present(co_k  )) co_k  (1,istore(i,1),istore(i,2)) = co_k  (1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%co (1))
          if (present(ch4_k )) ch4_k (1,istore(i,1),istore(i,2)) = ch4_k (1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%ch4(1))
        endif
        ! K (surface variables)
        t2m_k    (istore(i,1),istore(i,2)) = dble(profiles_k(i)% s2m % t    )
        q2m_k    (istore(i,1),istore(i,2)) = dble(profiles_k(i)% s2m % q    )
        stemp_k  (istore(i,1),istore(i,2)) = dble(profiles_k(i)% skin% t    )
#if defined(RTTOV12)
        emissiv  (istore(i,1),istore(i,2)) = dble(emissivity  (i)%emis_out  )
        emissiv_k(istore(i,1),istore(i,2)) = dble(emissivity_k(i)%emis_out  )
#else
        emissiv  (istore(i,1),istore(i,2)) = dble(emissivity_out(i)         )
        emissiv_k(istore(i,1),istore(i,2)) = dble(emissivity_out_k(i)       )
#endif
        if (present(psurf_k    )) psurf_k    (istore(i,1),istore(i,2)) = dble(profiles_k(i)% s2m % p    )
        if (present(u10m_k     )) u10m_k     (istore(i,1),istore(i,2)) = dble(profiles_k(i)% s2m % u    )
        if (present(u10m_k     )) v10m_k     (istore(i,1),istore(i,2)) = dble(profiles_k(i)% s2m % v    )
        if (present(o3_surf_k  )) o3_surf_k  (istore(i,1),istore(i,2)) = dble(profiles_k(i)% s2m % o    )
        if (present(wfetc_k    )) wfetc_k    (istore(i,1),istore(i,2)) = dble(profiles_k(i)% s2m % wfetc)
        if (present(ctp_k      )) ctp_k      (istore(i,1),istore(i,2)) = dble(profiles_k(i)% ctp        )
        if (present(cfraction_k)) cfraction_k(istore(i,1),istore(i,2)) = dble(profiles_k(i)% cfraction  )
      end do
    end if
    FTRACE_END('rttov_k_ifc:store')

    if (chk_reg_lims /= 0 .and. nprof > 0) then
      stat = realloc_rttov_arrays(nprof_aux,profiles(1)%nlevels,0,profs=profiles_aux)
      if (stat /= NO_ERROR) then
        rttov_k_ifc = int(stat)
        return
      endif
      do irun1 = iprof_s, iprof_e
        if (prof_used(irun1)) then
          call check_prof(irun1, coefs(instrIdx), lims_flag, lpio)
          if (chk_reg_lims > 1 .and. present(reg_lim)) then
            if (.not.present(istore)) then
              do i = 1, nprof
                if (lprofs((i-1)*nchans+1) ==irun1) then
                  reg_lim(:,:,i) = lims_flag(:,:)
                  exit
                end if
              end do
            else
              do i = 1, nchansprofs
                if (lprofs(i)==irun1) then
                  reg_lim(:,:,istore(i,2)) = lims_flag(:,:)
                  exit
                end if
              end do
            end if
          end if
        end if
      end do
      stat = dealloc_rttov_arrays(profs=profiles_aux)
      if (stat /= NO_ERROR) then
        rttov_k_ifc = int(stat)
        return
      endif
    end if

#if defined(RTTOV12)
    if (chk_god /= 0 .and. nprof > 0 .and. present(rflag)) then
#ifndef RADSHARE
      do i = 1, nchansprofs
        stat = 0
        if (associated(coefs(instrIdx)%coef%god)) then
          if (lprofs(i) == ipr) then
            print*,'check_god_ifc prof',profiles(ipr)%id,chans(i),coefs(instrIdx)%coef%ff_ori_chn(chans(i))
            print*,'check_god_ifc ntr',chans(i),any(coefs(instrIdx)%coef%god(:,chans(i))%ntr > 0)
          end if
          if (any(coefs(instrIdx)%coef%god(:,chans(i))%ntr > 0)) then
            call check_god_infl(transmission%tau_levels(:,i),         &
                                coefs(instrIdx)%coef%god(:,chans(i)), &
                                stat, msg_)
            if (stat /= 0) then
              if (btest(chk_god, 0)) then
                p => profiles(lprofs(i))
                write(msg,'(2x,A,I5.5," god-smoothing too large in ",A," (",f7.3,"N",f8.3,"E)")') &
                     'profile '//trim(p%id)//' chan ',coefs(instrIdx)%coef%ff_ori_chn(chans(i)),&
                     trim(msg_),p%latitude,p%longitude
                WR
              end if
            end if
          end if
        end if
        if ((btest(chk_god,1) .or. btest(chk_god,2)) .and. all(shape(rflag) > 0)) then
          if (present(istore)) then
            if (lprofs(i) == ipr) then
              print*,'check_god_ifc stat',stat,istore(i,:)
            end if
            rflag(istore(i,1), istore(i,2)) = stat
          else
            ! TODO !!!
          end if
        end if
      end do
#else
      rttov_k_ifc = ERR_INPUT
      return
#endif
    end if
#endif


#ifdef RTTOV10
    if (any(rttov_errorstatus == errorstatus_warning)) then
       rttov_k_ifc = WARN_RTTOV_DIR_ANY
    endif
    if (all(rttov_errorstatus == errorstatus_warning)) then
       rttov_k_ifc = WARN_RTTOV_DIR_ALL
    endif
#endif

FTRACE_END('rttov_k_ifc')

#endif
#ifndef NO_RTTOV
    if (l_dealloc) then
       stat = dealloc_rttov_arrays( profs      = profiles_k, &
                                    rads       = radiance_k, &
                                    transm     = transmission_k)
       stat = dealloc_rttov_arrays( rads       = radiance,   &
                                    transm     = transmission)
    end if
#endif
#endif
  end function rttov_k_ifc

#ifndef RADSHARE
 subroutine retrieve_emissivity(insidx, lprofs, chans, obs, emis, stat, pe)
   !--------------------------------------------
   ! Dynamically retrieves emissivity
   !--------------------------------------------
   integer,    intent(in)  :: insidx    ! rttov instrument index
   integer,    intent(in)  :: lprofs(:) ! list of profile indices
   integer,    intent(in)  :: chans(:)  ! list of channel indices
   real(wp),   intent(in)  :: obs(:)    ! observed brightness temperature
   real(wp),   intent(out) :: emis(:)   ! computed emissivities
   integer,    intent(out) :: stat      ! error status
   integer,    intent(in) ,optional:: pe

#ifdef RTTOV12
   real(jprb), allocatable :: t_b(:,:)
   real(jprb), allocatable :: emis_in(:,:)
   real(jprb), allocatable :: radclear(:,:)
   real(jprb), allocatable :: radupclear(:,:)
   real(jprb), allocatable :: raddnclear(:,:)
   real(jprb), allocatable :: gamma(:,:)
   real(jprb), allocatable :: obsrad(:)
   real(jprb), allocatable :: radup(:)
   real(jprb), allocatable :: rademi(:)
   integer                 :: k, iprint
#endif

   if (RTTOV_IFC_VERSION < 12) then
     stat = ERR_NO_ATLAS
     return
   end if

#ifdef RTTOV12
   !--------------------------------------------
   ! Check inputs consistency
   !--------------------------------------------
   if (coefs(insidx)% coef% id_sensor /= 2) then
     stat = ERR_INVALID_INSTR
     return
   end if
   if (size(lprofs) /= size(chans) .or. size(lprofs) /= size(chans) .or. &
        size(lprofs) /= size(obs) .or. size(lprofs) /= size(emis) ) then
     stat = ERR_DIM
     return
   end if
   if (any(lprofs(:) /= 1 )) then
     stat = ERR_DIM
     return
   end if     

   !--------------------------------------------
   ! Allocate auxiliary arrays
   !--------------------------------------------
   allocate(t_b       (size(lprofs),1))
   allocate(emis_in   (size(lprofs),1))
   emis_in = 0._wp
   allocate(radclear  (size(lprofs),1))
   allocate(radupclear(size(lprofs),1))
   allocate(raddnclear(size(lprofs),1))
   allocate(gamma     (size(lprofs),1))
   allocate(obsrad    (size(lprofs)))
   allocate(radup     (size(lprofs)))
   allocate(rademi    (size(lprofs)))

   !------------------------------------------------
   ! Call rttov_direct_ifc for first guesses needed
   !-----------------------------------------------
   stat = rttov_direct_ifc (            &
        instrIdx     = insidx,          & ! <--  rttov instrument index
        lProfs       = lprofs,          & ! <--  list of profile indices
        chans        = chans,           & ! <--  list of channel indices
        emissiv      = emis_in,         & ! <--> emissivities -
        t_b          = t_b,             & !  --> calculated brightness
        radclear     = radclear,        & !  --> TOA radiance
        radupclear   = radupclear,      & !  --> TOA upweeling radiance(with emissivity term)
        raddnclear   = raddnclear,      & !  --> downelling radiance at surface
        transmtotal  = gamma )!,        & !  --> calculated transmittance
!        pe           = pe,             & ! <--  pe, debug
!        iprint       = iprint )          ! <--  debug

   !------------------------------------------------
   ! Compute emissivity
   !-----------------------------------------------
   do k = 1, size(chans)
     ! Compute radiance associated to observed brightness temperatures
     call planck(coefs(insidx)%coef% planck1(chans(k)),coefs(insidx)%coef% planck2(chans(k)), &
          obs(k),obsrad(k))
     ! Compute radiance associated to skin temperature
     call planck(coefs(insidx)%coef% planck1(chans(k)),coefs(insidx)%coef% planck2(chans(k)), &
          profiles(1)% skin% t, rademi(k))
     ! Compute upwelling radiance without emission term
     radup(k) = radupclear(k,1)-rademi(K)*emis_in(k,1)*gamma(k,1)
     ! Compute emissivity
     emis(k) = (obsrad(k)-radup(k)-raddnclear(k,1)*gamma(k,1))/(gamma(k,1)*(rademi(k)-raddnclear(k,1)))
   end do

   !--------------------------------------------
   ! Clean up
   !--------------------------------------------
   deallocate(t_b       )
   deallocate(emis_in   )
   deallocate(radclear  )
   deallocate(radupclear)
   deallocate(raddnclear)
   deallocate(gamma     )
   deallocate(obsrad    )
   deallocate(radup     )
   deallocate(rademi    )
#endif

 end subroutine retrieve_emissivity

 subroutine init_rttov_mw_atlas(telsem, cnrm, month, path, my_proc_id, n_proc, io_proc_id, mpi_comm_type, stat)
   !--------------------------------------------
   ! Loads MW emissivity atlas
   !--------------------------------------------
   logical,            intent(in) :: telsem         ! load TELSEM
   integer,            intent(in) :: cnrm(3)        ! instruments for CNRM (AMSU-A, ATMS, MHS)
   integer,            intent(in) :: month          ! Month number
   character(len=128), intent(in) :: path           ! path to atlases
   integer,            intent(in) :: my_proc_id     ! ID of local processors
   integer,            intent(in) :: n_proc         ! no. of processors in mpi communication domain
   integer,            intent(in) :: io_proc_id     ! ID of IO processor
   integer,            intent(in) :: mpi_comm_type  ! mpi communicator type
   integer,            intent(out):: stat           ! error status

   logical   :: tel, cnr(3)
   integer   :: k, j, n_atlas
   integer   :: ierr(4)

   if (RTTOV_IFC_VERSION < 12) then
     stat = ERR_NO_ATLAS
     return
   end if

#ifdef RTTOV12
   tel = .false.
   cnr = .false.
   n_atlas = 0
   ierr = 0

   mpi_my_proc_id = my_proc_id
#if defined(_RTTOV_DO_DISTRIBCOEF)
   if (io_proc_id == mpi_my_proc_id) then
#endif
     do k = 1, size(coefs)
       if (coefs(k)% coef% id_sensor /= 2) cycle! then ! only microwave
       ! TELSEM needs only be loaded once and can the be used by all MW instruments
       if (telsem .and. .not. tel) then
         n_atlas = n_atlas + 1
         if (n_atlas > size(mw_atlas)) then
           ! mw_atlas array too small
           stat = ERR_DIM
           return
         end if
         call rttov_setup_emis_atlas(ierr(1), opts, month, 1, mw_atlas(n_atlas), atlas_id = 1, &
              path = path, coefs = coefs(k))
         tel = .true.
         if ( ierr(1) == 0 .and. dace% lpio ) write(stdout,*) 'TELSEM emissivity atlas initialized.'
       end if
       ! CNRM needs to be loaded for each instument that requires it (AMSU-A/-B, ATMS and/or MHS)
       do j = 1, size(cnrm)
         if (cnrm(j) == coefs(k)% coef% id_inst .and. .not. cnr(j)) then
           n_atlas = n_atlas + 1
           if (n_atlas > size(mw_atlas)) then
             !mw_atlas array too small
             stat = ERR_DIM
             return
           end if
           call rttov_setup_emis_atlas(ierr(j+1), opts, month, 1, mw_atlas(n_atlas), atlas_id = 2, &
                path = path, coefs = coefs(k))
           cnr(j) = .true.
           if ( ierr(j+1) == 0 .and. dace% lpio ) write(stdout,'(1x,A,I3,A)') 'CNRM emissivity atlas for instrument ', &
                coefs(k)% coef% id_inst, ' initialized.'
         end if
       end do
       if ( tel .and. all(cnr)) exit
     end do
#if defined(_RTTOV_DO_DISTRIBCOEF)
   else
     do j = 1, size(mw_atlas)
       call rttov_deallocate_emis_atlas(mw_atlas(j))
     end do
   endif
   ! distribute atlases to all PEs
   if (n_proc > 1 ) then
     call p_bcast(n_atlas,io_proc_id,mpi_comm_type)
     if (dace% lpio ) write(stdout,'(1x,A,I2,A)') 'Distribute ',n_atlas,' emissivity atlases.'
     do j = 1, n_atlas
       call p_bcast(mw_atlas(j),io_proc_id,mpi_comm_type)
     end do
   end if
#endif
   if (present(stat)) stat = sum(ierr)
#endif

 end subroutine init_rttov_mw_atlas

 subroutine get_emis_atlas(insidx, atlas_id, lprofs, chan, emis, stat)
   integer,           intent(in)  :: insidx   !instrument index
   integer,           intent(in)  :: atlas_id ! 1: TELSEM, 2: CNRM
   integer,           intent(in)  :: lprofs(:)! list of profile indices
   integer,           intent(in)  :: chan(:)  ! list of channels
   real(wp),          intent(out) :: emis(:)  ! emissivities
   integer,           intent(out) :: stat     ! error status

   integer                           :: irun, nchansprofs
#ifdef RTTOV12
   type(rttov_chanprof), allocatable :: chanprof(:)
   type(rttov_emis_atlas_data)       :: atlas
#endif

   !--------------------------------------------
   ! Check inputs consistency
   !--------------------------------------------
   if (RTTOV_IFC_VERSION < 12) then
     stat = ERR_NO_ATLAS
     return
   end if

#ifdef RTTOV12
   if (coefs(insidx)% coef% id_sensor /= 2) then
     stat = ERR_INVALID_INSTR
     return
   end if
   if (size(lprofs) /= size(chan) .or. size(lprofs) /= size(chan) .or. &
        size(lprofs) /= size(emis) ) then
     stat = ERR_DIM
     return
   end if
   if (any(lprofs(:) /= 1 )) then
     ! Can only deal with one profile at a time
     stat = ERR_DIM
     return
   end if

   nchansprofs   = size(chan)
   allocate(chanprof(nchansprofs))
   do irun=1,nchansprofs
      chanprof(irun)% chan = chan(irun)
      chanprof(irun)% prof = lprofs(irun)
   enddo

   do irun = 1, size(mw_atlas)
      if (.not. mw_atlas(irun)% init) cycle
      if (.not. mw_atlas(irun)% is_mw) cycle
      if (mw_atlas(irun)% atlas_id /= atlas_id) cycle
      if (atlas_id == 2 .and. mw_atlas(irun)% cnrm_mw_atlas% inst_id /= coefs(insidx)% coef% id_inst) cycle
      atlas = mw_atlas(irun)
   end do
   if ( .not. atlas% init ) then
     stat = ERR_ATLAS_INIT
     return
   end if
   call rttov_get_emis(stat, opts, chanprof, profiles, coefs(insidx), atlas, emis)
   deallocate(chanprof)
#endif

 end subroutine get_emis_atlas
#endif /* RADSHARE */

!------------------------------------------------------------------------------

  function rttov_dealloc_profiles(pe)
  integer, intent(in), OPTIONAL :: pe
  integer                       :: rttov_dealloc_profiles
#ifndef NO_RTTOV
    rttov_dealloc_profiles = NO_ERROR
    if (associated(profiles)) &
         rttov_dealloc_profiles = dealloc_rttov_arrays(profs=profiles)
#else
    rttov_dealloc_profiles = NO_ERROR
#endif
  end function rttov_dealloc_profiles

  !======================================================================
#ifndef NO_RTTOV

  subroutine check_prof(iprof, coefs, flg, lpio)
  integer,           intent(in)  :: iprof
  type(rttov_coefs), intent(in)  :: coefs
  integer,           intent(out) :: flg(:,:)
  logical,           intent(in)  :: lpio

  character(len=256)           :: msg = ''
  type(rttov_profile), pointer :: p1, p2
  integer                      :: ilev, ierr
  logical                      :: assoc1, assoc2, apply_reg_lims_aux

    flg = 0
    if (chk_reg_lims /= 0) then
      apply_reg_lims_aux = apply_reg_lims_p
      apply_reg_lims_p   = .true.
      call rttov_copy_prof(profiles_aux(1:1), profiles(iprof:iprof), larray=.true., lscalar=.true.)
#if defined(RTTOV12)
      call rttov_convert_profile_units(opts, coefs, profiles(iprof:iprof), profiles_aux(1:1))
      call rttov_copy_prof(profiles_aux(2:2), profiles_aux(1:1), larray=.true., lscalar=.true.)
      call rttov_apply_reg_limits(opts, profiles_aux(2:2), profiles_aux(1:1), coefs%coef, coefs%coef_pccomp)
      p1 => profiles_aux(2)
      p2 => profiles_aux(1)
#elif defined(RTTOV10)
      call rttov_checkinput(opts, profiles_aux(1:1), coefs%coef, coefs%coef_pccomp, ierr)
      if (ierr /= 0) then ; flg = 1 ; return ; end if
      p1 => profiles    (iprof)
      p2 => profiles_aux(1)
#endif
      apply_reg_lims_p   = apply_reg_lims_aux
      ! T
      do ilev = 1+nlevs_top, p1%nlevels
        flg(ilev-nlevs_top,1) = 0
        if (mask_lims_t(ilev)) then
          if (p1%t(ilev) < p2%t(ilev)) then
            if (btest(chk_reg_lims,0)) then
              write(msg,'(2x,A,I3," : ",f7.3,"<",f7.3," (",f7.3,"N",f8.3,"E)")') &
                   'profile '//trim(p1%id)//' t exceeds lower bound in level ',ilev,p1%t(ilev),p2%t(ilev),p1%latitude,p1%longitude
              WR
            end if
            flg(ilev-nlevs_top,1) = 1
          elseif (p1%t(ilev) > p2%t(ilev)) then
            if (btest(chk_reg_lims,0)) then
              write(msg,'(2x,A,I3," : ",f7.3,">",f7.3," (",f7.3,"N",f8.3,"E)")') &
                   'profile '//trim(p1%id)//' t exceeds upper bound in level ',ilev,p1%t(ilev),p2%t(ilev),p1%latitude,p1%longitude
              WR
            end if
            flg(ilev-nlevs_top,1) = -1
          end if
        end if
      end do
      ! q
      do ilev = 1+nlevs_top, p1%nlevels
        flg(ilev-nlevs_top,2) = 0
        if (mask_lims_q(ilev)) then
          if (p1%q(ilev) < p2%q(ilev)) then
            if (btest(chk_reg_lims,0)) then
              write(msg,'(2x,A,I3," : ",e13.6,"<",e13.6," (",f7.3,"N",f8.3,"E)")') &
                   'profile '//trim(p1%id)//' q exceeds lower bound in level ',ilev,p1%q(ilev),p2%q(ilev),p1%latitude,p1%longitude
              WR
            end if
            flg(ilev-nlevs_top,2) = 1
          elseif (p1%q(ilev) > p2%q(ilev)) then
            if (btest(chk_reg_lims,0)) then
              write(msg,'(2x,A,I3," : ",e13.6,">",e13.6," (",f7.3,"N",f8.3,"E)")') &
                   'profile '//trim(p1%id)//' q exceeds upper bound in level ',ilev,p1%q(ilev),p2%q(ilev),p1%latitude,p1%longitude
              WR
            end if
            flg(ilev-nlevs_top,2) = -1
          end if
        end if
      end do
      ! o3
!       assoc1 = associated(p1%o3)
!       assoc2 = associated(p2%o3)
!       if (.not.(assoc1 .eqv. assoc2)) flg = 1
!       if (assoc1 .and. assoc2) then
!         do ilev = 1, p1%nlevels
!           if (mask_lims_o3(ilev) .and. p1%o3(ilev) /= p2%o3(ilev)) then
!             print*,'profile '//trim(p1%id)//' o3 out of bounds in level ',ilev
!             flg = 1
!           end if
!         end do
!       end if
      ! TODO: other variables ...
    end if

  end subroutine check_prof

  !======================================================================

  subroutine rttov_clear_prof(profiles)
  type(rttov_profile),intent(inout) :: profiles(:) ! initialized profiles
  !-----------------------------------------------------------------------
  ! initialize rttov rttov_profile structures.
  !-----------------------------------------------------------------------
    integer :: i

    ! Initialize scalar components
    do i=1,size(profiles)
       profiles(i)% skin% surftype  = -1       ! no meaning
       profiles(i)% skin% watertype = -1       ! no meaning
       profiles(i)% skin% t         =  0._jprb ! on temperature
       profiles(i)% skin% fastem    =  0._jprb ! Fastem
#ifdef RTTOV12
       profiles(i)% skin% salinity  =  0._jprb ! ?
       profiles(i)% skin% soil_moisture =  0._jprb ! ?
       profiles(i)% skin% snow_fraction =  0._jprb ! ?
       profiles(i)% skin% foam_fraction =  0._jprb ! ?
#endif
       profiles(i)% s2m% t          =  0._jprb ! temperature
       profiles(i)% s2m% q          =  0._jprb ! WV
       profiles(i)% s2m% o          =  0._jprb ! O3
       profiles(i)% s2m% p          =  0._jprb ! pressure
       profiles(i)% s2m% u          =  0._jprb ! wind components
       profiles(i)% s2m% v          =  0._jprb ! wind components
       profiles(i)% s2m% wfetc      =  0._jprb ! wind fetc
       profiles(i)% zenangle        =  0._jprb ! no meaning
       profiles(i)% azangle         =  0._jprb ! no meaning
       profiles(i)% sunzenangle     =  0._jprb ! no meaning
       profiles(i)% sunazangle      =  0._jprb ! no meaning
       profiles(i)% elevation       =  0._jprb ! no meaning
       profiles(i)% latitude        =  0._jprb ! no meaning
       profiles(i)% ctp             =  0._jprb ! cloud top pressure
       profiles(i)% cfraction       =  0._jprb ! cloud fraction
#ifdef RTTOV12
       profiles(i)% id              =  ''
       profiles(i)% date            =  -1
       profiles(i)% time            =  -1
       profiles(i)% gas_units       = default_gas_units
       profiles(i)% mmr_cldaer      = .false.
       profiles(i)% idg             = default_idg
       profiles(i)% ice_scheme      = default_ice_scheme
       profiles(i)% longitude       =  0._jprb ! no meaning
       profiles(i)% Be              =  0._jprb ! no meaning
       profiles(i)% cosbk           =  0._jprb ! no meaning
#endif
    enddo

    ! Initialize profile components
    do i=1,size(profiles)
!CDIR SHORTLOOP
!CDIR ARRAYCOMB
       profiles(i)% p  (:)          =  0._jprb ! no AD on pressure levels
       profiles(i)% t  (:)          =  0._jprb ! temperature
       profiles(i)% q  (:)          =  0._jprb ! water vapour (ppmv)
       if (associated(profiles(i)%o3))  &
         profiles(i)% o3 (:)        =  0._jprb ! O3 ppmv
       if (associated(profiles(i)%co2)) &
         profiles(i)% co2(:)        =  0._jprb ! CO2 ppmv
       if (associated(profiles(i)%n2o)) &
         profiles(i)% n2o(:)        =  0._jprb ! N2O ppmv
       if (associated(profiles(i)%co))  &
         profiles(i)% co (:)         =  0._jprb ! CO ppmv
       if (associated(profiles(i)%ch4)) &
         profiles(i)% ch4(:)         =  0._jprb ! CH4 ppmv
       if (associated(profiles(i)%clw)) &
         profiles(i)% clw(:)         =  0._jprb ! cloud liquid water (kg/kg)
!CDIR END ARRAYCOMB
    end do

    do i=1,size(profiles)
       if (addaerosl_p) then
          profiles(i)% aerosols(:,:) = 0._jprb ! (iaer,   nlevels)
       end if
       if (addclouds_p) then
          profiles(i)% cloud(:,:)   =  0._jprb ! (ncldtyp,nlevels)
#if defined(RTTOV12)
          profiles(i)% cfrac(:)   =  0._jprb ! (icld   ,nlevels)
#else
          profiles(i)% cfrac(:,:)   =  0._jprb ! (icld   ,nlevels)
#endif
       endif
    end do

  end subroutine rttov_clear_prof
  !======================================================================

  subroutine rttov_clear_rad_var(nchannels,rad)
  integer(jpim),intent(in)           :: nchannels
  type(rttov_radiance),intent(inout) :: rad
  !-----------------------------------------------------------------------
  ! subroutine clears the per event variable part
  ! of an rttov rttov_radiance structure.
  !-----------------------------------------------------------------------
    rad% clear     (1:nchannels)   = 0.0_jprb
    rad% cloudy    (1:nchannels)   = 0.0_jprb
    rad% total     (1:nchannels)   = 0.0_jprb
    rad% bt        (1:nchannels)   = 0.0_jprb
    rad% bt_clear  (1:nchannels)   = 0.0_jprb
    rad% overcast  (:,1:nchannels) = 0.0_jprb
#if defined(RTTOV12)
    rad% refl_clear(1:nchannels)   = 0.0_jprb
    rad% refl      (1:nchannels)   = 0.0_jprb
#else
    rad% upclear   (1:nchannels)   = 0.0_jprb
    rad% dnclear   (1:nchannels)   = 0.0_jprb
    rad% reflclear (1:nchannels)   = 0.0_jprb
#endif
  end subroutine rttov_clear_rad_var
  !======================================================================
#if defined(RTTOV12)
  subroutine rttov_clear_rad2_var(nchannels,rad)
  integer(jpim),intent(in)           :: nchannels
  type(rttov_radiance2),intent(inout) :: rad
  !-----------------------------------------------------------------------
  ! subroutine clears the per event variable part
  ! of an rttov rttov_radiance structure.
  !-----------------------------------------------------------------------
    rad% upclear    (1:nchannels)   = 0.0_jprb
    rad% dnclear    (1:nchannels)   = 0.0_jprb
    rad% refldnclear(1:nchannels)   = 0.0_jprb
    rad% up         (:,1:nchannels) = 0.0_jprb
    rad% down       (:,1:nchannels) = 0.0_jprb
    rad% surf       (:,1:nchannels) = 0.0_jprb
  end subroutine rttov_clear_rad2_var
#endif
  !======================================================================
#if defined(RTTOV12)
  function alloc_rttov_arrays(n_profs,n_levs,n_channels,l_init,&
                              profs,rads,rads2,transm,height)
#else
  function alloc_rttov_arrays(n_profs,n_levs,n_channels,l_init,&
                              profs,rads,transm,height)
#endif
  integer(jpim)            ,intent(in)             :: n_profs        !Number of profiles to allocate
  integer(jpim)            ,intent(in)             :: n_levs
  integer(jpim)            ,intent(in)             :: n_channels     !Number of channels &
                                                                 !(e.g channels/prof * profiles)
  logical(jplm)            ,intent(in)             :: l_init
  type (rttov_profile)     ,pointer      ,optional :: profs(:)
  type (rttov_radiance)    ,intent(inout),optional :: rads
#if defined(RTTOV12)
  type (rttov_radiance2)   ,intent(inout),optional :: rads2
#endif
  type (rttov_transmission),intent(inout),optional :: transm
  real(kind=jprb)          ,pointer      ,optional :: height(:)
  integer(jpim)                                    :: alloc_rttov_arrays
  !-------------------------------------------------------------
  ! function allocates the RTTOV profile and radiance structures
  !-------------------------------------------------------------
    integer(jpim):: n_lev_layers ! This variable is used to switch between
                                 ! RTTOV9 levels and RTTOV10 layers
                                 ! Currently, we just cur ot the surface level
                                 ! as in the italian RTTOV10 modification

FTRACE_BEGIN('alloc_rttov_arrays')

    n_lev_layers    = n_levs - 1

    alloc_rttov_arrays  = ini_rttov_ifc_errHandle ()
    if (alloc_rttov_arrays /= NO_ERROR) return

    alloc_status(:) = 0_jpim

    if (present(profs )) call alloc_profs(profs )

    if (present(rads)) then
      ! allocate radiance structure
#if defined(RTTOV10)
      call rttov_alloc_rad(alloc_status(1),&
           n_channels,     &
           rads,           &
           n_lev_layers,   &
           RTTOV_ALLOC,    &
           init=.true.     )
#elif defined(RTTOV12)
      if (present(rads2)) then
        call rttov_alloc_rad(alloc_status(1),&
             n_channels,     &
             rads,           &
             n_levs,         &
             RTTOV_ALLOC,    &
             rads2,          &
             init=.true.     )
      else
        call rttov_alloc_rad(alloc_status(1),&
             n_channels,     &
             rads,           &
             n_levs,         &
             RTTOV_ALLOC,    &
             init=.true.     )
      end if
#endif
      if( any (alloc_status /= 0)) then
        alloc_rttov_arrays = ERR_ALLOC
        return
      endif
    end if

    if (present(transm)) then
      ! allocate parts of transmission structure and error status
      call rttov_alloc_transmission(alloc_status(1) ,&
           transm          ,&
#if defined(RTTOV12)
           n_levs          ,&
#else
           n_levs - 1      ,&
#endif
           n_channels      ,&
           RTTOV_ALLOC     ,&
           .true.)           !Initialize the transmission

      if(any(alloc_status /= 0)) then
        alloc_rttov_arrays = ERR_ALLOC
        return
      end if
    end if

    if (present(height)) then
      allocate(height(n_channels), stat=alloc_status(1))
      if(alloc_status(1) /= 0) then
        alloc_rttov_arrays = ERR_ALLOC
        return
      end if
    end if

FTRACE_END('alloc_rttov_arrays')

  contains

    subroutine alloc_profs(profs)
      type (rttov_profile)     ,pointer :: profs(:)
      integer(jpim) :: loc1 (1), loc2 (1)
      logical       :: l_req        ! Whether allocation of profiles is required

      ! allocate profile structure
      l_req = .not.associated(profs)
      if (.not.l_req) l_req = (size(profs) < n_profs)
      if (.not.l_req) l_req = (size(profs(1)%t) < n_levs)

      if (l_req) then
        if (associated(profs)) then
          call rttov_alloc_prof(alloc_status(1),size(profs),profs,         &
               size(profs(1)%p),opts,RTTOV_DEALLOC,coefs(1))
          deallocate(profs, stat = alloc_status(2))
          if(any (alloc_status(1:2) /= 0)) then
            alloc_rttov_arrays = ERR_ALLOC
            return
          end if
        END if
        allocate(profs(n_profs), stat = alloc_status(1))
        if(alloc_status(1) /= 0) then
          alloc_rttov_arrays = ERR_ALLOC
          return
        end if

        loc1 = maxloc(coefs(:)%coef_scatt_ir%fmv_aer_comp)
        loc2 = maxloc(coefs(:)%coef_scatt_ir%fmv_wcl_comp)

        if (loc1(1) /= loc2(1)) then
          alloc_rttov_arrays = ERR_CLOUD_AERO_MISSM
          return
        endif

        call rttov_alloc_prof (  &
             alloc_status(1) , &
             n_profs         , &
             profs           , &
             n_levs          , &
             opts            , &
             RTTOV_ALLOC     , &
             coefs(loc1(1))  , &
             init      = l_init)

        if(alloc_status(1) /= 0) then
          alloc_rttov_arrays  = ERR_ALLOC
          return
        endif
      else
        call rttov_init_prof(profs(1:n_profs))
      end if

#if defined(RTTOV12)
      profs(:)% gas_units  = default_gas_units
      profs(:)% mmr_cldaer = .false.
#endif
    end subroutine alloc_profs


  end function alloc_rttov_arrays

  !======================================================================
#if defined(RTTOV12)
  function dealloc_rttov_arrays(profs,rads,rads2,transm,height,pe)
#else
  function dealloc_rttov_arrays(profs,rads,transm,height,pe)
#endif
  type(rttov_profile)      ,pointer       ,optional :: profs(:)
  type(rttov_radiance)     ,intent(inout) ,optional :: rads
#if defined(RTTOV12)
  type(rttov_radiance2)    ,intent(inout) ,optional :: rads2
#endif
  type(rttov_transmission) ,intent(inout) ,optional :: transm
  real(kind=jprb)          ,pointer       ,optional :: height(:)
  integer, OPTIONAL :: pe

  integer(jpim)   :: dealloc_rttov_arrays
  !---------------------------------------------------------------
  ! function deallocates the RTTOV profile and radiance structures
  !---------------------------------------------------------------
    integer                    :: errorstatus

FTRACE_BEGIN('dealloc_rttov_arrays')

    dealloc_rttov_arrays = ini_rttov_ifc_errHandle()
    if (dealloc_rttov_arrays /= NO_ERROR) return
    errorstatus = errorstatus_success

    if (present(rads)) then
#if defined(RTTOV9)
      if (present(rads2)) then
        if (associated(rads%clear) .and. associated(rads2%upclear)) then
          call rttov_alloc_rad  (errorstatus,size(rads%clear),             &
               rads,size(rads% overcast, 1),RTTOV_DEALLOC, rads2)
          if(errorstatus /= errorstatus_success) then
            dealloc_rttov_arrays = ERR_ALLOC
            return
          endif
        end if
      else
#endif
        if (associated(rads%clear)) then
          call rttov_alloc_rad  (errorstatus,size(rads%clear),             &
               rads,size(rads% overcast, 1),RTTOV_DEALLOC)
          if(errorstatus /= errorstatus_success) then
            dealloc_rttov_arrays = ERR_ALLOC
            return
          endif
        end if
#if defined(RTTOV9)
      end if
#endif
    endif
    if (present(profs)) then
      if (associated(profs)) then
        call rttov_alloc_prof(errorstatus,size(profs),profs,         &
             size(profs(1)%p),opts,RTTOV_DEALLOC,coefs(1))
        if(errorstatus /= errorstatus_success) then
          dealloc_rttov_arrays = ERR_ALLOC
          return
        endif
        deallocate(profs)
        nullify(profs)
      endif
    end if

    if (present(transm)) then
#if defined(RTTOV10)
      if (associated(transm%tau_levels)) then
        call rttov_alloc_transmission(alloc_status(1),transm, &
                                      size(transm%tau_levels,1) - 1, &
                                      size(transm%tau_levels,2),     &
                                      RTTOV_DEALLOC)
      end if
#elif defined(RTTOV12)
      if (associated(transm%tau_levels)) then
        call rttov_alloc_transmission(alloc_status(1),transm, &
                                      size(transm%tau_levels,1) - 1, &
                                      size(transm%tau_levels,2),     &
                                      RTTOV_DEALLOC)
      end if
#endif

      if (present(height)) then
        deallocate(height, stat=dealloc_rttov_arrays)
      end if

    end if

FTRACE_END('dealloc_rttov_arrays')

  end function dealloc_rttov_arrays

#if defined(RTTOV12)
  function realloc_rttov_arrays(n_profiles, n_levels, n_channels, &
       profs,rads,rads2,transm,height)
#else
  function realloc_rttov_arrays(n_profiles, n_levels, n_channels, &
       profs,rads,transm,height)
#endif
  integer(jpim),      intent(in)           :: n_profiles     ! requested number of profiles
  integer(jpim),      intent(in)           :: n_levels       ! requested levels
  integer(jpim),      intent(in)           :: n_channels     ! requested channels
  type(rttov_profile),      pointer,       optional :: profs(:)       ! profile array to realloc
  type(rttov_radiance),     intent(inout), optional :: rads           ! radiance array to realloc
#if defined(RTTOV12)
  type(rttov_radiance2),    intent(inout), optional :: rads2          ! radiance array to realloc
#endif
  type(rttov_transmission), intent(inout), optional :: transm         ! transmission array to realloc
  real(kind=jprb),          pointer,       optional :: height(:)

  integer(jpim)                         :: realloc_rttov_arrays
  integer                               :: n_p
  !----------------------------------------------------------------------------------
  ! function reallocates the profile, radiance and transmission structure if necessary
  !----------------------------------------------------------------------------------
    logical :: l_req, l_req2

    realloc_rttov_arrays = NO_ERROR

    if (present(profs)) then
      l_req = .not.associated(profs)
      if (.not.l_req) then
        n_p = size(profs)
        l_req = (n_profiles > n_p)
        if (.not.l_req .and. n_p > 0) then
          l_req = n_levels /= profs(1)%nlevels
          ! RTTOV10 does not allocate gas profiles by default
          if (.not.l_req) l_req = ozone_data_p .neqv. associated(profs(1)%o3)
          if (.not.l_req) l_req = co2_data_p   .neqv. associated(profs(1)%co2)
          if (.not.l_req) l_req = n2o_data_p   .neqv. associated(profs(1)%n2o)
          if (.not.l_req) l_req = co_data_p    .neqv. associated(profs(1)%co)
          if (.not.l_req) l_req = ch4_data_p   .neqv. associated(profs(1)%ch4)
          if (.not.l_req) l_req = clw_data_p   .neqv. associated(profs(1)%clw)
          if (.not.l_req) l_req = addaerosl_p  .neqv. associated(profs(1)%aerosols)
          if (.not.l_req) l_req = addclouds_p  .neqv. associated(profs(1)%cloud)
        end if
      end if
      if (l_req) then
        realloc_rttov_arrays = dealloc_rttov_arrays(profs =profs)
        realloc_rttov_arrays = alloc_rttov_arrays (n_profiles,&
                                                   n_levels  ,&
                                                   n_channels,&
                                                   .true.    ,&
                                                   profs =profs)
      end if
    end if

    if (present(transm)) then
      l_req = .not.associated(transm%tau_total)
      if (.not.l_req) l_req = n_channels /= size(transm%tau_total)
      if (l_req) then
        realloc_rttov_arrays = dealloc_rttov_arrays(transm=transm)
        realloc_rttov_arrays = alloc_rttov_arrays (n_profiles,&
                                                   n_levels  ,&
                                                   n_channels,&
                                                   .true.    ,&
                                                   transm =transm)
      end if
    end if

    if (present(rads)) then
#if defined(RTTOV12)
      if (present(rads2)) then
        l_req  = .not.associated(rads%clear)
        l_req2 = .not.associated(rads2%upclear)
        if (.not.l_req)  l_req  = n_channels /= size(rads%clear)
        if (.not.l_req2) l_req2 = n_channels /= size(rads2%upclear)
        if (l_req .and. l_req2) then
          realloc_rttov_arrays = dealloc_rttov_arrays(rads=rads, &
                                                      rads2=rads2)
          realloc_rttov_arrays = alloc_rttov_arrays (n_profiles,&
                                                     n_levels  ,&
                                                     n_channels,&
                                                     .true.    ,&
                                                     rads=rads ,&
                                                     rads2=rads2)
        else if (l_req) then
          realloc_rttov_arrays = dealloc_rttov_arrays(rads=rads)
          realloc_rttov_arrays = alloc_rttov_arrays (n_profiles,&
                                                     n_levels  ,&
                                                     n_channels,&
                                                     .true.    ,&
                                                     rads=rads)
       end if
      else
#endif
        l_req = .not.associated(rads%clear)
        if (.not.l_req) l_req = n_channels /= size(rads%clear)
        if (l_req) then
          realloc_rttov_arrays = dealloc_rttov_arrays(rads=rads)
          realloc_rttov_arrays = alloc_rttov_arrays (n_profiles,&
                                                   n_levels  ,&
                                                   n_channels,&
                                                   .true.    ,&
                                                   rads=rads)
        end if
      end if
#if defined(RTTOV12)
    end if
#endif
    if (present(height)) then
      l_req = .not.associated(height)
      if (.not.l_req) l_req = n_channels /= size(height)
      if (l_req) then
        realloc_rttov_arrays = dealloc_rttov_arrays(height=height)
        realloc_rttov_arrays = alloc_rttov_arrays (n_profiles,&
                                                   n_levels  ,&
                                                   n_channels,&
                                                   .true.    ,&
                                                   height=height)
      end if
    end if

  end function realloc_rttov_arrays

#if defined(RTTOV10)
!! || defined(RTTOV12)
!-------------------------------------------------------
! RTTOV10 Coefficient handling functions (including mpi)
!-------------------------------------------------------
  function rttov_find_coef_file(path,form,instrument,coef_kind,pathCoefs)
  character(len=*),intent(out)    :: path          ! the path to the file
  character(len=*),intent(out)    :: form          ! the path to the file
  integer,intent(in)              :: instrument(:) ! instrument triplet
  character(len=*),intent(in)     :: coef_kind     ! the kind of the coefficients
  character(*),intent(in),optional:: pathCoefs     ! path to coefficient files
  integer                         :: rttov_find_coef_file! function result
  !--------------------------------------------------------
  ! Get the path to the requested coefficient file
  !--------------------------------------------------------
    logical                    :: fileexists

    rttov_find_coef_file = NO_ERROR
    fileexists = .false.

    form = "netcdf"
    call rttov_coeffname(rttov_find_coef_file,   &! -> return code
                         instrument             ,&! <- instrument triplet
#if defined(RTTOV12)
                         trim(coef_kind),        &! <- coef type
                         path                    )! -> filename
#else
                         path                   ,&! -> filename
                         trim(coef_kind),        &! <- coef type
                         form                    )! <- binary file
#endif
    if (rttov_find_coef_file /= 0) return

    if (present(pathCoefs)) then
      if (pathCoefs /= "") &
        path = trim(pathCoefs) // "/" // trim(path)
    endif

! RTTOV12 finds the file independent of the suffix/format
#if ! defined(RTTOV12)
    Inquire (file = path, exist = fileexists)

    if (.not. fileexists) then
      form = "unformatted"
      call rttov_coeffname(rttov_find_coef_file,   &! -> return code
                           instrument             ,&! <- instrument triplet
                           path                   ,&! -> filename
                           trim(coef_kind),        &! <- coef type
                           form                    )! <- binary file

      if (rttov_find_coef_file /= 0) return

      if (present(pathCoefs)) then
        if (pathCoefs /= "") &
          path = trim(pathCoefs) // "/" // trim(path)
      endif

      Inquire (file = path, exist = fileexists)
    endif

    If(.not.fileexists) then
      form = "formatted"
      call rttov_coeffname(rttov_find_coef_file  ,&! -> return code
                           instrument            ,&! <- instrument triplet
                           path                  ,&! -> filename
                           trim(coef_kind),       &! <- coef type
                           form                   )! <- ascii file

      if (rttov_find_coef_file /= 0) return

      if (present(pathCoefs)) then
        if (pathCoefs /= "") &
          path = trim(pathCoefs) // "/" // trim(path)
      endif
    endif
#else
    inquire (file = trim(path)//'.H5', exist = fileexists)
    if (.not.fileexists) then
       inquire (file = trim(path)//'.h5', exist = fileexists)
       if (.not.fileexists) then
          inquire (file = trim(path)//'.dat', exist = fileexists)
          if (.not.fileexists) then
             rttov_find_coef_file = 13
          else
             form = 'formatted'
             path = trim(path)//'.dat'
          end if
       else
          form = 'hdf5'
          path = trim(path)//'.h5'
       end if
    else
       form = 'hdf5'
       path = trim(path)//'.H5'
    end if
#endif
  end function rttov_find_coef_file


  function rttov_read_coef_path(coefs,instrument,channels,pathCoefs)
  type(rttov_coefs),intent(out)   :: coefs            ! the coefficient structure
  integer,intent(in)              :: instrument(:)    ! instrument triplet
  integer,intent(in)              :: channels(:)      ! channels to extract
  character(*),intent(in),optional:: pathCoefs        ! path to coefficient files
  integer                         :: rttov_read_coef_path! function result
  !--------------------------------------------------------
  ! Read in the coefficient files using the given base path
  !--------------------------------------------------------
    character(256)              :: filename_coef   = ''
    character(256)              :: filename_scaer  = ''
    character(256)              :: filename_sccld  = ''
    character(256)              :: filename_pccoef = ''
    character(16)               :: form_coef
    character(16)               :: form_scaer
    character(16)               :: form_sccld
    character(16)               :: form_pccoef

    logical :: l_aer, l_cld, l_pc

    rttov_read_coef_path = rttov_find_coef_file(filename_coef, &
                                                form_coef,     &
                                                instrument,    &
                                                'rtcoef',      &
                                                pathCoefs)

    if (rttov_read_coef_path /= 0) return

    if (addaerosl_p) then
       rttov_read_coef_path = rttov_find_coef_file(filename_scaer,&
                                                   form_scaer,     &
                                                   instrument,    &
                                                   'scaercoef',   &
                                                   pathCoefs)
       if (rttov_read_coef_path /= 0) return
    end if

    if (addclouds_p) then
       rttov_read_coef_path = rttov_find_coef_file(filename_sccld,&
                                                   form_sccld,     &
                                                   instrument,    &
                                                   'sccldcoef',   &
                                                   pathCoefs)
       if (rttov_read_coef_path /= 0) return
    end if

    if (addpc_p) then
       rttov_read_coef_path = rttov_find_coef_file(filename_pccoef,&
                                                   form_pccoef,      &
                                                   instrument,     &
                                                   'pccoef',       &
                                                   pathCoefs)
       if (rttov_read_coef_path /= 0) return
    end if

    if (usedRttov_info_Level >= RTTOV_ERR_INFO) then
      write (*,*) 'Reading RTTOV coefs from ', trim(filename_coef)

      if (addclouds_p) &
        write (*,*) 'Reading RTTOV cloud coefs from ', trim(filename_sccld)

      if (addaerosl_p) &
        write (*,*) 'Reading RTTOV aerosol coefs from ', trim(filename_scaer)
    endif

    call rttov_read_coefs(          &
      rttov_read_coef_path,         &! -> return code
      coefs,                        &! -> coefficients
      opts,                         &! <- options
      channels = channels,          &! <- channels to load
      file_coef   = filename_coef,  &! <- filenames
      file_scaer  = filename_scaer, &
      file_sccld  = filename_sccld, &
      file_pccoef = filename_pccoef,&
      form_coef   = form_coef,      &! <- file formats
      form_scaer  = form_scaer,     &
      form_sccld  = form_sccld,     &
      form_pccoef = form_pccoef     )
  end function
#endif


#if defined(_RTTOV_DO_DISTRIBCOEF)
!--------------------------------
! RTTOV IFC MPI Transfer routines
!--------------------------------
  subroutine p_bcast_rttov_coefs (coefs, source, comm)
  type(rttov_coefs), intent(inout) :: coefs
  integer,           intent(in)    :: source
  integer, optional, intent(in)    :: comm
  !-------------------------------------------------------------------
  ! Broadcast an rttov_coefs structure across all available processors
  !-------------------------------------------------------------------
    call p_bcast(coefs % coef, source, comm)

    if (ANY(instr_addaerosl) .or. ANY(instr_addclouds)) then
      call p_bcast(coefs % coef_scatt_ir,source,comm)
      call p_bcast(coefs % optp,source,comm)
    endif

    if (addpc_p) then
      call p_bcast(coefs % coef_pccomp,source,comm)
    endif
  end subroutine

  subroutine p_bcast_rttov_coef(coef, source, comm)
  type(rttov_coef),  intent(inout) :: coef
  integer,           intent(in)    :: source
  integer, optional, intent(in)    :: comm
  integer :: i
  !------------------------------------------------------------------
  ! Broadcast an rttov_coef structure across all available processors
  !------------------------------------------------------------------
    call p_bcast_rttov_container(coef, source, comm)

    if (mpi_my_proc_id /= source) then
      if (associated(coef%fmv_gas_id       )) allocate(coef%fmv_gas_id       (coef%fmv_gas))
      if (associated(coef%fmv_gas_pos      )) allocate(coef%fmv_gas_pos      (ngases_max))
      if (associated(coef%fmv_var          )) allocate(coef%fmv_var          (coef%fmv_gas))
      if (associated(coef%fmv_lvl          )) allocate(coef%fmv_lvl          (coef%fmv_gas))
      if (associated(coef%fmv_coe          )) allocate(coef%fmv_coe          (coef%fmv_gas))
#if defined(RTTOV10)
      if (associated(coef%gaz_units        )) allocate(coef%gaz_units        (coef%fmv_gas))
#endif

      if (associated(coef%ff_ori_chn       )) allocate(coef%ff_ori_chn       (coef%fmv_chn))
      if (associated(coef%ff_val_chn       )) allocate(coef%ff_val_chn       (coef%fmv_chn))
      if (associated(coef%ff_cwn           )) allocate(coef%ff_cwn           (coef%fmv_chn))
      if (associated(coef%ff_bco           )) allocate(coef%ff_bco           (coef%fmv_chn))
      if (associated(coef%ff_bcs           )) allocate(coef%ff_bcs           (coef%fmv_chn))
      if (associated(coef%ff_gam           )) allocate(coef%ff_gam           (coef%fmv_chn))

#if defined(RTTOV10)
      if (associated(coef%tt_chn           )) allocate(coef%tt_chn           (coef%fmv_chn))
      if (associated(coef%tt_cwn           )) allocate(coef%tt_cwn           (coef%fmv_chn))
#endif
      if (associated(coef%tt_val_chn       )) allocate(coef%tt_val_chn       (coef%fmv_chn))
      if (associated(coef%tt_a0            )) allocate(coef%tt_a0            (coef%fmv_chn))
      if (associated(coef%tt_a1            )) allocate(coef%tt_a1            (coef%fmv_chn))

#if defined(RTTOV10)
      if (associated(coef%ss_chn           )) allocate(coef%ss_chn           (coef%fmv_chn))
      if (associated(coef%ss_cwn           )) allocate(coef%ss_cwn           (coef%fmv_chn))
#elif defined(RTTOV12)
      if (associated(coef%pw_val_chn       )) allocate(coef%pw_val_chn       (coef%fmv_chn))
#endif
      if (associated(coef%ss_val_chn       )) allocate(coef%ss_val_chn       (coef%fmv_chn))
      if (associated(coef%ss_solar_spectrum)) allocate(coef%ss_solar_spectrum(coef%fmv_chn))

#if defined(RTTOV10)
      if (associated(coef%woc_chn          )) allocate(coef%woc_chn          (coef%fmv_chn))
      if (associated(coef%woc_cwn          )) allocate(coef%woc_cwn          (coef%fmv_chn))
#elif defined(RTTOV12)
      if (associated(coef%refl_visnir_ow   )) allocate(coef%refl_visnir_ow   (coef%fmv_chn))
      if (associated(coef%refl_visnir_fw   )) allocate(coef%refl_visnir_fw   (coef%fmv_chn))
#endif
      if (associated(coef%woc_waopc_ow     )) allocate(coef%woc_waopc_ow     (coef%fmv_chn))
      if (associated(coef%woc_waopc_fw     )) allocate(coef%woc_waopc_fw     (coef%fmv_chn))

      if (associated(coef%ws_npoint        )) allocate(coef%ws_npoint        (coef%ws_nomega))
      if (associated(coef%ws_k_omega       )) allocate(coef%ws_k_omega       (coef%ws_nomega))

      if (associated(coef%fastem_polar     )) allocate(coef%fastem_polar     (coef%fmv_chn))

#if defined(RTTOV10)
      if (associated(coef%ssirem_chn       )) allocate(coef%ssirem_chn       (coef%fmv_chn))
#endif
      if (associated(coef%ssirem_a0        )) allocate(coef%ssirem_a0        (coef%fmv_chn))
      if (associated(coef%ssirem_a1        )) allocate(coef%ssirem_a1        (coef%fmv_chn))
      if (associated(coef%ssirem_a2        )) allocate(coef%ssirem_a2        (coef%fmv_chn))
      if (associated(coef%ssirem_xzn1      )) allocate(coef%ssirem_xzn1      (coef%fmv_chn))
      if (associated(coef%ssirem_xzn2      )) allocate(coef%ssirem_xzn2      (coef%fmv_chn))
#if defined(RTTOV12)
      if (associated(coef%iremis_coef      )) allocate(coef%iremis_coef      (coef%iremis_ncoef,coef%fmv_chn))
#endif
    endif

    if (associated(coef%fmv_gas_id       )) call p_bcast(coef%fmv_gas_id       ,source,comm)
    if (associated(coef%fmv_gas_pos      )) call p_bcast(coef%fmv_gas_pos      ,source,comm)
    if (associated(coef%fmv_var          )) call p_bcast(coef%fmv_var          ,source,comm)
    if (associated(coef%fmv_lvl          )) call p_bcast(coef%fmv_lvl          ,source,comm)
    if (associated(coef%fmv_coe          )) call p_bcast(coef%fmv_coe          ,source,comm)
#if defined(RTTOV10)
    if (associated(coef%gaz_units        )) call p_bcast(coef%gaz_units        ,source,comm)
#endif

    if (associated(coef%ff_ori_chn       )) call p_bcast(coef%ff_ori_chn       ,source,comm)
    if (associated(coef%ff_val_chn       )) call p_bcast(coef%ff_val_chn       ,source,comm)
    if (associated(coef%ff_cwn           )) call p_bcast(coef%ff_cwn           ,source,comm)
    if (associated(coef%ff_bco           )) call p_bcast(coef%ff_bco           ,source,comm)
    if (associated(coef%ff_bcs           )) call p_bcast(coef%ff_bcs           ,source,comm)
    if (associated(coef%ff_gam           )) call p_bcast(coef%ff_gam           ,source,comm)

#if defined(RTTOV10)
    if (associated(coef%tt_chn           )) call p_bcast(coef%tt_chn           ,source,comm)
    if (associated(coef%tt_cwn           )) call p_bcast(coef%tt_cwn           ,source,comm)
#endif
    if (associated(coef%tt_val_chn       )) call p_bcast(coef%tt_val_chn       ,source,comm)
    if (associated(coef%tt_a0            )) call p_bcast(coef%tt_a0            ,source,comm)
    if (associated(coef%tt_a1            )) call p_bcast(coef%tt_a1            ,source,comm)

#if defined(RTTOV10)
    if (associated(coef%ss_chn           )) call p_bcast(coef%ss_chn           ,source,comm)
    if (associated(coef%ss_cwn           )) call p_bcast(coef%ss_cwn           ,source,comm)
#elif defined(RTTOV12)
    if (associated(coef%pw_val_chn       )) call p_bcast(coef%pw_val_chn       ,source,comm)
#endif
    if (associated(coef%ss_val_chn       )) call p_bcast(coef%ss_val_chn       ,source,comm)
    if (associated(coef%ss_solar_spectrum)) call p_bcast(coef%ss_solar_spectrum,source,comm)

#if defined(RTTOV10)
    if (associated(coef%woc_chn          )) call p_bcast(coef%woc_chn          ,source,comm)
    if (associated(coef%woc_cwn          )) call p_bcast(coef%woc_cwn          ,source,comm)
#elif defined(RTTOV12)
    if (associated(coef%refl_visnir_ow   )) call p_bcast(coef%refl_visnir_ow   ,source,comm)
    if (associated(coef%refl_visnir_fw   )) call p_bcast(coef%refl_visnir_fw   ,source,comm)
#endif
    if (associated(coef%woc_waopc_ow     )) call p_bcast(coef%woc_waopc_ow     ,source,comm)
    if (associated(coef%woc_waopc_fw     )) call p_bcast(coef%woc_waopc_fw     ,source,comm)

    if (associated(coef%ws_npoint        )) call p_bcast(coef%ws_npoint        ,source,comm)
    if (associated(coef%ws_k_omega       )) call p_bcast(coef%ws_k_omega       ,source,comm)

    if (associated(coef%fastem_polar     )) call p_bcast(coef%fastem_polar     ,source,comm)

#if defined(RTTOV10)
    if (associated(coef%ssirem_chn       )) call p_bcast(coef%ssirem_chn       ,source,comm)
#endif
    if (associated(coef%ssirem_a0        )) call p_bcast(coef%ssirem_a0        ,source,comm)
    if (associated(coef%ssirem_a1        )) call p_bcast(coef%ssirem_a1        ,source,comm)
    if (associated(coef%ssirem_a2        )) call p_bcast(coef%ssirem_a2        ,source,comm)
    if (associated(coef%ssirem_xzn1      )) call p_bcast(coef%ssirem_xzn1      ,source,comm)
    if (associated(coef%ssirem_xzn2      )) call p_bcast(coef%ssirem_xzn2      ,source,comm)
#if defined(RTTOV12)
    if (associated(coef%iremis_coef      )) call p_bcast(coef%iremis_coef      ,source,comm)
#endif

    if (mpi_my_proc_id /= source) then
#if defined(RTTOV12)
      if (associated(coef%bkg_prfl_mr   )) allocate(coef%bkg_prfl_mr   (coef%fmv_lvl(gas_id_mixed), coef%fmv_gas))
      if (associated(coef%env_prfl_tmax )) allocate(coef%env_prfl_tmax (coef%fmv_lvl(gas_id_mixed)              ))
      if (associated(coef%env_prfl_tmin )) allocate(coef%env_prfl_tmin (coef%fmv_lvl(gas_id_mixed)              ))
      if (associated(coef%env_prfl_gmax )) allocate(coef%env_prfl_gmax (coef%fmv_lvl(gas_id_mixed), coef%fmv_gas))
      if (associated(coef%env_prfl_gmin )) allocate(coef%env_prfl_gmin (coef%fmv_lvl(gas_id_mixed), coef%fmv_gas))
#endif

      if (associated(coef%ref_prfl_p    )) allocate(coef%ref_prfl_p    (coef%fmv_lvl(gas_id_mixed)              ))
      if (associated(coef%ref_prfl_t    )) allocate(coef%ref_prfl_t    (coef%fmv_lvl(gas_id_mixed), coef%fmv_gas))
      if (associated(coef%ref_prfl_mr   )) allocate(coef%ref_prfl_mr   (coef%fmv_lvl(gas_id_mixed), coef%fmv_gas))

      if (associated(coef%lim_prfl_p    )) allocate(coef%lim_prfl_p    (coef%fmv_lvl(gas_id_mixed)              ))
      if (associated(coef%lim_prfl_tmax )) allocate(coef%lim_prfl_tmax (coef%fmv_lvl(gas_id_mixed)              ))
      if (associated(coef%lim_prfl_tmin )) allocate(coef%lim_prfl_tmin (coef%fmv_lvl(gas_id_mixed)              ))
      if (associated(coef%lim_prfl_gmax )) allocate(coef%lim_prfl_gmax (coef%fmv_lvl(gas_id_mixed), coef%fmv_gas))
      if (associated(coef%lim_prfl_gmin )) allocate(coef%lim_prfl_gmin (coef%fmv_lvl(gas_id_mixed), coef%fmv_gas))

#if defined(RTTOV10)
      if (associated(coef%mixedgas      )) allocate(coef%mixedgas      (coef%nlayers, coef%fmv_chn, coef%ncmixed ))
      if (associated(coef%watervapour   )) allocate(coef%watervapour   (coef%nlayers, coef%fmv_chn, coef%ncwater ))
      if (associated(coef%ozone         )) allocate(coef%ozone         (coef%nlayers, coef%fmv_chn, coef%ncozone ))
      if (associated(coef%wvcont        )) allocate(coef%wvcont        (coef%nlayers, coef%fmv_chn, coef%ncwvcont))
      if (associated(coef%co2           )) allocate(coef%co2           (coef%nlayers, coef%fmv_chn, coef%ncco2   ))
      if (associated(coef%n2o           )) allocate(coef%n2o           (coef%nlayers, coef%fmv_chn, coef%ncn2o   ))
      if (associated(coef%co            )) allocate(coef%co            (coef%nlayers, coef%fmv_chn, coef%ncco    ))
      if (associated(coef%ch4           )) allocate(coef%ch4           (coef%nlayers, coef%fmv_chn, coef%ncch4   ))
      if (associated(coef%mixedgasint   )) allocate(coef%mixedgasint   (2, coef % nintmixed))
      if (associated(coef%watervapourint)) allocate(coef%watervapourint(2, coef % nintwater))
      if (associated(coef%ozoneint      )) allocate(coef%ozoneint      (2, coef % nintozone))
      if (associated(coef%wvcontint     )) allocate(coef%wvcontint     (2, coef % nintwvcont))
      if (associated(coef%co2int        )) allocate(coef%co2int        (2, coef % nintco2))
      if (associated(coef%n2oint        )) allocate(coef%n2oint        (2, coef % nintn2o))
      if (associated(coef%coint         )) allocate(coef%coint         (2, coef % nintco))
      if (associated(coef%ch4int        )) allocate(coef%ch4int        (2, coef % nintch4))
#elif defined(RTTOV12)
!       l_sol2therm = .false.
!       if (associated(coef%solar,coef%thermal)) then
!          l_sol2therm = .true.
!       elseif (associated(coef%solar)) then
!          allocate(coef%solar(coef%fmv_chn))
!       end if
!       if (associated(coef%thermal)) then
!          allocate(coef%thermal(coef%fmv_chn))
!          if (l_sol2therm) coef%solar => coef%thermal
!       end if
      if (associated(coef%thermal       ))  allocate(coef%thermal      (coef%fmv_chn))
      if (coef%solarcoef) then
         if (associated(coef%solar      )) allocate(coef%solar        (coef%fmv_chn))
      else
         coef%solar => coef%thermal
      endif
      if (associated(coef%nlte_coef     )) allocate(coef%nlte_coef                  )
      if (associated(coef%pmc_pnominal  )) allocate(coef%pmc_pnominal (coef%fmv_chn))
      if (associated(coef%pmc_coef      )) allocate(coef%pmc_coef     (coef%pmc_nlay, coef%fmv_chn, coef%pmc_nvar))
      if (associated(coef%pmc_ppmc      )) allocate(coef%pmc_ppmc     (coef%fmv_chn))
#endif

      if (associated(coef%planck1       )) allocate(coef%planck1      (coef%fmv_chn))
      if (associated(coef%planck2       )) allocate(coef%planck2      (coef%fmv_chn))
      if (associated(coef%frequency_ghz )) allocate(coef%frequency_ghz(coef%fmv_chn))

      if (associated(coef%dp            )) allocate(coef%dp           (coef%nlayers))
      if (associated(coef%dpp           )) allocate(coef%dpp          (0:coef%nlayers))
      if (associated(coef%tstar         )) allocate(coef%tstar        (coef%nlayers))
      if (associated(coef%to3star       )) allocate(coef%to3star      (coef%nlayers))
      if (associated(coef%wstar         )) allocate(coef%wstar        (coef%nlayers))
      if (associated(coef%ostar         )) allocate(coef%ostar        (coef%nlayers))
      if (associated(coef%co2star       )) allocate(coef%co2star      (coef%nlayers))
      if (associated(coef%n2ostar       )) allocate(coef%n2ostar      (coef%nlayers))
      if (associated(coef%costar        )) allocate(coef%costar       (coef%nlayers))
      if (associated(coef%ch4star       )) allocate(coef%ch4star      (coef%nlayers))
#if defined(RTTOV12)
      if (associated(coef%so2star       )) allocate(coef%so2star      (coef%nlayers))
      if (associated(coef%bounds        )) allocate(coef%bounds       (2, coef%fmv_gas, coef%fmv_chn, 2))
#endif
    endif

#if defined(RTTOV12)
    if (associated(coef%bkg_prfl_mr   )) call p_bcast(coef%bkg_prfl_mr   ,source,comm)
    if (associated(coef%env_prfl_tmax )) call p_bcast(coef%env_prfl_tmax ,source,comm)
    if (associated(coef%env_prfl_tmin )) call p_bcast(coef%env_prfl_tmin ,source,comm)
    if (associated(coef%env_prfl_gmax )) call p_bcast(coef%env_prfl_gmax ,source,comm)
    if (associated(coef%env_prfl_gmin )) call p_bcast(coef%env_prfl_gmin ,source,comm)
#endif

    if (associated(coef%ref_prfl_p    )) call p_bcast(coef%ref_prfl_p    ,source,comm)
    if (associated(coef%ref_prfl_t    )) call p_bcast(coef%ref_prfl_t    ,source,comm)
    if (associated(coef%ref_prfl_mr   )) call p_bcast(coef%ref_prfl_mr   ,source,comm)

    if (associated(coef%lim_prfl_p    )) call p_bcast(coef%lim_prfl_p    ,source,comm)
    if (associated(coef%lim_prfl_tmax )) call p_bcast(coef%lim_prfl_tmax ,source,comm)
    if (associated(coef%lim_prfl_tmin )) call p_bcast(coef%lim_prfl_tmin ,source,comm)
    if (associated(coef%lim_prfl_gmax )) call p_bcast(coef%lim_prfl_gmax ,source,comm)
    if (associated(coef%lim_prfl_gmin )) call p_bcast(coef%lim_prfl_gmin ,source,comm)

#if defined(RTTOV10)
    if (associated(coef%mixedgas      )) call p_bcast(coef%mixedgas      ,source,comm)
    if (associated(coef%watervapour   )) call p_bcast(coef%watervapour   ,source,comm)
    if (associated(coef%ozone         )) call p_bcast(coef%ozone         ,source,comm)
    if (associated(coef%wvcont        )) call p_bcast(coef%wvcont        ,source,comm)
    if (associated(coef%co2           )) call p_bcast(coef%co2           ,source,comm)
    if (associated(coef%n2o           )) call p_bcast(coef%n2o           ,source,comm)
    if (associated(coef%co            )) call p_bcast(coef%co            ,source,comm)
    if (associated(coef%ch4           )) call p_bcast(coef%ch4           ,source,comm)

    if (associated(coef%mixedgasint   )) call p_bcast(coef%mixedgasint   ,source,comm)
    if (associated(coef%watervapourint)) call p_bcast(coef%watervapourint,source,comm)
    if (associated(coef%ozoneint      )) call p_bcast(coef%ozoneint      ,source,comm)
    if (associated(coef%wvcontint     )) call p_bcast(coef%wvcontint     ,source,comm)
    if (associated(coef%co2int        )) call p_bcast(coef%co2int        ,source,comm)
    if (associated(coef%n2o           )) call p_bcast(coef%n2o           ,source,comm)
    if (associated(coef%coint         )) call p_bcast(coef%coint         ,source,comm)
    if (associated(coef%ch4int        )) call p_bcast(coef%ch4int        ,source,comm)
#else
    if (associated(coef%thermal       )) then
       do i = 1, size(coef%thermal)
          call p_bcast(coef%thermal(i) ,source,comm)
       end do
    end if
    if (coef%solarcoef) then
       if (associated(coef%solar      )) then
          do i = 1, size(coef%solar)
             call p_bcast(coef%solar(i) ,source,comm)
          end do
       end if
    end if

    if (associated(coef%nlte_coef     )) call p_bcast(coef%nlte_coef     ,source,comm)
    if (associated(coef%pmc_pnominal  )) call p_bcast(coef%pmc_pnominal  ,source,comm)
    if (associated(coef%pmc_coef      )) call p_bcast(coef%pmc_coef      ,source,comm)
    if (associated(coef%pmc_ppmc      )) call p_bcast(coef%pmc_ppmc      ,source,comm)
#endif

    if (associated(coef%planck1       )) call p_bcast(coef%planck1       ,source,comm)
    if (associated(coef%planck2       )) call p_bcast(coef%planck2       ,source,comm)
    if (associated(coef%frequency_ghz )) call p_bcast(coef%frequency_ghz ,source,comm)

    if (associated(coef%dp            )) call p_bcast(coef%dp            ,source,comm)
    if (associated(coef%dpp           )) call p_bcast(coef%dpp           ,source,comm)
    if (associated(coef%tstar         )) call p_bcast(coef%tstar         ,source,comm)
    if (associated(coef%to3star       )) call p_bcast(coef%to3star       ,source,comm)
    if (associated(coef%wstar         )) call p_bcast(coef%wstar         ,source,comm)
    if (associated(coef%ostar         )) call p_bcast(coef%ostar         ,source,comm)
    if (associated(coef%co2star       )) call p_bcast(coef%co2star       ,source,comm)
    if (associated(coef%n2ostar       )) call p_bcast(coef%n2ostar       ,source,comm)
    if (associated(coef%costar        )) call p_bcast(coef%costar        ,source,comm)
    if (associated(coef%ch4star       )) call p_bcast(coef%ch4star       ,source,comm)
#if defined(RTTOV12)
    if (associated(coef%so2star       )) call p_bcast(coef%so2star       ,source,comm)
    if (associated(coef%bounds        )) call p_bcast(coef%bounds        ,source,comm)
#endif
  end subroutine

  subroutine p_bcast_rttov_coef_scatt_ir(coef, source, comm)
  type(rttov_coef_scatt_ir),intent(inout) :: coef
  integer,                  intent(in)    :: source
  integer, optional,        intent(in)    :: comm
  !------------------------------------------------------------------
  ! Broadcast an rttov_coef structure across all available processors
  !------------------------------------------------------------------
    integer :: count
    integer :: dimensions(21)

    call p_bcast_rttov_container (coef, source, comm)

    dimensions(:) = -1

    if (associated(coef%fmv_aer_rh_val)) dimensions(1:1)   = shape(coef%fmv_aer_rh_val)
    if (associated(coef%fmv_wcl_rh_val)) dimensions(2:2)   = shape(coef%fmv_wcl_rh_val)
#if defined(RTTOV10)
    if (associated(coef%channels_solar)) dimensions(3:3)   = shape(coef%channels_solar)
#endif
    if (associated(coef%abs))            dimensions(4:5)   = shape(coef%abs)
    if (associated(coef%sca))            dimensions(6:7)   = shape(coef%sca)
    if (associated(coef%bpr))            dimensions(8:9)   = shape(coef%bpr)
    if (associated(coef%pha))            dimensions(10:12) = shape(coef%pha)
#if defined(RTTOV12)
    if (associated(coef%fmv_icl_deff ))  dimensions(13:13) = shape(coef%fmv_icl_deff)
    if (associated(coef%aer_pha_index))  dimensions(14:14) = shape(coef%aer_pha_index)
    if (associated(coef%wcl_pha_index))  dimensions(15:15) = shape(coef%wcl_pha_index)
    if (associated(coef%icl_pha_index))  dimensions(16:16) = shape(coef%icl_pha_index)
    if (associated(coef%nmom         ))  dimensions(17:18) = shape(coef%nmom)
    if (associated(coef%legcoef      ))  dimensions(19:21) = shape(coef%legcoef)
#endif

    call p_bcast(dimensions,source,comm)

    if (mpi_my_proc_id /= source) then
#if defined(RTTOV12)
      if (associated(coef%fmv_aer_comp_name )) allocate(coef%fmv_aer_comp_name (coef%fmv_aer_comp))
      if (associated(coef%fmv_wcl_comp_name )) allocate(coef%fmv_wcl_comp_name (coef%fmv_wcl_comp))
#endif
      if (associated(coef%fmv_aer_rh        )) allocate(coef%fmv_aer_rh        (coef%fmv_aer_comp))
      if (associated(coef%fmv_wcl_rh        )) allocate(coef%fmv_wcl_rh        (coef%fmv_wcl_comp))
      if (associated(coef%fmv_aer_rh_val    )) allocate(coef%fmv_aer_rh_val    (dimensions(1)))
      if (associated(coef%fmv_wcl_rh_val    )) allocate(coef%fmv_wcl_rh_val    (dimensions(2)))
#if defined(RTTOV10)
      if (associated(coef%fmv_wcl_ph_val    )) allocate(coef%fmv_wcl_ph_val    (coef%fmv_wcl_ph))
      if (associated(coef%fmv_wcl_ph_val_cos)) allocate(coef%fmv_wcl_ph_val_cos(coef%fmv_wcl_ph))
      if (associated(coef%fmv_aer_ph_val    )) allocate(coef%fmv_aer_ph_val    (coef%fmv_aer_ph))
      if (associated(coef%fmv_aer_ph_val_cos)) allocate(coef%fmv_aer_ph_val_cos(coef%fmv_aer_ph))
      if (associated(coef%fmv_icl_ph_val    )) allocate(coef%fmv_icl_ph_val    (coef%fmv_icl_ph))
      if (associated(coef%fmv_icl_ph_val_cos)) allocate(coef%fmv_icl_ph_val_cos(coef%fmv_icl_ph))
      if (associated(coef%fmv_icl_dg        )) allocate(coef%fmv_icl_dg        (coef%fmv_icl_comp, coef%fmv_icl_ishp))
      if (associated(coef%channels_solar    )) allocate(coef%channels_solar    (dimensions(3)))
#elif defined(RTTOV12)
      if (associated(coef%fmv_icl_deff      )) allocate(coef%fmv_icl_deff      (dimensions(13)))
      if (associated(coef%aer_pha_chanlist  )) allocate(coef%aer_pha_chanlist  (coef%fmv_aer_pha_chn))
      if (associated(coef%wcl_pha_chanlist  )) allocate(coef%wcl_pha_chanlist  (coef%fmv_wcl_pha_chn))
      if (associated(coef%icl_pha_chanlist  )) allocate(coef%icl_pha_chanlist  (coef%fmv_icl_pha_chn))
      if (associated(coef%aer_pha_index     )) allocate(coef%aer_pha_index     (dimensions(14)))
      if (associated(coef%wcl_pha_index     )) allocate(coef%wcl_pha_index     (dimensions(15)))
      if (associated(coef%icl_pha_index     )) allocate(coef%icl_pha_index     (dimensions(16)))
      if (associated(coef%aer_phangle       )) allocate(coef%aer_phangle       (coef%aer_nphangle))
      if (associated(coef%wcl_phangle       )) allocate(coef%wcl_phangle       (coef%wcl_nphangle))
      if (associated(coef%icl_phangle       )) allocate(coef%icl_phangle       (coef%icl_nphangle))
      if (associated(coef%nmom              )) allocate(coef%nmom              (dimensions(17),dimensions(18)))
      if (associated(coef%legcoef           )) allocate(coef%legcoef           (dimensions(19),dimensions(20),dimensions(21)))
      if (associated(coef%aer_mmr2nd        )) allocate(coef%aer_mmr2nd        (coef%fmv_aer_comp))
#endif
      if (associated(coef%abs               )) allocate(coef%abs               (dimensions(4),dimensions(5)))
      if (associated(coef%sca               )) allocate(coef%sca               (dimensions(6),dimensions(7)))
      if (associated(coef%bpr               )) allocate(coef%bpr               (dimensions(8),dimensions(9)))
      if (associated(coef%pha               )) allocate(coef%pha               (dimensions(10),dimensions(11),dimensions(12)))
      if (associated(coef%confac            )) allocate(coef%confac            (coef%fmv_wcl_comp))
    endif

#if defined(RTTOV12)
    if (associated(coef%fmv_aer_comp_name )) call p_bcast(coef%fmv_aer_comp_name ,source,comm)
    if (associated(coef%fmv_wcl_comp_name )) call p_bcast(coef%fmv_wcl_comp_name ,source,comm)
#endif
    if (associated(coef%fmv_aer_rh        )) call p_bcast(coef%fmv_aer_rh        ,source,comm)
    if (associated(coef%fmv_wcl_rh        )) call p_bcast(coef%fmv_wcl_rh        ,source,comm)
    if (associated(coef%fmv_aer_rh_val    )) call p_bcast(coef%fmv_aer_rh_val    ,source,comm)
    if (associated(coef%fmv_wcl_rh_val    )) call p_bcast(coef%fmv_wcl_rh_val    ,source,comm)
#if defined(RTTOV10)
    if (associated(coef%fmv_wcl_ph_val    )) call p_bcast(coef%fmv_wcl_ph_val    ,source,comm)
    if (associated(coef%fmv_wcl_ph_val_cos)) call p_bcast(coef%fmv_wcl_ph_val_cos,source,comm)
    if (associated(coef%fmv_aer_ph_val    )) call p_bcast(coef%fmv_aer_ph_val    ,source,comm)
    if (associated(coef%fmv_aer_ph_val_cos)) call p_bcast(coef%fmv_aer_ph_val_cos,source,comm)
    if (associated(coef%fmv_icl_ph_val    )) call p_bcast(coef%fmv_icl_ph_val    ,source,comm)
    if (associated(coef%fmv_icl_ph_val_cos)) call p_bcast(coef%fmv_icl_ph_val_cos,source,comm)
    if (associated(coef%fmv_icl_dg        )) call p_bcast(coef%fmv_icl_dg        ,source,comm)
    if (associated(coef%channels_solar    )) call p_bcast(coef%channels_solar    ,source,comm)
#elif defined(RTTOV12)
    if (associated(coef%fmv_icl_deff      )) call p_bcast(coef%fmv_icl_deff      ,source,comm)
    if (associated(coef%aer_pha_chanlist  )) call p_bcast(coef%aer_pha_chanlist  ,source,comm)
    if (associated(coef%wcl_pha_chanlist  )) call p_bcast(coef%wcl_pha_chanlist  ,source,comm)
    if (associated(coef%icl_pha_chanlist  )) call p_bcast(coef%icl_pha_chanlist  ,source,comm)
    if (associated(coef%aer_pha_index     )) call p_bcast(coef%aer_pha_index     ,source,comm)
    if (associated(coef%wcl_pha_index     )) call p_bcast(coef%wcl_pha_index     ,source,comm)
    if (associated(coef%icl_pha_index     )) call p_bcast(coef%icl_pha_index     ,source,comm)
    if (associated(coef%aer_phangle       )) call p_bcast(coef%aer_phangle       ,source,comm)
    if (associated(coef%wcl_phangle       )) call p_bcast(coef%wcl_phangle       ,source,comm)
    if (associated(coef%icl_phangle       )) call p_bcast(coef%icl_phangle       ,source,comm)
    if (associated(coef%nmom              )) call p_bcast(coef%nmom              ,source,comm)
    if (associated(coef%legcoef           )) call p_bcast(coef%legcoef           ,source,comm)
    if (associated(coef%aer_mmr2nd        )) call p_bcast(coef%aer_mmr2nd        ,source,comm)
#endif
    if (associated(coef%abs               )) call p_bcast(coef%abs               ,source,comm)
    if (associated(coef%sca               )) call p_bcast(coef%sca               ,source,comm)
    if (associated(coef%bpr               )) call p_bcast(coef%bpr               ,source,comm)
    if (associated(coef%pha               )) call p_bcast(coef%pha               ,source,comm)
    if (associated(coef%confac            )) call p_bcast(coef%confac            ,source,comm)

#if defined(RTTOV10)
    if (mpi_my_proc_id /= source) then
      if (associated(coef%ifmv_icl_ph_val)) then
        count = int (coef%fmv_icl_ph_val(coef%fmv_icl_ph) / &
                     coef%fmv_icl_ph_val_min, jpim)
        allocate(coef%ifmv_icl_ph_val(count))
      endif

      if (associated(coef%ifmv_wcl_ph_val)) then
        count = int(coef%fmv_wcl_ph_val(coef%fmv_wcl_ph) / &
                    coef%fmv_wcl_ph_val_min, jpim)
        allocate(coef%ifmv_wcl_ph_val(count))
      endif

      if (associated(coef%ifmv_aer_ph_val)) then
        count = int(coef%fmv_aer_ph_val(coef%fmv_aer_ph) / &
                    coef%fmv_aer_ph_val_min, jpim)
        allocate(coef%ifmv_aer_ph_val(count))
      endif
    endif
    if (associated(coef%ifmv_icl_ph_val)) call p_bcast(coef%ifmv_icl_ph_val,source,comm)
    if (associated(coef%ifmv_wcl_ph_val)) call p_bcast(coef%ifmv_wcl_ph_val,source,comm)
    if (associated(coef%ifmv_aer_ph_val)) call p_bcast(coef%ifmv_aer_ph_val,source,comm)
#elif defined(RTTOV12)
    call p_bcast(coef%aer_phfn_int,source,comm)
    call p_bcast(coef%wcl_phfn_int,source,comm)
    call p_bcast(coef%icl_phfn_int,source,comm)
#endif

  end subroutine

  subroutine p_bcast_rttov_optpar_ir(optp, source, comm)
  type(rttov_optpar_ir),intent(inout) :: optp
  integer,              intent(in)    :: source
  integer, optional,    intent(in)    :: comm
  !-----------------------------------------------------------------------
  ! Broadcast an rttov_optpar_ir structure across all available processors
  !-----------------------------------------------------------------------
    integer :: run1
    integer :: dimensions(3)


    call p_bcast_rttov_container (optp, source, comm)

    dimensions(:) = -1

    if(associated(optp % optpaer))  &
      dimensions(1:1) = shape(optp % optpaer)
    if(associated(optp % optpwcl))  &
      dimensions(2:2) = shape(optp % optpwcl)
#if defined(RTTOV10)
    if(associated(optp % optpicl))  &
      dimensions(3:3) = shape(optp % optpicl)
#endif

    call p_bcast(dimensions,source,comm)

    if (mpi_my_proc_id /= source) then
      if(associated(optp % optpaer)) allocate(optp % optpaer(dimensions(1)))
      if(associated(optp % optpwcl)) allocate(optp % optpwcl(dimensions(2)))
#if defined(RTTOV10)
      if(associated(optp % optpicl)) allocate(optp % optpicl(dimensions(3)))
#elif defined(RTTOV12)
      if(associated(optp % optpicl )) allocate(optp % optpicl )
      if(associated(optp % optpiclb)) allocate(optp % optpiclb)
#endif
    endif

    if(associated(optp % optpaer)) then
      do run1=1,dimensions(1)
        call p_bcast(optp % optpaer(run1),source,comm)
      enddo
    endif

    if(associated(optp % optpwcl)) then
      do run1=1,dimensions(2)
        call p_bcast(optp % optpwcl(run1),source,comm)
      enddo
    endif

    if(associated(optp % optpicl)) then
#if defined(RTTOV10)
      do run1=1,dimensions(3)
        call p_bcast(optp % optpicl(run1),source,comm)
      enddo
#elif defined(RTTOV12)
      call p_bcast(optp % optpicl ,source,comm)
      call p_bcast(optp % optpiclb,source,comm)
#endif
    endif
  end subroutine

#if defined(RTTOV12)
  subroutine p_bcast_rttov_coef_optpiclb(coef, source, comm)
  type(rttov_coef_optpiclb),intent(inout) :: coef
  integer,                  intent(in)    :: source
  integer,     optional,    intent(in)    :: comm

    integer :: dimensions(1)

    call p_bcast_rttov_container (coef, source, comm)

    dimensions = -1

    if (associated(coef%iwn2014)) dimensions(1:1) = shape(coef%iwn2014)

    call p_bcast(dimensions,source,comm)

    if (mpi_my_proc_id /= source) then
      if (associated(coef%iwn2014   )) allocate(coef%iwn2014   (dimensions(1)))
      if (associated(coef%jwn2014   )) allocate(coef%jwn2014   (dimensions(1)))
      if (associated(coef%dx_dwn2014)) allocate(coef%dx_dwn2014(dimensions(1)))
      if (associated(coef%q         )) allocate(coef%q         (baran_ngauss))
      if (associated(coef%w         )) allocate(coef%w         (baran_ngauss))
   end if

    if (associated(coef%iwn2014   )) call p_bcast(coef%iwn2014   ,source,comm)
    if (associated(coef%jwn2014   )) call p_bcast(coef%jwn2014   ,source,comm)
    if (associated(coef%dx_dwn2014)) call p_bcast(coef%dx_dwn2014,source,comm)
    if (associated(coef%q         )) call p_bcast(coef%q         ,source,comm)
    if (associated(coef%w         )) call p_bcast(coef%w         ,source,comm)

   call p_bcast(coef%baran_phfn_int, source, comm)

  end subroutine p_bcast_rttov_coef_optpiclb

  subroutine p_bcast_rttov_phasefn_int(phfn, source, comm)
  type(rttov_phasefn_int),intent(inout) :: phfn
  integer,                intent(in)    :: source
  integer,     optional,  intent(in)    :: comm

    integer :: dimensions(2)

    call p_bcast_rttov_container (phfn, source, comm)

    dimensions = -1

    if (associated(phfn%cosphangle)) dimensions(1:1) = shape(phfn%cosphangle)
    if (associated(phfn%iphangle  )) dimensions(2:2) = shape(phfn%iphangle)

    call p_bcast(dimensions,source,comm)

    if (mpi_my_proc_id /= source) then
      if (associated(phfn%cosphangle)) allocate(phfn%cosphangle(dimensions(1)))
      if (associated(phfn%iphangle  )) allocate(phfn%iphangle  (dimensions(2)))
    end if

    if (associated(phfn%cosphangle)) call p_bcast(phfn%cosphangle,source,comm)
    if (associated(phfn%iphangle  )) call p_bcast(phfn%iphangle  ,source,comm)

  end subroutine p_bcast_rttov_phasefn_int

  subroutine p_bcast_rttov_fast_coef(coef, source, comm)
  type(rttov_fast_coef),intent(inout) :: coef
  integer,              intent(in)    :: source
  integer,   optional,  intent(in)    :: comm

    integer :: dimensions(19)
    integer :: i

    call p_bcast_rttov_container (coef, source, comm)

    dimensions = -1

    if (associated(coef%mixedgas   )) dimensions( 1: 2) = shape (coef%mixedgas)
    if (associated(coef%watervapour)) dimensions( 3: 4) = shape (coef%watervapour)
    if (associated(coef%ozone      )) dimensions( 5: 6) = shape (coef%ozone)
    if (associated(coef%wvcont     )) dimensions( 7: 8) = shape (coef%wvcont)
    if (associated(coef%co2        )) dimensions( 9:10) = shape (coef%co2)
    if (associated(coef%n2o        )) dimensions(11:12) = shape (coef%n2o)
    if (associated(coef%co         )) dimensions(13:14) = shape (coef%co)
    if (associated(coef%ch4        )) dimensions(15:16) = shape (coef%ch4)
    if (associated(coef%co2        )) dimensions(17:18) = shape (coef%co2)
    if (associated(coef%gasarray   )) dimensions(19:19) = shape (coef%gasarray)

    call p_bcast(dimensions,source,comm)

    if (mpi_my_proc_id /= source) then
      if (associated(coef%mixedgas   )) allocate(coef%mixedgas   (dimensions(1),dimensions(2)))
      if (associated(coef%watervapour)) allocate(coef%watervapour(dimensions(3),dimensions(4)))
      if (associated(coef%ozone      )) allocate(coef%ozone      (dimensions(5),dimensions(6)))
      if (associated(coef%wvcont     )) allocate(coef%wvcont     (dimensions(7),dimensions(8)))
      if (associated(coef%co2        )) allocate(coef%co2        (dimensions(9),dimensions(10)))
      if (associated(coef%n2o        )) allocate(coef%n2o        (dimensions(11),dimensions(12)))
      if (associated(coef%co         )) allocate(coef%co         (dimensions(13),dimensions(14)))
      if (associated(coef%ch4        )) allocate(coef%ch4        (dimensions(15),dimensions(16)))
      if (associated(coef%co2        )) allocate(coef%co2        (dimensions(17),dimensions(18)))
      if (associated(coef%gasarray   )) allocate(coef%gasarray   (dimensions(19)))
   end if

    if (associated(coef%mixedgas   )) call p_bcast(coef%mixedgas   ,source,comm)
    if (associated(coef%watervapour)) call p_bcast(coef%watervapour,source,comm)
    if (associated(coef%ozone      )) call p_bcast(coef%ozone      ,source,comm)
    if (associated(coef%wvcont     )) call p_bcast(coef%wvcont     ,source,comm)
    if (associated(coef%co2        )) call p_bcast(coef%co2        ,source,comm)
    if (associated(coef%n2o        )) call p_bcast(coef%n2o        ,source,comm)
    if (associated(coef%co         )) call p_bcast(coef%co         ,source,comm)
    if (associated(coef%ch4        )) call p_bcast(coef%ch4        ,source,comm)
    if (associated(coef%co2        )) call p_bcast(coef%co2        ,source,comm)
    if (associated(coef%gasarray   )) then
       do i = 1, dimensions(19)
          call p_bcast(coef%gasarray(i),source,comm)
       end do
    end if

  end subroutine p_bcast_rttov_fast_coef

  subroutine p_bcast_rttov_fast_coef_gas(coef, source, comm)
  type(rttov_fast_coef_gas),intent(inout) :: coef
  integer,                  intent(in)    :: source
  integer,       optional,  intent(in)    :: comm

    integer :: dimensions(2)

    call p_bcast_rttov_container (coef, source, comm)

    dimensions = -1

    if (associated(coef%coef)) dimensions(1:2) = shape(coef%coef)

    call p_bcast(dimensions,source,comm)

    if (mpi_my_proc_id /= source) then
      if (associated(coef%coef)) allocate(coef%coef(dimensions(1),dimensions(2)))
   end if

    if (associated(coef%coef)) call p_bcast(coef%coef,source,comm)

  end subroutine p_bcast_rttov_fast_coef_gas

  subroutine p_bcast_rttov_nlte_coef(coef, source, comm)
  type(rttov_nlte_coef),intent(inout) :: coef
  integer,              intent(in)    :: source
  integer,   optional,  intent(in)    :: comm

    call p_bcast_rttov_container (coef, source, comm)

    if (mpi_my_proc_id /= source) then
      if (associated(coef%coef         )) allocate(coef%coef         (coef%ncoef,coef%nsat,coef%nsol,coef%nchan))
      if (associated(coef%sol_zen_angle)) allocate(coef%sol_zen_angle(coef%nsol))
      if (associated(coef%sat_zen_angle)) allocate(coef%sat_zen_angle(coef%nsat))
      if (associated(coef%cos_sol      )) allocate(coef%cos_sol      (coef%nsol))
      if (associated(coef%sec_sat      )) allocate(coef%sec_sat      (coef%nsat))
   end if

   if (associated(coef%coef         )) call p_bcast(coef%coef         ,source,comm)
   if (associated(coef%sol_zen_angle)) call p_bcast(coef%sol_zen_angle,source,comm)
   if (associated(coef%sat_zen_angle)) call p_bcast(coef%sat_zen_angle,source,comm)
   if (associated(coef%cos_sol      )) call p_bcast(coef%cos_sol      ,source,comm)
   if (associated(coef%sec_sat      )) call p_bcast(coef%sec_sat      ,source,comm)

  end subroutine p_bcast_rttov_nlte_coef


#endif


  subroutine p_bcast_rttov_coef_pccomp(coef_pccomp,source,comm)
  type(rttov_coef_pccomp),intent(inout) :: coef_pccomp
  integer,              intent(in)      :: source
  integer, optional,    intent(in)      :: comm
  !-------------------------------------------------------------------------
  ! Broadcast an rttov_coef_pccomp structure across all available processors
  !-------------------------------------------------------------------------
    integer :: run1, run2
    integer :: dimensions(4)

    call p_bcast_rttov_container (coef_pccomp, source, comm)

    dimensions = -1

    if (associated(coef_pccomp%co2_pc_ref)) dimensions(1:1) = shape(coef_pccomp%co2_pc_ref)
    if (associated(coef_pccomp%n2o_pc_ref)) dimensions(2:2) = shape(coef_pccomp%n2o_pc_ref)
    if (associated(coef_pccomp%co_pc_ref))  dimensions(3:3) = shape(coef_pccomp%co_pc_ref)
    if (associated(coef_pccomp%ch4_pc_ref)) dimensions(4:4) = shape(coef_pccomp%ch4_pc_ref)

    call p_bcast(dimensions,source,comm)

    if (mpi_my_proc_id /= source) then
#if defined(RTTOV10)
      if (associated(coef_pccomp%eigenvectors     )) allocate(coef_pccomp%eigenvectors     (coef_pccomp%fmv_pc_nchn,&
                                                                                            coef_pccomp%fmv_pc_mnum))
#elif defined(RTTOV12)
      if (associated(coef_pccomp%fmv_pc_sets      )) allocate(coef_pccomp%emiss_chn        (coef_pccomp%fmv_pc_bands))
#endif
      if (associated(coef_pccomp%emiss_chn        )) allocate(coef_pccomp%emiss_chn        (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c1         )) allocate(coef_pccomp%emiss_c1         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c2         )) allocate(coef_pccomp%emiss_c2         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c3         )) allocate(coef_pccomp%emiss_c3         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c4         )) allocate(coef_pccomp%emiss_c4         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c5         )) allocate(coef_pccomp%emiss_c5         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c6         )) allocate(coef_pccomp%emiss_c6         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c7         )) allocate(coef_pccomp%emiss_c7         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c8         )) allocate(coef_pccomp%emiss_c8         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c9         )) allocate(coef_pccomp%emiss_c9         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%ref_pc_prfl_p    )) allocate(coef_pccomp%ref_pc_prfl_p    (coef_pccomp%fmv_pc_nlev))
      if (associated(coef_pccomp%ref_pc_prfl_mr   )) allocate(coef_pccomp%ref_pc_prfl_mr   (coef_pccomp%fmv_pc_nlev,&
                                                                                            coef_pccomp%fmv_pc_gas))
      if (associated(coef_pccomp%lim_pc_prfl_tmin )) allocate(coef_pccomp%lim_pc_prfl_tmin (coef_pccomp%fmv_pc_nlev))
      if (associated(coef_pccomp%lim_pc_prfl_tmax )) allocate(coef_pccomp%lim_pc_prfl_tmax (coef_pccomp%fmv_pc_nlev))
      if (associated(coef_pccomp%lim_pc_prfl_qmin )) allocate(coef_pccomp%lim_pc_prfl_qmin (coef_pccomp%fmv_pc_nlev))
      if (associated(coef_pccomp%lim_pc_prfl_qmax )) allocate(coef_pccomp%lim_pc_prfl_qmax (coef_pccomp%fmv_pc_nlev))
      if (associated(coef_pccomp%lim_pc_prfl_ozmin)) allocate(coef_pccomp%lim_pc_prfl_ozmin(coef_pccomp%fmv_pc_nlev))
      if (associated(coef_pccomp%lim_pc_prfl_ozmax)) allocate(coef_pccomp%lim_pc_prfl_ozmax(coef_pccomp%fmv_pc_nlev))
      if (associated(coef_pccomp%co2_pc_ref       )) allocate(coef_pccomp%co2_pc_ref       (dimensions(1)))
      if (associated(coef_pccomp%n2o_pc_ref       )) allocate(coef_pccomp%n2o_pc_ref       (dimensions(2)))
      if (associated(coef_pccomp%co_pc_ref        )) allocate(coef_pccomp%co_pc_ref        (dimensions(3)))
      if (associated(coef_pccomp%ch4_pc_ref       )) allocate(coef_pccomp%ch4_pc_ref       (dimensions(4)))
      if (associated(coef_pccomp%noise            )) allocate(coef_pccomp%noise            (coef_pccomp%fmv_pc_nchn_noise))
      if (associated(coef_pccomp%noise_in         )) allocate(coef_pccomp%noise_in         (coef_pccomp%fmv_pc_nchn))
      if (associated(coef_pccomp%ff_ori_chn_in    )) allocate(coef_pccomp%ff_ori_chn_in    (coef_pccomp%fmv_pc_nchn))
      if (associated(coef_pccomp%ff_cwn_in        )) allocate(coef_pccomp%ff_cwn_in        (coef_pccomp%fmv_pc_nchn))
      if (associated(coef_pccomp%ff_bco_in        )) allocate(coef_pccomp%ff_bco_in        (coef_pccomp%fmv_pc_nchn))
      if (associated(coef_pccomp%ff_bcs_in        )) allocate(coef_pccomp%ff_bcs_in        (coef_pccomp%fmv_pc_nchn))
      if (associated(coef_pccomp%planck1_in       )) allocate(coef_pccomp%planck1_in       (coef_pccomp%fmv_pc_nchn))
      if (associated(coef_pccomp%planck2_in       )) allocate(coef_pccomp%planck2_in       (coef_pccomp%fmv_pc_nchn))
#if defined(RTTOV12)
      if (associated(coef_pccomp%pcreg            )) allocate(coef_pccomp%pcreg            (coef_pccomp%fmv_pc_bands,&
                                                                                            coef_pccomp%fmv_pc_msets))
      if (associated(coef_pccomp%eigen            )) allocate(coef_pccomp%eigen            (coef_pccomp%fmv_pc_bands))
#else
      if (associated(coef_pccomp%pcreg            )) allocate(coef_pccomp%pcreg            (coef_pccomp%fmv_pc_sets))
#endif
    endif

#if defined(RTTOV10)
    if (associated(coef_pccomp%eigenvectors     )) call p_bcast(coef_pccomp%eigenvectors     ,source,comm)
#elif defined(RTTOV12)
    if (associated(coef_pccomp%fmv_pc_sets      )) call p_bcast(coef_pccomp%fmv_pc_sets      ,source,comm)
#endif
    if (associated(coef_pccomp%emiss_chn        )) call p_bcast(coef_pccomp%emiss_chn        ,source,comm)
    if (associated(coef_pccomp%emiss_c1         )) call p_bcast(coef_pccomp%emiss_c1         ,source,comm)
    if (associated(coef_pccomp%emiss_c2         )) call p_bcast(coef_pccomp%emiss_c2         ,source,comm)
    if (associated(coef_pccomp%emiss_c3         )) call p_bcast(coef_pccomp%emiss_c3         ,source,comm)
    if (associated(coef_pccomp%emiss_c4         )) call p_bcast(coef_pccomp%emiss_c4         ,source,comm)
    if (associated(coef_pccomp%emiss_c5         )) call p_bcast(coef_pccomp%emiss_c5         ,source,comm)
    if (associated(coef_pccomp%emiss_c6         )) call p_bcast(coef_pccomp%emiss_c6         ,source,comm)
    if (associated(coef_pccomp%emiss_c7         )) call p_bcast(coef_pccomp%emiss_c7         ,source,comm)
    if (associated(coef_pccomp%emiss_c8         )) call p_bcast(coef_pccomp%emiss_c8         ,source,comm)
    if (associated(coef_pccomp%emiss_c9         )) call p_bcast(coef_pccomp%emiss_c9         ,source,comm)
    if (associated(coef_pccomp%ref_pc_prfl_p    )) call p_bcast(coef_pccomp%ref_pc_prfl_p    ,source,comm)
    if (associated(coef_pccomp%ref_pc_prfl_mr   )) call p_bcast(coef_pccomp%ref_pc_prfl_mr   ,source,comm)
    if (associated(coef_pccomp%lim_pc_prfl_tmin )) call p_bcast(coef_pccomp%lim_pc_prfl_tmin ,source,comm)
    if (associated(coef_pccomp%lim_pc_prfl_tmax )) call p_bcast(coef_pccomp%lim_pc_prfl_tmax ,source,comm)
    if (associated(coef_pccomp%lim_pc_prfl_qmin )) call p_bcast(coef_pccomp%lim_pc_prfl_qmin ,source,comm)
    if (associated(coef_pccomp%lim_pc_prfl_qmax )) call p_bcast(coef_pccomp%lim_pc_prfl_qmax ,source,comm)
    if (associated(coef_pccomp%lim_pc_prfl_ozmin)) call p_bcast(coef_pccomp%lim_pc_prfl_ozmin,source,comm)
    if (associated(coef_pccomp%lim_pc_prfl_ozmax)) call p_bcast(coef_pccomp%lim_pc_prfl_ozmax,source,comm)
    if (associated(coef_pccomp%co2_pc_ref       )) call p_bcast(coef_pccomp%co2_pc_ref        ,source,comm)
    if (associated(coef_pccomp%n2o_pc_ref       )) call p_bcast(coef_pccomp%n2o_pc_ref        ,source,comm)
    if (associated(coef_pccomp%co_pc_ref        )) call p_bcast(coef_pccomp%co_pc_ref         ,source,comm)
    if (associated(coef_pccomp%ch4_pc_ref       )) call p_bcast(coef_pccomp%ch4_pc_ref        ,source,comm)
    if (associated(coef_pccomp%noise            )) call p_bcast(coef_pccomp%noise             ,source,comm)
    if (associated(coef_pccomp%noise_in         )) call p_bcast(coef_pccomp%noise_in          ,source,comm)
    if (associated(coef_pccomp%ff_ori_chn_in    )) call p_bcast(coef_pccomp%ff_ori_chn_in     ,source,comm)
    if (associated(coef_pccomp%ff_cwn_in        )) call p_bcast(coef_pccomp%ff_cwn_in         ,source,comm)
    if (associated(coef_pccomp%ff_bco_in        )) call p_bcast(coef_pccomp%ff_bco_in         ,source,comm)
    if (associated(coef_pccomp%ff_bcs_in        )) call p_bcast(coef_pccomp%ff_bcs_in         ,source,comm)
    if (associated(coef_pccomp%planck1_in       )) call p_bcast(coef_pccomp%planck1_in        ,source,comm)
    if (associated(coef_pccomp%planck2_in       )) call p_bcast(coef_pccomp%planck2_in        ,source,comm)

#if defined(RTTOV12)
    if (associated(coef_pccomp%eigen)) then
      do run1 = 1,coef_pccomp%fmv_pc_bands
        call p_bcast(coef_pccomp%eigen(run1),source,comm)
      enddo
    endif
    if (associated(coef_pccomp%pcreg)) then
      do run1 = 1,coef_pccomp%fmv_pc_bands
         do run2 = 1,coef_pccomp%fmv_pc_msets
            call p_bcast(coef_pccomp%pcreg(run1,run2),source,comm)
         end do
      enddo
    endif
#else
    if (associated(coef_pccomp%pcreg)) then
      do run1 = 1,coef_pccomp%fmv_pc_sets
        call p_bcast(coef_pccomp%pcreg(run1),source,comm)
      enddo
    endif
#endif
  end subroutine p_bcast_rttov_coef_pccomp

  subroutine p_bcast_rttov_coef_pccomp1(coef_pccomp1,source,comm)
  type(rttov_coef_pccomp1),intent(inout) :: coef_pccomp1
  integer,              intent(in)      :: source
  integer, optional,    intent(in)      :: comm
  !-------------------------------------------------------------------------
  ! Broadcast an rttov_coef_pccomp structure across all available processors
  !-------------------------------------------------------------------------
    integer :: run1
    integer :: dimensions(2)

    call p_bcast_rttov_container (coef_pccomp1, source, comm)

    if (associated(coef_pccomp1%coefficients)) &
      dimensions(1:2) = shape(coef_pccomp1%coefficients)

    call p_bcast(dimensions,source,comm)

    if (mpi_my_proc_id /= source) then
      if (associated(coef_pccomp1%predictindex)) &
        allocate(coef_pccomp1%predictindex(coef_pccomp1%fmv_pc_npred))
      if (associated(coef_pccomp1%coefficients)) &
        allocate(coef_pccomp1%coefficients(dimensions(1),dimensions(2)))
    endif

    if (associated(coef_pccomp1%predictindex)) &
      call p_bcast(coef_pccomp1%predictindex, source, comm)
    if (associated(coef_pccomp1%coefficients)) &
      call p_bcast(coef_pccomp1%coefficients, source, comm)
  end subroutine p_bcast_rttov_coef_pccomp1

#if defined(RTTOV12)
  subroutine p_bcast_rttov_coef_pccomp2(coef_pccomp2,source,comm)
  type(rttov_coef_pccomp2),intent(inout) :: coef_pccomp2
  integer,              intent(in)      :: source
  integer, optional,    intent(in)      :: comm
  !-------------------------------------------------------------------------
  ! Broadcast an rttov_coef_pccomp structure across all available processors
  !-------------------------------------------------------------------------
    integer :: run1
    integer :: dimensions(4)

    call p_bcast_rttov_container (coef_pccomp2, source, comm)

    if (associated(coef_pccomp2%eigenvectors  )) dimensions(1:2) = shape(coef_pccomp2%eigenvectors  )
    if (associated(coef_pccomp2%eigenvectors_t)) dimensions(3:4) = shape(coef_pccomp2%eigenvectors_t)

    call p_bcast(dimensions,source,comm)

    if (mpi_my_proc_id /= source) then
      if (associated(coef_pccomp2%eigenvectors  )) allocate(coef_pccomp2%eigenvectors  (dimensions(1),dimensions(2)))
      if (associated(coef_pccomp2%eigenvectors_t)) allocate(coef_pccomp2%eigenvectors_t(dimensions(1),dimensions(2)))
    endif

    if (associated(coef_pccomp2%eigenvectors  )) call p_bcast(coef_pccomp2%eigenvectors  , source, comm)
    if (associated(coef_pccomp2%eigenvectors_t)) call p_bcast(coef_pccomp2%eigenvectors_t, source, comm)

  end subroutine p_bcast_rttov_coef_pccomp2
#endif

  subroutine p_bcast_rttov_cnt_coef(buffer,source,comm)
  type(rttov_coef), intent(inout)   :: buffer
  integer,          intent(in)      :: source
  integer, optional,intent(in)      :: comm
  !------------------------------------------------------------------
  ! Broadcast an rttov_coef container across all available processors
  !------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#endif
  end subroutine p_bcast_rttov_cnt_coef

  subroutine p_bcast_rttov_cnt_coef_scatt_ir(buffer,source,comm)
  type(rttov_coef_scatt_ir), intent(inout) :: buffer
  integer,          intent(in)             :: source
  integer, optional,intent(in)             :: comm
  !---------------------------------------------------------------------------
  ! Broadcast an rttov_coef_scatt_ir container across all available processors
  !---------------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#endif
  end subroutine p_bcast_rttov_cnt_coef_scatt_ir

  subroutine p_bcast_rttov_cnt_optpar_ir(buffer,source,comm)
  type(rttov_optpar_ir), intent(inout) :: buffer
  integer,          intent(in)         :: source
  integer, optional,intent(in)         :: comm
  !-----------------------------------------------------------------------
  ! Broadcast an rttov_optpar_ir container across all available processors
  !-----------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#endif
  end subroutine p_bcast_rttov_cnt_optpar_ir

#if defined(RTTOV12)
  subroutine p_bcast_rttov_cnt_coef_optpiclb(buffer,source,comm)
  type(rttov_coef_optpiclb), intent(inout) :: buffer
  integer,                   intent(in)    :: source
  integer,         optional, intent(in)    :: comm
  !-----------------------------------------------------------------------
  ! Broadcast an rttov_optpar_ir container across all available processors
  !-----------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#endif
  end subroutine p_bcast_rttov_cnt_coef_optpiclb

  subroutine p_bcast_rttov_cnt_phasefn_int(buffer,source,comm)
  type(rttov_phasefn_int),  intent(inout) :: buffer
  integer,                  intent(in)    :: source
  integer,         optional,intent(in)    :: comm
  !-----------------------------------------------------------------------
  ! Broadcast an rttov_optpar_ir container across all available processors
  !-----------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#endif
  end subroutine p_bcast_rttov_cnt_phasefn_int

  subroutine p_bcast_rttov_cnt_fast_coef(buffer,source,comm)
  type(rttov_fast_coef),  intent(inout) :: buffer
  integer,                intent(in)    :: source
  integer,       optional,intent(in)    :: comm
  !-----------------------------------------------------------------------
  ! Broadcast an rttov_optpar_ir container across all available processors
  !-----------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#endif
  end subroutine p_bcast_rttov_cnt_fast_coef

  subroutine p_bcast_rttov_cnt_fast_coef_gas(buffer,source,comm)
  type(rttov_fast_coef_gas),  intent(inout) :: buffer
  integer,                    intent(in)    :: source
  integer,           optional,intent(in)    :: comm
  !-----------------------------------------------------------------------
  ! Broadcast an rttov_optpar_ir container across all available processors
  !-----------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#endif
  end subroutine p_bcast_rttov_cnt_fast_coef_gas

  subroutine p_bcast_rttov_cnt_nlte_coef(buffer,source,comm)
  type(rttov_nlte_coef),  intent(inout) :: buffer
  integer,                intent(in)    :: source
  integer,       optional,intent(in)    :: comm
  !-----------------------------------------------------------------------
  ! Broadcast an rttov_optpar_ir container across all available processors
  !-----------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#endif
  end subroutine p_bcast_rttov_cnt_nlte_coef

#endif

  subroutine p_bcast_rttov_cnt_coef_pccomp(buffer,source,comm)
  type(rttov_coef_pccomp), intent(inout) :: buffer
  integer,          intent(in)           :: source
  integer, optional,intent(in)           :: comm
  !-------------------------------------------------------------------------
  ! Broadcast an rttov_coef_pccomp container across all available processors
  !-------------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#endif
  end subroutine p_bcast_rttov_cnt_coef_pccomp

  subroutine p_bcast_rttov_cnt_coef_pccomp1(buffer,source,comm)
  type(rttov_coef_pccomp1), intent(inout) :: buffer
  integer,          intent(in)            :: source
  integer, optional,intent(in)            :: comm
  !--------------------------------------------------------------------------
  ! Broadcast an rttov_coef_pccomp1 container across all available processors
  !--------------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#endif
  end subroutine p_bcast_rttov_cnt_coef_pccomp1

  subroutine p_bcast_rttov_complex_1d(buffer,source,comm)
  complex(jprb),    intent(inout)   :: buffer(:)
  integer,          intent(in)      :: source
  integer, optional,intent(in)      :: comm
  !-----------------------------------------------------------
  ! Broadcast a complex vector across all available processors
  !-----------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_DOUBLE_COMPLEX, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
  end subroutine p_bcast_rttov_complex_1d

#if defined(RTTOV12)
  subroutine p_bcast_rttov_cnt_coef_pccomp2(buffer,source,comm)
  type(rttov_coef_pccomp2), intent(inout) :: buffer
  integer,          intent(in)            :: source
  integer, optional,intent(in)            :: comm
  !--------------------------------------------------------------------------
  ! Broadcast an rttov_coef_pccomp2 container across all available processors
  !--------------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#endif
  end subroutine p_bcast_rttov_cnt_coef_pccomp2
#endif

  subroutine p_bcast_rttov_char_1d(buffer,source,comm)
  character(len=4), intent(inout)   :: buffer(:)
  integer,          intent(in)      :: source
  integer, optional,intent(in)      :: comm
  !------------------------------------------------------------
  ! Broadcast an integer vector across all available processors
  !------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF) || defined(RTTOV_IFC_USE_MPI_DACE)
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_CHARACTER, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
#endif
  end subroutine p_bcast_rttov_char_1d

#ifdef RTTOV_IFC_USE_MPIF
  subroutine p_bcast_rttov_integer_1d(buffer,source,comm)
  integer(jpim),    intent(inout)   :: buffer(:)
  integer,          intent(in)      :: source
  integer, optional,intent(in)      :: comm
  !------------------------------------------------------------
  ! Broadcast an integer vector across all available processors
  !------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_INTEGER, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
  end subroutine p_bcast_rttov_integer_1d

  subroutine p_bcast_rttov_integer_2d(buffer,source,comm)
  integer(jpim),    intent(inout)   :: buffer(:,:)
  integer,          intent(in)      :: source
  integer, optional,intent(in)      :: comm
  !--------------------------------------------------------------
  ! Broadcast a real rank-2 array across all available processors
  !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_INTEGER, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
  end subroutine p_bcast_rttov_integer_2d

  subroutine p_bcast_rttov_integer_3d(buffer,source,comm)
  integer(jpim),    intent(inout)   :: buffer(:,:,:)
  integer,          intent(in)      :: source
  integer, optional,intent(in)      :: comm
  !--------------------------------------------------------------
  ! Broadcast a real rank-2 array across all available processors
  !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_INTEGER, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
  end subroutine p_bcast_rttov_integer_3d

  subroutine p_bcast_rttov_integer_4d(buffer,source,comm)
  integer(jpim),    intent(inout)   :: buffer(:,:,:,:)
  integer,          intent(in)      :: source
  integer, optional,intent(in)      :: comm
  !--------------------------------------------------------------
  ! Broadcast a real rank-2 array across all available processors
  !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_INTEGER, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
  end subroutine p_bcast_rttov_integer_4d

  subroutine p_bcast_rttov_real_1d(buffer,source,comm)
  real(jprb),    intent(inout)   :: buffer(:)
  integer,          intent(in)      :: source
  integer, optional,intent(in)      :: comm
  !--------------------------------------------------------
  ! Broadcast a real vector across all available processors
  !--------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_DOUBLE_PRECISION, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
  end subroutine p_bcast_rttov_real_1d

  subroutine p_bcast_rttov_real_2d(buffer,source,comm)
  real(jprb),    intent(inout)   :: buffer(:,:)
  integer,          intent(in)      :: source
  integer, optional,intent(in)      :: comm
  !--------------------------------------------------------------
  ! Broadcast a real rank-2 array across all available processors
  !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_DOUBLE_PRECISION, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
  end subroutine p_bcast_rttov_real_2d

  subroutine p_bcast_rttov_real_3d(buffer,source,comm)
  real(jprb),    intent(inout)   :: buffer(:,:,:)
  integer,          intent(in)      :: source
  integer, optional,intent(in)      :: comm
  !--------------------------------------------------------------
  ! Broadcast a real rank-3 array across all available processors
  !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_DOUBLE_PRECISION, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
  end subroutine p_bcast_rttov_real_3d

  subroutine p_bcast_rttov_real_4d(buffer,source,comm)
  real(jprb),    intent(inout)   :: buffer(:,:,:,:)
  integer,          intent(in)      :: source
  integer, optional,intent(in)      :: comm
  !--------------------------------------------------------------
  ! Broadcast a real rank-3 array across all available processors
  !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_DOUBLE_PRECISION, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
  end subroutine p_bcast_rttov_real_4d

  subroutine p_bcast_rttov_bool(buffer,source,comm)
  logical,          intent(inout)   :: buffer
  integer,          intent(in)      :: source
  integer, optional,intent(in)      :: comm
  !------------------------------------------------------------
  ! Broadcast an integer vector across all available processors
  !------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,1, MPI_LOGICAL, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
  end subroutine p_bcast_rttov_bool
#endif /* RTTOV_IFC_USE_MPIF */

#if defined(RTTOV12)

  subroutine p_bcast_rttov_sinteger_4d(buffer,source,comm)
  integer(jpis),    intent(inout)   :: buffer(:,:,:,:)
  integer,          intent(in)      :: source
  integer, optional,intent(in)      :: comm
  !--------------------------------------------------------------
  ! Broadcast a small integer rank-4 array across all available processors
  !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
        source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
  end subroutine p_bcast_rttov_sinteger_4d

#ifndef RADSHARE
 !--------------------------------
 ! RTTOV IFC MPI Transfer routines
 !--------------------------------
 subroutine p_bcast_rttov_atlas (atlas, source, comm)
   type(rttov_emis_atlas_data), intent(inout) :: atlas
   integer,                     intent(in)    :: source
   integer,           optional, intent(in)    :: comm
   !-----------------------------------------------------------------------------
   ! Broadcast an rttov_emis_atlas_data structure across all available processors
   !-----------------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
   integer :: lcom, errorcode

   lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

   call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
        source, lcom, errorcode)

   if (errorcode /= MPI_SUCCESS) then
     print *, 'MPI ERROR in MPI_Bcast: ', errorcode
     stop 'MPI ERROR'
   endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)

   call p_bcast(atlas% init,     source, comm)
   call p_bcast(atlas% is_mw,    source, comm)
   call p_bcast(atlas% atlas_id, source, comm)
#endif

   if (.not. atlas% is_mw) return

   if (atlas% atlas_id == 1) then
     call p_bcast(atlas% telsem2_atlas, source, comm)
   else
     call p_bcast(atlas% cnrm_mw_atlas, source, comm)
   endif

 end subroutine p_bcast_rttov_atlas

 subroutine p_bcast_rttov_telsem(telsem, source, comm)
   type(telsem2_atlas_data) ,intent(inout) :: telsem
   integer,                  intent(in)    :: source
   integer, optional,        intent(in)    :: comm
   !--------------------------------------------------------------------------
   ! Broadcast a telsem2_atlas_data structure across all available processors
   !--------------------------------------------------------------------------
   integer :: dimensions(13)

   call p_bcast_rttov_container (telsem, source, comm)

   dimensions(:) = -1

   if (associated(telsem% ncells))         dimensions(1:1)   = shape(telsem% ncells)
   if (associated(telsem% firstcell))      dimensions(2:2)   = shape(telsem% firstcell)
   if (associated(telsem% emis))           dimensions(3:4)   = shape(telsem% emis)
   if (associated(telsem% correl))         dimensions(5:7)   = shape(telsem% correl)
   if (associated(telsem% emis_err))       dimensions(8:9)   = shape(telsem% emis_err)
   if (associated(telsem% class1))         dimensions(10:10) = shape(telsem% class1)
   if (associated(telsem% class2))         dimensions(11:11) = shape(telsem% class2)
   if (associated(telsem% cellnum))        dimensions(12:12) = shape(telsem% cellnum)
   if (associated(telsem% correspondance)) dimensions(13:13) = shape(telsem% correspondance)

   call p_bcast(dimensions,source,comm)

   if (mpi_my_proc_id /= source) then
     if (associated(telsem% ncells))         allocate(telsem% ncells    (dimensions(1)))
     if (associated(telsem% firstcell))      allocate(telsem% firstcell (dimensions(2)))
     if (associated(telsem% emis))           allocate(telsem% emis      (dimensions(3),dimensions(4)))
     if (associated(telsem% ncells))         allocate(telsem% correl    (dimensions(5),dimensions(6),dimensions(7)))
     if (associated(telsem% emis_err))       allocate(telsem% emis_err  (dimensions(8),dimensions(9)))
     if (associated(telsem% class1))         allocate(telsem% class1    (dimensions(10)))
     if (associated(telsem% class2))         allocate(telsem% class2    (dimensions(11)))
     if (associated(telsem% cellnum))        allocate(telsem% cellnum   (dimensions(12)))
     if (associated(telsem% correspondance)) allocate(telsem% correspondance (dimensions(13)))
   endif

   if (associated(telsem% ncells))          call p_bcast(telsem% ncells,source,comm)
   if (associated(telsem% firstcell))       call p_bcast(telsem% firstcell,source,comm)
   if (associated(telsem% emis))            call p_bcast(telsem% emis,source,comm)
   if (associated(telsem% correl))          call p_bcast(telsem% correl,source,comm)
   if (associated(telsem% emis_err))        call p_bcast(telsem% emis_err,source,comm)
   if (associated(telsem% class1))          call p_bcast(telsem% class1,source,comm)
   if (associated(telsem% class2))          call p_bcast(telsem% class2,source,comm)
   if (associated(telsem% cellnum))         call p_bcast(telsem% cellnum,source,comm)
   if (associated(telsem% correspondance))  call p_bcast(telsem% correspondance,source,comm)

 end subroutine p_bcast_rttov_telsem

 subroutine p_bcast_rttov_cnrm(cnrm, source, comm)
   type(cnrm_mw_atlas_data) ,intent(inout) :: cnrm
   integer,                  intent(in)    :: source
   integer, optional,        intent(in)    :: comm
   !--------------------------------------------------------------------------
   ! Broadcast a cnrm2_atlas_data structure across all available processors
   !--------------------------------------------------------------------------
   integer :: dimensions(9)

   call p_bcast_rttov_container (cnrm, source, comm)

   dimensions(:) = -1

   if (associated(cnrm% emissivity))     dimensions(1:4)   = shape(cnrm% emissivity)
   if (associated(cnrm% frequencies))    dimensions(5:5)   = shape(cnrm% frequencies)
   if (associated(cnrm% polarisation))   dimensions(6:6)   = shape(cnrm% polarisation)
   if (associated(cnrm% freq_range_min)) dimensions(7:7)   = shape(cnrm% freq_range_min)
   if (associated(cnrm% freq_range_max)) dimensions(8:8)   = shape(cnrm% freq_range_max)
   if (associated(cnrm% zenith_angles))  dimensions(9:9)   = shape(cnrm% zenith_angles)

   call p_bcast(dimensions,source,comm)

   if (mpi_my_proc_id /= source) then
     if (associated(cnrm% emissivity))     allocate(cnrm% emissivity(dimensions(1),dimensions(2),dimensions(3),dimensions(4)))
     if (associated(cnrm% frequencies))    allocate(cnrm% frequencies(dimensions(5)))
     if (associated(cnrm% polarisation))   allocate(cnrm% polarisation(dimensions(6)))
     if (associated(cnrm% freq_range_min)) allocate(cnrm% freq_range_min(dimensions(7)))
     if (associated(cnrm% freq_range_max)) allocate(cnrm% freq_range_max(dimensions(8)))
     if (associated(cnrm% zenith_angles))  allocate(cnrm% zenith_angles(dimensions(9)))
   endif

   if (associated(cnrm% emissivity))      call p_bcast(cnrm% emissivity,source,comm)
   if (associated(cnrm% frequencies))     call p_bcast(cnrm% frequencies,source,comm)
   if (associated(cnrm% polarisation))    call p_bcast(cnrm% polarisation,source,comm)
   if (associated(cnrm% freq_range_min))  call p_bcast(cnrm% freq_range_min,source,comm)
   if (associated(cnrm% freq_range_max))  call p_bcast(cnrm% freq_range_max,source,comm)
   if (associated(cnrm% zenith_angles))   call p_bcast(cnrm% zenith_angles,source,comm)
 end subroutine p_bcast_rttov_cnrm

 subroutine p_bcast_rttov_cnt_telsem(buffer,source,comm)
   type(telsem2_atlas_data), intent(inout)   :: buffer
   integer,                  intent(in)      :: source
   integer,         optional,intent(in)      :: comm
   !------------------------------------------------------------------
   ! Broadcast an rttov_coef container across all available processors
   !------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
   integer :: lcom, errorcode

   lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

   call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
        source, lcom, errorcode)

   if (errorcode /= MPI_SUCCESS) then
     print *, 'MPI ERROR in MPI_Bcast: ', errorcode
     stop 'MPI ERROR'
   endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)
   call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#endif
 end subroutine p_bcast_rttov_cnt_telsem

 subroutine p_bcast_rttov_cnt_cnrm(buffer,source,comm)
   type(cnrm_mw_atlas_data), intent(inout)   :: buffer
   integer,                  intent(in)      :: source
   integer,         optional,intent(in)      :: comm
   !------------------------------------------------------------------
   ! Broadcast an rttov_coef container across all available processors
   !------------------------------------------------------------------
#if defined(RTTOV_IFC_USE_MPIF)
   integer :: lcom, errorcode

   lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

   call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
        source, lcom, errorcode)

   if (errorcode /= MPI_SUCCESS) then
     print *, 'MPI ERROR in MPI_Bcast: ', errorcode
     stop 'MPI ERROR'
   endif
#elif defined(RTTOV_IFC_USE_MPI_DACE)
   call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#endif
 end subroutine p_bcast_rttov_cnt_cnrm
#endif /* RADSHARE */

#endif /* RTTOV12 */

#endif /* _RTTOV_DO_DISTRIBCOEF */

!------------------------------------------------------------------------------
  subroutine rttov_print_profiles(unit)
    integer, intent(in) :: unit
    integer             :: j
    do j = 1, size(profiles)
      call print_profile(j,unit=unit)
    end do
  end subroutine rttov_print_profiles
!------------------------------------------------------------------------------

  subroutine print_profile(i, unit, form)
  !--------------------------------------------
  ! subroutine print_profile for debugging only
  !--------------------------------------------
    integer, intent(in)           :: i
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: form

#if defined(RTTOV12)
    real(jprb), parameter :: RDRD = mh2o/mair
    real(jprb) :: q, ppmv_dry
#endif
    integer :: iunit, j
    logical :: form_loc

    if (present(unit)) then
      iunit=unit
    else
      iunit=stdout
    end if
    form_loc=.true.
    if (present(form)) form_loc = form

    if (form_loc) then
      if (i > size(profiles)) then
        write(iunit,*) 'rttov input profile: invalid index'
        return
      end if
      write(iunit,*) 'rttov input profile ',i,trim(profiles(i)%id)
      write(iunit,*) profiles(i)%nlevels
      if (associated(profiles(i)%p) .and. associated(profiles(i)%t) .and. associated(profiles(i)%q)) then
#if defined(RTTOV12)
        write(iunit,'(3x,9x,"p",9x,"t",13x,"q",4x,"ppmv(dry)")')
        do j = 1, size(profiles(i)%p)
          select case (profiles(i)%gas_units)
          case(gas_unit_specconc)
            q = profiles(i)%q(j)
            ppmv_dry = (1.d6 * q) / (RDRD * (1 - q))
          case(gas_unit_ppmvdry)
            ppmv_dry = profiles(i)%q(j)
            q = (RDRD * ppmv_dry) / (1.d6 + RDRD * ppmv_dry)
          case(gas_unit_ppmv)
            ppmv_dry = profiles(i)%q(j) / (1._jprb - 1.E-6_jprb * profiles(i)%q(j))
            q = (RDRD * ppmv_dry) / (1.d6 + RDRD * ppmv_dry)
          end select
          write(iunit,'(I3,2(1x,f9.4),2(1x,E13.6))') j, profiles(i)%p(j), profiles(i)%t(j), q, ppmv_dry
        end do
#else
        write(iunit,'(3x,9x,"p",9x,"t",13x,"q")')
        do j = 1, size(profiles(i)%p)
          write(iunit,'(I3,2(1x,f9.4),1x,E13.6)') j, profiles(i)%p(j), profiles(i)%t(j), profiles(i)%q(j)
        end do
#endif
      end if
      write(iunit,*) 'cloud'
      !    if (associated(profiles(i)%cloud)) write(iunit,'(6(1x,e25.17))'),profiles(i)%cloud(:,:)
      if (associated(profiles(i)%cloud)) write(iunit,*) profiles(i)%cloud(:,:)
      write(iunit,*) 'cfrac'
      !    if (associated(profiles(i)%cfrac)) write(iunit,'(6(1x,e25.17))'),profiles(i)%cfrac(:,:)
#if defined(RTTOV12)
      if (associated(profiles(i)%cfrac)) write(iunit,*) profiles(i)%cfrac(:)
#else
      if (associated(profiles(i)%cfrac)) write(iunit,*) profiles(i)%cfrac(:,:)
#endif
      write(iunit,*) 'idg=',profiles(i)%idg
#if defined(RTTOV12)
      write(iunit,*) 'ice_scheme=',profiles(i)%ice_scheme
#else
      write(iunit,*) 'ice_scheme=',profiles(i)%ish
#endif
      write(iunit,*) 'skin'
      write(iunit,*) 'surftype=',profiles(i)%skin%surftype
      write(iunit,*) 'watertype=',profiles(i)%skin%watertype
      write(iunit,*) 't=',profiles(i)%skin%t
      write(iunit,*) 'fastem=',profiles(i)%skin%fastem
      write(iunit,*) '2m'
      write(iunit,*) 't=',profiles(i)%s2m%t
      write(iunit,*) 'q=',profiles(i)%s2m%q
      write(iunit,*) 'o=',profiles(i)%s2m%o
      write(iunit,*) 'p=',profiles(i)%s2m%p
      write(iunit,*) 'u=',profiles(i)%s2m%u
      write(iunit,*) 'v=',profiles(i)%s2m%v
      write(iunit,*) 'wfetch =',profiles(i)%s2m%wfetc
      write(iunit,*) 'zenangle=',profiles(i)%zenangle
      write(iunit,*) 'azangle=',profiles(i)%azangle
      write(iunit,*) 'sunzenangle=',profiles(i)%sunzenangle
      write(iunit,*) 'sunazangle=',profiles(i)%sunazangle
      write(iunit,*) 'elevation=',profiles(i)%elevation
      write(iunit,*) 'latitude=',profiles(i)%latitude
      write(iunit,*) 'longitude=',profiles(i)%longitude
      write(iunit,*) 'Be=',profiles(i)%Be
      write(iunit,*) 'cosbk=',profiles(i)%cosbk
      write(iunit,*) 'ctp=',profiles(i)%ctp
      write(iunit,*) 'cfraction=',profiles(i)%cfraction
    else
      if (i > size(profiles)) then
        write(iunit) 'rttov input profile: invalid index'
        return
      end if
      write(iunit) profiles(i)%nlevels
      if (associated(profiles(i)%p) .and. associated(profiles(i)%t) .and. associated(profiles(i)%q)) then
        do j = 1, size(profiles(i)%p)
          write(iunit) j, profiles(i)%p(j), profiles(i)%t(j), profiles(i)%q(j)
        end do
      end if
      if (associated(profiles(i)%cloud)) write(iunit) profiles(i)%cloud(:,:)
#if defined(RTTOV12)
      if (associated(profiles(i)%cfrac)) write(iunit) profiles(i)%cfrac(:)
#else
      if (associated(profiles(i)%cfrac)) write(iunit) profiles(i)%cfrac(:,:)
#endif
      write(iunit) profiles(i)%idg
#if defined(RTTOV12)
      write(iunit) profiles(i)%ice_scheme
#else
      write(iunit) profiles(i)%ish
#endif
      write(iunit) profiles(i)%skin%surftype
      write(iunit) profiles(i)%skin%watertype
      write(iunit) profiles(i)%skin%t
      write(iunit) profiles(i)%skin%fastem
      write(iunit) profiles(i)%s2m%t
      write(iunit) profiles(i)%s2m%q
      write(iunit) profiles(i)%s2m%o
      write(iunit) profiles(i)%s2m%p
      write(iunit) profiles(i)%s2m%u
      write(iunit) profiles(i)%s2m%v
      write(iunit) profiles(i)%s2m%wfetc
      write(iunit) profiles(i)%zenangle
      write(iunit) profiles(i)%azangle
      write(iunit) profiles(i)%sunzenangle
      write(iunit) profiles(i)%sunazangle
      write(iunit) profiles(i)%elevation
      write(iunit) profiles(i)%latitude
      write(iunit) profiles(i)%longitude
      write(iunit) profiles(i)%Be
      write(iunit) profiles(i)%cosbk
      write(iunit) profiles(i)%ctp
      write(iunit) profiles(i)%cfraction
    end if
  end subroutine print_profile

  subroutine read_profile(i, iunit, ind)
  !--------------------------------------------
  ! subroutine print_profile for debugging only
  !--------------------------------------------
    integer, intent(in) :: i
    integer, intent(in) :: iunit
    integer, intent(in), optional :: ind

    integer :: j,k,ii

    do
      read(iunit) j,ii
      read(iunit) profiles(i)%nlevels
      if (associated(profiles(i)%p) .and. associated(profiles(i)%t) .and. associated(profiles(i)%q)) then
        do j = 1, size(profiles(i)%p)
          read(iunit) k, profiles(i)%p(j), profiles(i)%t(j), profiles(i)%q(j)
          if (k /= j) stop
        end do
      end if
      if (associated(profiles(i)%cloud)) read(iunit) profiles(i)%cloud(:,:)
#if defined(RTTOV12)
      if (associated(profiles(i)%cfrac)) read(iunit) profiles(i)%cfrac(:)
#else
      if (associated(profiles(i)%cfrac)) read(iunit) profiles(i)%cfrac(:,:)
#endif
      read(iunit) profiles(i)%idg
#if defined(RTTOV12)
      read(iunit) profiles(i)%ice_scheme
#else
      read(iunit) profiles(i)%ish
#endif
      read(iunit) profiles(i)%skin%surftype
      read(iunit) profiles(i)%skin%watertype
      read(iunit) profiles(i)%skin%t
      read(iunit) profiles(i)%skin%fastem
      read(iunit) profiles(i)%s2m%t
      read(iunit) profiles(i)%s2m%q
      read(iunit) profiles(i)%s2m%o
      read(iunit) profiles(i)%s2m%p
      read(iunit) profiles(i)%s2m%u
      read(iunit) profiles(i)%s2m%v
      read(iunit) profiles(i)%s2m%wfetc
      read(iunit) profiles(i)%zenangle
      read(iunit) profiles(i)%azangle
      read(iunit) profiles(i)%sunzenangle
      read(iunit) profiles(i)%sunazangle
      read(iunit) profiles(i)%elevation
      read(iunit) profiles(i)%latitude
      read(iunit) profiles(i)%longitude
      read(iunit) profiles(i)%Be
      read(iunit) profiles(i)%cosbk
      read(iunit) profiles(i)%ctp
      read(iunit) profiles(i)%cfraction

      if (.not.present(ind)) then
        exit
      else
        if (ii == ind) exit
      end if
    end do
  end subroutine read_profile

  subroutine print_transmission(i, unit)
  !--------------------------------------------
  ! subroutine print_profile for debugging only
  !--------------------------------------------
    integer, intent(in)           :: i
    integer, intent(in), optional :: unit

    integer :: iunit

    if (present(unit)) then
      iunit=unit
    else
      iunit=stdout
    end if
    write(iunit,*) 'rttov transmission output ',i
    write(iunit,*) 'tau_total   = ',transmission%tau_total (i)
    write(iunit,*) 'tau_levels  = ',transmission%tau_levels(:,i)
  end subroutine print_transmission

  subroutine print_radiance(i, unit)
  !--------------------------------------------
  ! subroutine print_profile for debugging only
  !--------------------------------------------
    integer, intent(in)           :: i
    integer, intent(in), optional :: unit

    integer :: iunit

    if (present(unit)) then
      iunit=unit
    else
      iunit=stdout
    end if
    write(iunit,*) 'rttov radiance output ',i
    write(iunit,*) 'clear      = ',radiance%clear(i)
    write(iunit,*) 'total      = ',radiance%total(i)
    write(iunit,*) 'bt_clear   = ',radiance%bt_clear(i)
    write(iunit,*) 'bt         = ',radiance%bt(i)
#if defined(RTTOV12)
    write(iunit,*) 'refl_clear = ',radiance%refl_clear(i)
    write(iunit,*) 'refl       = ',radiance%refl(i)
#elif defined(RTTOV10)
    write(iunit,*) 'refl_clear = ',radiance%reflclear(i)
#endif
    write(iunit,*) 'overcast   = ',radiance%overcast(i,:)
    write(iunit,*) 'cloudy     = ',radiance%cloudy(i)
 end subroutine print_radiance

#endif /* !NO_RTTOV */

end module mo_rttov_ifc
