!
!+ Basis for RTTOV interface modules
!
MODULE mo_rtifc_base
!
! Description:
!   This module contains basic routines, parameters and variables.
!   for the interface modules. This might be stuff, that does not
!   depend on the RTTOV version, or DWD specific features, that will
!   not change (quite likely) in future RTTOV versions.

!
! Current Code Owner: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    fax:   +49 69 8062 3721
!    email:  robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_X         2020/XX/XX Robin Faulwetter
!  initial version

! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Robin Faulwetter  DWD 2020-05-11 Initial release
!=================================================


!---------------------
! MACRO SETTINGS BEGIN
!---------------------

#include "mo_rtifc_macros.incf"

  !-------------
  ! Modules used
  !-------------
  
  use iso_fortran_env,    only: stdout => output_unit, &!
                                iostat_end              !

#if (_RTTOV_VERSION > 0)
  use rttov_const,        only: version,               &!
                                release,               &!
                                minor_version,         &!
                                fastem_sp,             &! max. number of fastem surface parameters
                                qmin_rttov => qmin,    &! minimum allowed water vapour
                                qmax_rttov => qmax,    &! maximum allowed water vapour
                                tmin_rttov => tmin,    &! minimum allowed temperature
                                tmax_rttov => tmax,    &! maximum allowed temperature
                                min_od,                &!
                                def_gas_unit =>        &
                                   gas_unit_specconc
  use parkind1,           only: jprb, jpim, jplm, jpis
#endif

#if defined(_RTTOV_GOD)
  use rttov_types,        only: rttov_god_par
  use rttov_god,          only: rttov_o2god
#endif

#if defined(_RTIFC_DISTRIBCOEF) && defined(HAVE_MPI_MOD) 
  ! prefer MPI module over mpif.h
  use mpi
#endif

  implicit none

  public

  !-----------
  ! Interfaces
  !-----------

  ! MPI routines
#if defined(_RTIFC_DISTRIBCOEF)

#if !defined(HAVE_MPI_MOD)
  include "mpif.h"
#endif


#if (_RTTOV_VERSION >= 12) || defined(_RTIFC_USE_MPIF)
  interface p_bcast
#if defined(_RTIFC_USE_MPIF)
    module procedure p_bcast_rttov_integer
    module procedure p_bcast_rttov_integer_1d
    module procedure p_bcast_rttov_integer_2d
    module procedure p_bcast_rttov_integer_3d
    module procedure p_bcast_rttov_integer_4d
    module procedure p_bcast_rttov_real_1d
    module procedure p_bcast_rttov_real_2d
    module procedure p_bcast_rttov_real_3d
    module procedure p_bcast_rttov_real_4d
    module procedure p_bcast_rttov_bool
    module procedure p_bcast_rttov_char_1d
#endif
#if (_RTTOV_VERSION >= 12)
    module procedure p_bcast_rttov_sinteger_4d
    module procedure p_bcast_rttov_sinteger_3d
    module procedure p_bcast_rttov_sinteger_2d
    module procedure p_bcast_rttov_complex_1d
#endif
  end interface
#endif

#endif /* _RTIFC_DISTRIBCOEF */

  !-----------
  ! Parameters
  !-----------

  ! data type precision
  integer,         parameter :: dp              = selected_real_kind(13)     ! Double Precision
  !integer,        parameter :: sp              = selected_real_kind(6)      ! Single Precision
  integer,         parameter :: wp              = dp                         ! Working precision

  ! Values for verbosity
  integer,         parameter :: silent          =  0
  integer,         parameter :: production      =  1
  integer,         parameter :: verbose         =  2
  integer,         parameter :: debug           =  3

  ! Flag variables for output arrays
  integer,         parameter :: OUT_ASB         = 0                          ! all sky brightness temp.
  integer,         parameter :: OUT_CSB         = 1                          ! clear sky brightness temp.
  integer,         parameter :: OUT_ASR         = 2                          ! all sky radiances
  integer,         parameter :: OUT_CSR         = 3                          ! clear sky radiances
  integer,         parameter :: OUT_VIS         = 4                          ! reflectances (all sky and clear sky)

  ! Physical parameters required for weighting function calculation
  real(kind=wp),   parameter :: rd              = 287.05_wp
  real(kind=wp),   parameter :: g               = 9.80665_wp

  ! Possible coefficient level numbers
  integer,         parameter :: levels_rttov(4) = (/ 44, 51, 54, 101 /)

  ! Define variables for no-rttov case, that are used from RTTOV modules otherwise
#if (_RTTOV_VERSION <= 0)
!  integer,         parameter :: jpim            = selected_int_kind(9)       ! standard integer type
!  integer,         parameter :: jprb            = selected_real_kind(13,300) ! standard real type
!  integer,         parameter :: jplm            = kind(.true.)               ! standard logical type
  integer,         parameter :: version         = 0                          ! RTTOV version
  integer,         parameter :: fastem_sp       = 5                          ! Number of FASTEM parameters
  integer,         parameter :: def_gas_unit    = 0
  real(kind=wp),   parameter :: qmin_rttov      = -1._wp
  real(kind=wp),   parameter :: qmax_rttov      = -1._wp
  real(kind=wp),   parameter :: tmin_rttov      = -1._wp
  real(kind=wp),   parameter :: tmax_rttov      = -1._wp
#else
  ! Used from rttov modules
#endif



  ! error codes and messages
  integer,         parameter :: nerr            = 23                        ! number of different error messages
  character(len=100)         :: err_msg(0:nerr)
#define DEF_ERR(codename, code, msg) integer,parameter::codename=code;data err_msg(code)/msg/
  DEF_ERR(NO_ERROR            ,  0, 'No error. Everything okay.')
  DEF_ERR(ERR_ALLOC           ,  1, 'Allocation error.')
  DEF_ERR(ERR_DIM             ,  2, 'Wrong dimension size of input array')
  DEF_ERR(ERR_RTTOV_SETUP     ,  3, 'Error in RTTOV setup')
  DEF_ERR(ERR_CLOUD_AERO_MISSM,  4, 'Cloud/Aerosol class mismatch')
  DEF_ERR(ERR_RTTOV_CALL      ,  5, 'RTTOV call failed')
  DEF_ERR(WARN_RTTOV          ,  6, 'Warning: RTTOV error status /= 0')
  DEF_ERR(ERR_RTTOV_MPI       ,  7, 'MPI error')
  DEF_ERR(ERR_NO_RTTOV_LIB    ,  8, 'No RTTOV library available')
  DEF_ERR(ERR_RTTOV_PREC      ,  9, 'Mismatch of working precision and RTTOV precision')
  DEF_ERR(ERR_CLOUD_INCONS    , 10, 'Cloud cover and cloud water content arrays inconsistent')
  DEF_ERR(ERR_GOD_FILE        , 11, 'Failed to read god_par_file')
  DEF_ERR(ERR_WR_PROF         , 12, 'Failed to write hdf5 profile file')
  DEF_ERR(ERR_INVALID_TSKIN   , 13, 'Some invalid t_surf/T_G/ts_fg')
  DEF_ERR(ERR_INPUT           , 14, 'Invalid/unsupported input')
  DEF_ERR(ERR_NO_ATLAS        , 15, 'No atlas support in current configuration')
  DEF_ERR(ERR_ATLAS_INIT      , 16, 'Atlas was not initialized')
  DEF_ERR(ERR_INVALID_INSTR   , 17, 'Invalid instrument (e.g. not supported by atlas)')
  DEF_ERR(ERR_TRACEGAS_INCONS , 18, 'Trace gase options inconsistent with input')
  DEF_ERR(ERR_INVALID_NLEVS   , 19, 'Invalid number of levels')
  DEF_ERR(ERR_INVALID_VERSION , 20, 'Invalid RTTOV version')
  DEF_ERR(ERR_PRECISION_INCONS, 21, 'Real data precisions inconsistent (wp /= jprb).')
  DEF_ERR(ERR_NO_COEFS        , 22, 'No matching coefs found.')
  DEF_ERR(ERR_MULTIPLE_COEFS  , 23, 'Multiple matching coefs found.')
! DEF_ERR(ERR_                , XX, '')
#undef DEF_ERR



  !--------------------------
  ! Internal module variables
  !--------------------------
  ! Only for internal use
  integer                    :: pe_ifc                    = -1           ! mpi id of this processor
  ! For internal use, might be used by the user, but MUST not be modified by user!
  integer                    :: nlevs_top                 = 0            ! Number of coeff. levels above user levels (0 or 1)

  !------------------------------------
  ! External module variables
  ! Intended to be modified by the user
  !------------------------------------
  ! default profile values
  real(kind=wp)              :: default_wfetch            =  100000._wp ! wind fetch
  real(kind=wp)              :: default_fastem(fastem_sp) = (/3.0_wp,5.0_wp,15.0_wp,0.1_wp,0.3_wp/)
                                                                        ! fastem coefficients relevant for land/ice
  integer                    :: default_watertype         =       1     ! water type (fresh 0/ocean 1)
  real(kind=wp)              :: default_salinity          =       0._wp ! salinity
  real(kind=wp)              :: default_o3_surf           = 0.031438_wp ! o3 surface
  real(kind=wp)              :: default_satazim           =      0.0_wp ! satellite azimuth angle
  real(kind=wp)              :: default_sunzenangle       =      0.0_wp ! solar zenith angle
  real(kind=wp)              :: default_sunazangle        =      0.0_wp ! solar azimuth angle
  real(kind=wp)              :: default_ctp               =    500.0_wp ! cloud top pressure
  real(kind=wp)              :: default_cfraction         =      0.0_wp ! cloud fraction
#if (_RTTOV_VERSION >= 13)
  integer                    :: default_icede_param       =       4     ! Same as default_idg, but for RTTOVv13 !CSt
#else
  integer                    :: default_idg               =       4     ! Scheme for IWC to eff 
                                                                        ! shape of the ice crystals RTTOVv12
#endif
  integer                    :: default_ice_scheme        =       1     ! ice particle scheme
  integer                    :: default_clw_scheme        =       2     ! cloud liquid water scheme
  integer                    :: default_gas_units         = def_gas_unit! default gas unit

  integer                    :: verbosity                 = production
  logical                    :: read1pe                   = .false.      ! Read coeffs.only on I/O PE
                                                  ! (Only effective with -D_RTIFC_DISTRIBCOEF)


  ! T/q hard limits
  real(kind=wp)              :: qmin_ifc                  = qmin_rttov
  real(kind=wp)              :: qmax_ifc                  = qmax_rttov
  real(kind=wp)              :: tmin_ifc                  = tmin_rttov
  real(kind=wp)              :: tmax_ifc                  = tmax_rttov

#if (_RTTOV_VERSION <= 0)  
  real(kind=wp)              :: min_od                    = -1._wp
#else
  ! Used from rttov_const
#endif

  ! check on regularization limits
  logical, save, allocatable :: mask_lims_t(:)                          ! mask for t exceeding limit
  logical, save, allocatable :: mask_lims_q(:)                          ! mask for t exceeding limit
  integer                    :: chk_reg_lims              = 0           ! Check regularization limits in rtifc
                                                                        ! bit1 (1): print results (invalid profiles)
                                                                        ! bit2 (2): set flag for use in calling prog

  ! generalized optical depth
  character(len=300)         :: god_par_file              = ''
  logical                    :: wr_god                    = .false.
  character(len=300)         :: out_path                  = ''
  ! check influence of god smoothing
  real(kind=wp)              :: god_thresh                = 1._wp
  integer                    :: chk_god                   = 0           ! Check influence of god smoothing
                                                                        ! bit1 (1): print results (invalid profiles)
                                                                        ! bit2 (2): set flag for use in calling prog

contains


  subroutine rtifc_check_config(vers, nlevs, status)
    integer, intent(in)  :: vers
    integer, intent(in)  :: nlevs
    integer, intent(out) :: status
    !--------------------------------------------
    ! Check RTTOV version and number of levels.
    ! Additionally set nlevs_top (in check_nlevs)
    !--------------------------------------------

    call check_version(vers, status)
    if (status /= NO_ERROR) return
    call check_nlevs(nlevs, nlevs_top, status)
    if (status /= NO_ERROR) return
#if (_RTTOV_VERSION > 0)
    if (wp /= jprb) then
      status = ERR_PRECISION_INCONS
      return
    end if
#endif
  end subroutine rtifc_check_config


  subroutine check_version(vers, status)
    integer, intent(in)  :: vers
    integer, intent(out) :: status

    if (vers == version) then
      status = NO_ERROR
    else
      status = ERR_INVALID_VERSION
    end if

  end subroutine check_version


  subroutine check_nlevs(nlevs, nlevs_top, status)
    integer, intent(in)  :: nlevs
    integer, intent(out) :: nlevs_top
    integer, intent(out) :: status

    integer :: i

    status = ERR_INVALID_NLEVS

    do i = 1, size(levels_rttov)
      nlevs_top = (levels_rttov(i) - nlevs)
      if (any(nlevs_top == (/0, 1/))) then
        status = NO_ERROR
        return
      end if
    end do

  end subroutine check_nlevs


  function rttov_version() result(vers)
    character(len=11) :: vers
#if (_RTTOV_VERSION <= 0)    
    vers = 'NO_RTTOV'
#else
    write(vers,'("RTTOV",I2.2,".",I1,".",I1)') version, release, minor_version
#endif
  end function rttov_version


  function errmsg(code) result(msg)
    character(len=120)            :: msg
    integer, intent(in)           :: code

    write(msg, '("ERROR (",I4,"):")') code
    if (abs(code) >= lbound(err_msg,1) .and. abs(code) <= ubound(err_msg,1)) then
      msg = trim(msg)//' '//trim(err_msg(abs(code)))
    else
      msg = trim(msg)//' '//'Unknown error'
    end if

  end function errmsg


#if (_RTTOV_VERSION > 0)

  subroutine get_weighting_function(p, t, transm, wf, height)
    real(kind=jprb)           :: p(:)
    real(kind=jprb)           :: t(:)
    real(kind=jprb)           :: transm(:)
    real(kind=jprb), optional :: wf(:)
    real(kind=jprb), optional :: height

    integer :: wf_method = 1  ! 0: Method proposed by Lucio Torrisi
                              ! 1: Leapfrog-like method
                              ! 2: spline interpolation
    integer :: hgt_method = 0 ! 0: Method proposed by Lucio Torrisi
                              ! 1: Maximum of weighting function
                              ! 2: spline based maximum of weighting function
    real(kind=jprb) :: wf_(size(p))
    real(kind=jprb) :: p_, t_, psum, wsum, max_wf
    integer         :: i, i1, i2, i3, nlev

    nlev = size(p)

    select case(wf_method)
    case(0)
      do i = 1, nlev
        if (i==1) then
          i1 = 1
          i2 = 2
          i3 = i1
        else
          i1 = i - 1
          i2 = i
          i3 = i2
        end if
        ! Hydrostatic assumption
        wf_(i) = - g * p(i3) / (rd * t(i3)) * &
             (transm(i1) - transm(i2)) / (p(i1) - p(i2))
      end do
    case(1)
      do i = 1, nlev
        if (i==1) then
          i1 = 1
          i2 = 2
          p_ = (p(i1) * p(i2)) * 0.5
          t_ = (t(i1) * t(i2)) * 0.5
        elseif (i==nlev) then
          i1 = i - 1
          i2 = i
          p_ = (p(i1) * p(i2)) * 0.5
          t_ = (t(i1) * t(i2)) * 0.5
        else
          i1 = i - 1
          i2 = i + 1
          p_ = p(i)
          t_ = t(i)
        end if
        ! Hydrostatic assumption
        wf_(i) = - g * p_ / (rd * t_) * &
             (transm(i1) - transm(i2)) / (p(i1) - p(i2))
      end do
    case(2)
    end select
    
    if (present(wf)) then
      wf(1:nlev) = wf_(1:nlev)
    end if

    if (present(height)) then
      select case(hgt_method)
      case(0)
        max_wf = maxval(wf_(:))
        psum = 0._jprb
        wsum = 0._jprb
        do i = 1, nlev
          if (wf_(i) > 0.8 * max_wf) then
            psum = psum + p(i) * wf_(i)
            wsum = wsum + wf_(i)
          end if
        end do
        if (wsum /= 0._jprb) then
          height = psum / wsum
        else
          height = -999._jprb
        end if
      case(1)
        i1 = maxloc(wf_(:), 1)
        if (i1 > 1 .and. i1 < nlev) then
          psum = 0._jprb
          wsum = 0._jprb
          do i = i1 - 1, i1 + 1
            psum = psum + p(i) * wf_(i)
            wsum = wsum + wf_(i)
          end do
          if (wsum /= 0._jprb) then
            height = psum / wsum
          else
            height = -999._jprb
          end if
        else
          height = p(i1)
        end if
      case(2)
        !> \todo spline-based weighting function
      end select
    end if

  end subroutine get_weighting_function
#endif

#if defined(_RTTOV_GOD)
  subroutine check_god_infl(tau, god, istat, msg, debug)
    real(kind=jprb),     intent(in)           :: tau(:)
    type(rttov_god_par), intent(in)           :: god(:)
    integer,             intent(out)          :: istat
    character(len=*),    intent(out)          :: msg
    logical,             intent(in), optional :: debug
    ! Calculate transmission and optical depth profile (without god smoothing) on the
    ! basis on god-smoothed profiles.
    real(kind=jprb) :: tau0(size(tau))
    real(kind=jprb) :: infl(size(tau)-1), infl0(size(tau)-1)
    real(kind=jprb) :: s_infl, od, od0, od0_tot, rdiff
    integer :: nl, i, l
    logical :: l_debug = .false.

    if (present(debug)) then
      l_debug = debug
    else
      l_debug = .false.
    end if

    istat = 0
    msg   = ''

    nl = size(tau)-1

    s_infl = 0._jprb
    od0_tot = 0._jprb
    do i = 1,nl
      if (tau(i) > tau(i+1)) then
        od  = log(tau(i+1)/tau(i))
        od0 = min(rttov_o2god(od, p=god(i)), 0._jprb)
      else
        od  = 0._jprb
        od0 = 0._jprb
      end if
      tau0(i) = exp(od0_tot)
      od0_tot = od0_tot + od0
      if (l_debug) print*,i,od,od0,od0_tot, tau(i), tau0(i)
    end do
    tau0(nl+1) = exp(od0_tot)

    do i = 1,nl
      infl (i) = tau (i) * (tau (i) - tau (i+1))
      infl0(i) = tau0(i) * (tau0(i) - tau0(i+1))
      s_infl = s_infl + infl0(i)
    end do

    if (l_debug) print*,i,s_infl
    if (s_infl > 0._jprb) then
      do i = 1, nl
        rdiff = (infl(i)-infl0(i))/s_infl
        if (rdiff > god_thresh) then
          istat = 1
          l = len_trim(msg) + 1
          if (l > 1) then
            msg = trim(msg)//','
            l = l + 1
          end if
          write(msg(l:),'("layer ",I3," d(infl)/infl=",f6.4)') i, rdiff
        end if
        if (l_debug) print*,i,infl(i),infl0(i),infl(i)-infl0(i),(infl(i)-infl0(i))/s_infl
      end do
    else
      istat = 2
      write(msg,*) 'Unable to calculate influence, total influence==0 !!!'
    end if

  end subroutine check_god_infl  
#endif

!==============================
#if defined(_RTIFC_DISTRIBCOEF)
!==============================

#if defined(_RTIFC_USE_MPIF)
  subroutine p_bcast_rttov_integer(buffer,source,comm)
    integer(jpim),    intent(inout)   :: buffer
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !------------------------------------------------------------
    ! Broadcast an integer vector across all available processors
    !------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,1, MPI_INTEGER, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
  end subroutine p_bcast_rttov_integer


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


  subroutine p_bcast_rttov_char_1d(buffer,source,comm)
    character(len=4), intent(inout)   :: buffer(:)
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !------------------------------------------------------------
    ! Broadcast an integer vector across all available processors
    !------------------------------------------------------------
    integer :: lcom, errorcode
    
    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_CHARACTER, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) then
      print *, 'MPI ERROR in MPI_Bcast: ', errorcode
      stop 'MPI ERROR'
    endif
  end subroutine p_bcast_rttov_char_1d

#endif /* _RTIFC_USE_MPIF */

#if (_RTTOV_VERSION >= 12)
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


  subroutine p_bcast_rttov_sinteger_3d(buffer,source,comm)
    integer(jpis),    intent(inout)   :: buffer(:,:,:)
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
  end subroutine p_bcast_rttov_sinteger_3d


  subroutine p_bcast_rttov_sinteger_2d(buffer,source,comm)
    integer(jpis),    intent(inout)   :: buffer(:,:)
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
  end subroutine p_bcast_rttov_sinteger_2d


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
#endif


!==============================
#endif /* _RTIFC_DISTRIBCOEF */
!==============================

end module mo_rtifc_base
