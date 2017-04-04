module mo_psrad_general
  use, intrinsic :: iso_c_binding, only: c_float, c_double, c_long, c_int
  !use mo_kind, only : kdp => dp, ki8 => i8, ki4 => i4

  implicit none
  public
  integer, parameter :: &
    dp = c_double, wp = c_double, &
    i8 = c_long, i4 = c_int, &
    mg = 16, &!< number of original g-intervals per spectral band
    jpband = 29, & !< number of last band (lw and sw share band 16)
    nbndsw = 14, & !< number of spectral bands in sw model
    nbndlw = 16, & !< number of spectral bands in lw model
    maxperband = 16, & !< max num gpts per band
    max_ref = 10, &
    max_minor_species = 5, &
    ngptsw = 112, & !< total number of gpts 
    ngptlw = 140, & !< total number of reduced g-intervals for rrtmg_lw
    ngas = 8, ih2o  = 1, ico2  = 2, ich4  = 3, &
    io2   = 4, io3   = 5, in2o = 6, ico = 7, in2 = 8, &
    ncfc = 4, &! maximum number of cross-section molecules (cfcs)
    cfc_offset = ngas, ifirst_cfc = cfc_offset + 1, &
    iccl4 = cfc_offset+1, &
    icfc11 = cfc_offset+2, &
    icfc12 = cfc_offset+3, &
    icfc22 = cfc_offset+4, &
    ilast_cfc = cfc_offset+ncfc, &
    nreact = 6, ih2oco2 = 1, ih2oo3 = 2, ih2on2o = 3, ih2och4 = 4, &
    in2oco2 = 5, io3co2 = 6
  !
  ! These pressures are chosen such that the ln of the first pressure
  ! has only a few non-zero digits (i.e. ln(PREF(1)) = 6.96000) and
  ! each subsequent ln(pressure) differs from the previous one by 0.2.
  !
  real(wp), parameter :: pref(59) = (/ &
      1.05363e+03_wp, 8.62642e+02_wp, 7.06272e+02_wp, 5.78246e+02_wp, &
      4.73428e+02_wp, 3.87610e+02_wp, 3.17348e+02_wp, 2.59823e+02_wp, &
      2.12725e+02_wp, 1.74164e+02_wp, 1.42594e+02_wp, 1.16746e+02_wp, &
      9.55835e+01_wp, 7.82571e+01_wp, 6.40715e+01_wp, 5.24573e+01_wp, &
      4.29484e+01_wp, 3.51632e+01_wp, 2.87892e+01_wp, 2.35706e+01_wp, &
      1.92980e+01_wp, 1.57998e+01_wp, 1.29358e+01_wp, 1.05910e+01_wp, &
      8.67114e+00_wp, 7.09933e+00_wp, 5.81244e+00_wp, 4.75882e+00_wp, &
      3.89619e+00_wp, 3.18993e+00_wp, 2.61170e+00_wp, 2.13828e+00_wp, &
      1.75067e+00_wp, 1.43333e+00_wp, 1.17351e+00_wp, 9.60789e-01_wp, &
      7.86628e-01_wp, 6.44036e-01_wp, 5.27292e-01_wp, 4.31710e-01_wp, &
      3.53455e-01_wp, 2.89384e-01_wp, 2.36928e-01_wp, 1.93980e-01_wp, &
      1.58817e-01_wp, 1.30029e-01_wp, 1.06458e-01_wp, 8.71608e-02_wp, &
      7.13612e-02_wp, 5.84256e-02_wp, 4.78349e-02_wp, 3.91639e-02_wp, &
      3.20647e-02_wp, 2.62523e-02_wp, 2.14936e-02_wp, 1.75975e-02_wp, &
      1.44076e-02_wp, 1.17959e-02_wp, 9.65769e-03_wp /), &
    preflog(59) = (/ &
      6.9600e+00_wp, 6.7600e+00_wp, 6.5600e+00_wp, 6.3600e+00_wp, &
      6.1600e+00_wp, 5.9600e+00_wp, 5.7600e+00_wp, 5.5600e+00_wp, &
      5.3600e+00_wp, 5.1600e+00_wp, 4.9600e+00_wp, 4.7600e+00_wp, &
      4.5600e+00_wp, 4.3600e+00_wp, 4.1600e+00_wp, 3.9600e+00_wp, &
      3.7600e+00_wp, 3.5600e+00_wp, 3.3600e+00_wp, 3.1600e+00_wp, &
      2.9600e+00_wp, 2.7600e+00_wp, 2.5600e+00_wp, 2.3600e+00_wp, &
      2.1600e+00_wp, 1.9600e+00_wp, 1.7600e+00_wp, 1.5600e+00_wp, &
      1.3600e+00_wp, 1.1600e+00_wp, 9.6000e-01_wp, 7.6000e-01_wp, &
      5.6000e-01_wp, 3.6000e-01_wp, 1.6000e-01_wp, -4.0000e-02_wp, &
      -2.4000e-01_wp, -4.4000e-01_wp, -6.4000e-01_wp, -8.4000e-01_wp, &
      -1.0400e+00_wp, -1.2400e+00_wp, -1.4400e+00_wp, -1.6400e+00_wp, &
      -1.8400e+00_wp, -2.0400e+00_wp, -2.2400e+00_wp, -2.4400e+00_wp, &
      -2.6400e+00_wp, -2.8400e+00_wp, -3.0400e+00_wp, -3.2400e+00_wp, &
      -3.4400e+00_wp, -3.6400e+00_wp, -3.8400e+00_wp, -4.0400e+00_wp, &
      -4.2400e+00_wp, -4.4400e+00_wp, -4.6400e+00_wp /), &
    !
    ! These are the temperatures associated with the respective pressures
    !
    tref(59) = (/ &
      2.9420e+02_wp, 2.8799e+02_wp, 2.7894e+02_wp, 2.6925e+02_wp, &
      2.5983e+02_wp, 2.5017e+02_wp, 2.4077e+02_wp, 2.3179e+02_wp, &
      2.2306e+02_wp, 2.1578e+02_wp, 2.1570e+02_wp, 2.1570e+02_wp, &
      2.1570e+02_wp, 2.1706e+02_wp, 2.1858e+02_wp, 2.2018e+02_wp, &
      2.2174e+02_wp, 2.2328e+02_wp, 2.2479e+02_wp, 2.2655e+02_wp, &
      2.2834e+02_wp, 2.3113e+02_wp, 2.3401e+02_wp, 2.3703e+02_wp, &
      2.4022e+02_wp, 2.4371e+02_wp, 2.4726e+02_wp, 2.5085e+02_wp, &
      2.5457e+02_wp, 2.5832e+02_wp, 2.6216e+02_wp, 2.6606e+02_wp, &
      2.6999e+02_wp, 2.7340e+02_wp, 2.7536e+02_wp, 2.7568e+02_wp, &
      2.7372e+02_wp, 2.7163e+02_wp, 2.6955e+02_wp, 2.6593e+02_wp, &
      2.6211e+02_wp, 2.5828e+02_wp, 2.5360e+02_wp, 2.4854e+02_wp, &
      2.4348e+02_wp, 2.3809e+02_wp, 2.3206e+02_wp, 2.2603e+02_wp, &
      2.2000e+02_wp, 2.1435e+02_wp, 2.0887e+02_wp, 2.0340e+02_wp, &
      1.9792e+02_wp, 1.9290e+02_wp, 1.8809e+02_wp, 1.8329e+02_wp, &
      1.7849e+02_wp, 1.7394e+02_wp, 1.7212e+02_wp /) 

  ! spec_sampling config
  integer:: rad_perm = 0 ! Integer for perturbing random number seeds

  real(wp) :: pi = 3.14159265358979323846264338327950288_wp, &
    rhoh2o = 1000._wp,  & ! density of liquid water [kg/m3]
    grav = 9.80665_wp,  & ! average gravity [m/s2]
    amd = 28.970_wp, & ! molar weight of dry air [g/mol]
    amw = 18.0154_wp, & ! molar weight of water [g/mol]
    avo = 6.02214179e23_wp ! Avogadro constant [1/mo]
  
  interface 
    subroutine t_finish_cb(name, text, exit_no)
      character(len=*), intent(in) :: name
      character(len=*), intent(in), optional :: text
      integer, intent(in), optional :: exit_no
    end subroutine
  
    subroutine t_message_cb(name, text)
      character(len=*), intent(in) :: name, text
    end subroutine
  
  end interface

  public :: message, finish, warning, default_finish 

  procedure(t_finish_cb), pointer :: finish_cb
  procedure(t_message_cb), pointer :: warning_cb => null(), &
    message_cb => null()


contains

  subroutine default_finish(name, text, exit_no)
    character(len=*), intent(in) :: name
    character(len=*), optional, intent(in) :: text
    integer, optional, intent(in) :: exit_no
    
    write(*,*) name, ": (", exit_no, ") ", text
    stop 
  end subroutine

  subroutine finish(name, text, exit_no)
    character(len=*), intent(in) :: name
    character(len=*), intent(in), optional :: text
    integer, intent(in), optional :: exit_no
    if (associated(finish_cb)) then
      call finish_cb(name, text, exit_no)
    endif
  end subroutine
  subroutine warning(name, text)
    character(len=*), intent(in) :: name, text
    if (associated(warning_cb)) then
      call warning_cb(name, text)
    endif
  end subroutine
  subroutine message(name, text)
    character(len=*), intent(in) :: name, text
    if (associated(message_cb)) then
      call message_cb(name, text)
    endif
  end subroutine

end module mo_psrad_general
