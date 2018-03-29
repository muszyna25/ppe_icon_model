MODULE mo_psrad_general
  use, intrinsic :: iso_c_binding, only: c_float, c_double, c_long, c_int

  IMPLICIT NONE
  PUBLIC

  LOGICAL :: upwards
  INTEGER :: jTOA, jSFC, jINC, jABOVE, jBELOW
  INTEGER, PARAMETER :: &
    dp = c_double, wp = c_double, &
    i8 = c_long, i4 = c_int, &
    ngpt_orig = 16, &!< number of original g-intervals per spectral band
    jpband = 29, & !< number of last band (lw and sw share band 16)
    nbndsw = 14, & !< number of spectral bands in sw model
    nbndlw = 16, & !< number of spectral bands in lw model
    npressure(2) = (/13,47/), &
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
    nmixture = 6, ih2oco2 = 1, ih2oo3 = 2, ih2on2o = 3, ih2och4 = 4, &
    in2oco2 = 5, io3co2 = 6

  REAL(wp), PARAMETER :: ten20 = 1e20_wp, ten20inv = 1e-20_wp
  REAL(wp), ALLOCATABLE :: dummy4(:,:,:,:)

  REAL(wp), PARAMETER :: oneminus = 1.0_wp - 1.0e-06_wp

  ! spec_sampling config
  INTEGER:: rad_perm = 0 ! Integer for perturbing random number seeds

  REAL(wp) :: pi = 3.14159265358979323846264338327950288_wp, &
    rhoh2o = 1000._wp,  & ! density of liquid water [kg/m3]
    grav = 9.80665_wp,  & ! average gravity [m/s2]
    amd = 28.970_wp, & ! molar weight of dry air [g/mol]
    amw = 18.0154_wp, & ! molar weight of water [g/mol]
    avo = 6.02214179e23_wp ! Avogadro constant [1/mo]
  
  INTERFACE 
    SUBROUTINE t_finish_cb(name, text, exit_no)
      CHARACTER(len=*), INTENT(IN) :: name
      CHARACTER(len=*), INTENT(IN), optional :: text
      INTEGER, INTENT(IN), OPTIONAL :: exit_no
    END SUBROUTINE
  
    SUBROUTINE t_message_cb(name, text)
      CHARACTER(len=*), INTENT(IN) :: name, text
    END SUBROUTINE
  
  END INTERFACE

  PUBLIC :: message, finish, warning, default_finish 

  PROCEDURE(t_finish_cb), POINTER :: finish_cb
  PROCEDURE(t_message_cb), POINTER :: warning_cb => null(), &
    message_cb => null()


CONTAINS

  SUBROUTINE default_finish(name, text, exit_no)
    CHARACTER(len=*), INTENT(IN) :: name
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: text
    INTEGER, OPTIONAL, INTENT(IN) :: exit_no
    
    write(*,*) name, ": (", exit_no, ") ", text
    STOP 
  END SUBROUTINE

  SUBROUTINE finish(name, text, exit_no)
    CHARACTER(len=*), INTENT(IN) :: name
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: text
    INTEGER, INTENT(IN), OPTIONAL :: exit_no
    IF (associated(finish_cb)) THEN
      CALL finish_cb(name, text, exit_no)
    ENDIF
  END SUBROUTINE
  SUBROUTINE warning(name, text)
    CHARACTER(len=*), INTENT(IN) :: name, text
    IF (associated(warning_cb)) THEN
      CALL warning_cb(name, text)
    ENDIF
  END SUBROUTINE
  SUBROUTINE message(name, text)
    CHARACTER(len=*), INTENT(IN) :: name, text
    IF (associated(message_cb)) THEN
      CALL message_cb(name, text)
    ENDIF
  END SUBROUTINE

END MODULE mo_psrad_general
