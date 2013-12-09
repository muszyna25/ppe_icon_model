! **********************************************************************
!
! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL BUFLY
!
! THIS SUBMODEL IS USED TO ILLUSTRATE THE BUTTERFLY EFFECT
! VIA INITIAL SURFACE PRESSURE AND/OR TEMPERATURE PERTURBATIONS.
!
! MOREOVER, THIS SUBMODEL SERVES AS A TEMPLATE FOR NEW SUBMODELS.
!
! Author : Patrick Joeckel, DLR-IPA, October  2009
!
! References:
!
! * P. Jöckel, R. Sander, A. Kerkweg, H. Tost, and J. Lelieveld,
!   Technical Note: The Modular Earth Submodel System (MESSy) - a new
!   approach towards Earth System Modeling,
!   Atmos. Chem. Phys., 5, 433-444, 2005.
!   http://www.atmos-chem-phys.net/5/433 
! * Patrick Jöckel, Technical note: Recursive rediscretisation of
!   geo-scientific data in the Modular Earth Submodel System (MESSy),
!   Atmos. Chem. Phys., 6, 3557-3562, 2006.
!   http://www.atmos-chem-phys.net/6/3557
! * P. Jöckel, A. Kerkweg, J. Buchholz-Dietsch, H. Tost, R. Sander, and
!   A. Pozzer, Technical Note: Coupling of chemical processes with the
!   Modular Earth Submodel System (MESSy) submodel TRACER,
!   Atmos. Chem. Phys., 8, 1677-1687, 2008.
!   http://www.atmos-chem-phys.net/8/1677 
!
! **********************************************************************

! **********************************************************************
MODULE messy_bufly
! **********************************************************************

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'bufly'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.1'
 
  ! CTRL-NAMELIST PARAMETERS
  ! latitude +/- width
  REAL(DP), DIMENSION(2), PUBLIC :: r_lat = (/0.0_dp, 0.0_dp/)
  ! langitude +/- width
  REAL(DP), DIMENSION(2), PUBLIC :: r_lon = (/0.0_dp, 0.0_dp/)
  ! perturbation temperature [K]
  REAL(DP), PUBLIC               :: r_dt = 0.0_dp
  ! perturbation surface pressure [Pa]
  REAL(DP), PUBLIC               :: r_dps = 0.0_dp

  ! PUBLIC SUBROUTINES (to be called from messy_bufly_e5.f90)
  PUBLIC :: bufly_read_nml_ctrl
  ! ### add your own public subroutines here
  PUBLIC :: perturb

  ! PRIVATE SUBROUTINES
  ! ### add your own private subroutines here

CONTAINS

  ! =========================================================================
  ! ### add your own public subroutines here
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE perturb(temp, glon, glat)

    IMPLICIT NONE
    INTRINSIC :: SIZE, ABS

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(OUT) :: temp ! perturbation temperature
    REAL(DP), DIMENSION(:), INTENT(IN)  :: glon ! grid longitude
    REAL(DP), DIMENSION(:), INTENT(IN)  :: glat ! grid latitude
    
    ! LOCAL
    INTEGER :: jp
    REAL(DP) :: zlon

    ! INIT
    temp(:) = 0.0_dp

    ! allow longitude in [-180,...,360]
    IF (r_lon(1) < 0.0_dp) THEN
       zlon = r_lon(1) + 360.0_dp
    ELSE
       zlon = r_lon(1)
    END IF

    DO jp=1, SIZE(temp)
       IF ( (ABS(glon(jp) - zlon)     <= r_lon(2)) .AND. &
            (ABS(glat(jp) - r_lat(1)) <= r_lat(2)) ) THEN
          temp(jp) = r_dt
       ENDIF
    END DO

  END SUBROUTINE perturb
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE bufly_read_nml_ctrl(status, iou)

    ! ------------------------------------------------------------------
    ! This routine is used to read the CTRL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy INTERFACE
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/ r_lon, r_lat, r_dt, r_dps

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='bufly_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR
    
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    
    ! ### ADD HERE DIAGNOSTIC OUPUT FOR LOG-FILE
    WRITE(*,*) 'INITIAL PERTURBATION OF ',r_dt,' K AT ' &
         , r_lat(1),'+/-',r_lat(2),' (latitude) / '     &
         , r_lon(1),'+/-',r_lon(2),' (longitude)'
    WRITE(*,*) 'INITIAL PERTURBATION OF ',r_dps,' Pa AT ' &
         , r_lat(1),' (latitude) / '         &
         , r_lon(1),' (longitude)'

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR
    
  END SUBROUTINE bufly_read_nml_ctrl
  ! =========================================================================

  ! =========================================================================
  ! ### add your own private subroutines here
  ! =========================================================================

! **********************************************************************
END MODULE messy_bufly
! **********************************************************************

