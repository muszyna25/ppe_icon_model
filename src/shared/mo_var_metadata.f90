MODULE mo_var_metadata

  USE mo_kind, ONLY: dp
  USE mo_grib1, ONLY: grib1_var
  USE mo_grib2, ONLY: grib2_var
  USE mo_cf_convention, ONLY: cf_var

  IMPLICIT NONE

  TYPE var_metadata
    !
    INTEGER           :: key                   ! hash value of name
    CHARACTER(len=31) :: name                  ! variable name  
    !
    TYPE(cf_var)      :: cf                    ! CF convention information 
    TYPE(grib1_var)   :: grib1                 ! GRIB1 related information
    TYPE(grib2_var)   :: grib2                 ! GRIB2 related information
    !
    LOGICAL           :: allocated             ! allocation status
    INTEGER           :: allocated_dimensions  ! number of dimensions used
    INTEGER           :: used_dimensions(4)    ! final dimensions of variable
    ! 
    LOGICAL           :: lrestart              ! write field to restart
    LOGICAL           :: lpost                 ! write field to output
    LOGICAL           :: laccu                 ! accumulation flag
    REAL(dp)          :: reset                 ! reset value for accumulated fields
    LOGICAL           :: lmiss                 ! missing value flag
    REAL(dp)          :: missval               ! missing value
    LOGICAL           :: lrestart_cont         ! continue if not in restart file     
    LOGICAL           :: lrestart_read         ! field has been set from restart file
    INTEGER           :: grid_representation   ! CDI grid representation type
    !
    ! CDI handler
    !
    INTEGER           :: gridID
    INTEGER           :: zaxisID 
    !
  END type var_metadata

END MODULE mo_var_metadata
