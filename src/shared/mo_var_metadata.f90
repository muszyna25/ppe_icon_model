MODULE mo_var_metadata

  USE mo_kind,          ONLY: dp
  USE mo_grib1,         ONLY: t_grib1_var
  USE mo_grib2,         ONLY: t_grib2_var
  USE mo_cf_convention, ONLY: t_cf_var

  IMPLICIT NONE

  PRIVATE

  TYPE t_union_vals
    REAL(dp) :: rval
    INTEGER  :: ival
    LOGICAL  :: lval
  END type t_union_vals

  TYPE t_var_metadata
    !
    INTEGER            :: key                   ! hash value of name
    CHARACTER(len=31)  :: name                  ! variable name  
    !
    TYPE(t_cf_var)     :: cf                    ! CF convention information 
    TYPE(t_grib1_var)  :: grib1                 ! GRIB1 related information
    TYPE(t_grib2_var)  :: grib2                 ! GRIB2 related information
    !
    LOGICAL            :: allocated             ! allocation status
    INTEGER            :: ndims                 ! number of dimensions used
    INTEGER            :: used_dimensions(5)    ! final dimensions of variable
    ! 
    LOGICAL            :: lrestart              ! write field to restart
    LOGICAL            :: lpost                 ! write field to output
    LOGICAL            :: laccu                 ! accumulation flag
    TYPE(t_union_vals) :: resetval              ! reset value for accumulated fields
    LOGICAL            :: lmiss                 ! missing value flag
    TYPE(t_union_vals) :: missval               ! missing value
    LOGICAL            :: lrestart_cont         ! continue if not in restart file     
    LOGICAL            :: lrestart_read         ! field has been set from restart file
    TYPE(t_union_vals) :: initval               ! value if not in restart file
    !     
    INTEGER            :: hgrid                 ! CDI horizontal grid type
    INTEGER            :: vgrid                 ! CDI vertical grid type
    !
  END TYPE t_var_metadata

  PUBLIC :: t_union_vals
  PUBLIC :: t_var_metadata

END MODULE mo_var_metadata
