MODULE mo_cf_convention
  !!
  !!----------------------------------------------------------------------------
  !!
  !! File name space:
  !!
  !! title
  !!
  !!    A succint description of what is in the dataset. 
  !!
  !! institution
  !!
  !!    Specifies where the original data was produced. 
  !!
  !! source
  !!
  !!    The method of production of the original data. If it 
  !!    was model-generated, source should name the model and 
  !!    its version, as specifically as could be useful.
  !!
  !! history
  !!
  !!    Provides an audit trail for modifications to the original 
  !!    data. Well-behaved generic netCDF filters will automatically 
  !!    append their name and the parameters with which they were 
  !!    invoked to the global history attribute of an input netCDF file. 
  !!
  !! references
  !!
  !!    Published or web-based references that describe the data or methods 
  !!    used to produce it. 
  !!
  !! comment
  !!
  !!    Miscellaneous information about the data or methods used to produce it. 
  !!
  !!
  !!----------------------------------------------------------------------------
  !!
  !! Variable name space:
  !!
  !! standard name
  !!
  !!    The name used to identify the physical quantity. A standard name 
  !!    contains no whitespace and is case sensitive. 
  !!
  !! canonical units
  !!
  !!    Representative units of the physical quantity. Unless it is 
  !!    dimensionless, a variable with a standard_name attribute must have 
  !!    units which are physically equivalent (not necessarily identical) to 
  !!    the canonical units, possibly modified by an operation specified by 
  !!    either the standard name modifier.
  !!
  !!  description
  !!
  !!    The description is meant to clarify the qualifiers of the fundamental 
  !!    quantities such as which surface a quantity is defined on or what the 
  !!    flux sign conventions are. We don"t attempt to provide precise 
  !!    definitions of fundumental physical quantities (e.g., temperature) 
  !!    which may be found in the literature. 
  !!
  !!----------------------------------------------------------------------------

  IMPLICIT NONE 

  PRIVATE

  TYPE t_cf_global
    CHARACTER(len=128) :: title         = ''
    CHARACTER(len=128) :: institution   = ''
    CHARACTER(len=128) :: source        = ''
    CHARACTER(len=128) :: history       = ''
    CHARACTER(len=128) :: references    = ''
    CHARACTER(len=128) :: comment       = ''
  END TYPE t_cf_global

  TYPE t_cf_var
    CHARACTER(len=128) :: long_name     = ''
    CHARACTER(len=128) :: units         = ''
    CHARACTER(len=128) :: standard_name = ''
  END TYPE t_cf_var

  PUBLIC :: t_cf_global
  PUBLIC :: t_cf_var

  PUBLIC :: set_cf_global
  PUBLIC :: set_cf_var

CONTAINS

  FUNCTION set_cf_global(title, institution, source, history, references, comment) RESULT(cf_global_info)
    TYPE(t_cf_global) :: cf_global_info
    CHARACTER(len=*), OPTIONAL :: title
    CHARACTER(len=*), OPTIONAL :: institution
    CHARACTER(len=*), OPTIONAL :: source
    CHARACTER(len=*), OPTIONAL :: history
    CHARACTER(len=*), OPTIONAL :: references
    CHARACTER(len=*), OPTIONAL :: comment
   
    IF (PRESENT(title))       cf_global_info%title       = title
    IF (PRESENT(institution)) cf_global_info%institution = institution
    IF (PRESENT(source))      cf_global_info%source      = source 
    IF (PRESENT(history))     cf_global_info%history     = history  
    IF (PRESENT(references))  cf_global_info%references  = references
    IF (PRESENT(comment))     cf_global_info%comment     = comment

  END FUNCTION set_cf_global

  FUNCTION set_cf_var(long_name, units, standard_name) RESULT(cf_var_info)
    TYPE(t_cf_var) :: cf_var_info
    CHARACTER(len=*), OPTIONAL :: long_name
    CHARACTER(len=*), OPTIONAL :: units
    CHARACTER(len=*), OPTIONAL :: standard_name

    IF (PRESENT(long_name))     cf_var_info%long_name = long_name
    IF (PRESENT(units))         cf_var_info%units     = units
    IF (PRESENT(standard_name)) cf_var_info%standard_name = standard_name
    
  END FUNCTION set_cf_var

END MODULE mo_cf_convention
