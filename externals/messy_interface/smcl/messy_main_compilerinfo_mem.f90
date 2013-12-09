!****************************************************************************
!                Time-stamp: <2011-02-17 13:00:30 joec_pa>
!****************************************************************************

MODULE messy_main_compilerinfo_mem

  IMPLICIT NONE
  PUBLIC
  SAVE

  ! THIS WILL BE REPLACED BY configure and/or gmake
  CHARACTER(LEN=*), PARAMETER :: compiler_version  = '@F90VERS@'
  CHARACTER(LEN=*), PARAMETER :: compiler_call     = '@F90@'
  CHARACTER(LEN=*), PARAMETER :: compiler_flags    = '@F90FLAGS@'
  CHARACTER(LEN=*), PARAMETER :: compiler_cppdefs  = '@F90DEFS@'
  CHARACTER(LEN=*), PARAMETER :: compiler_includes = '@INCLUDES@'
  
END MODULE messy_main_compilerinfo_mem
