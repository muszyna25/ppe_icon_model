MODULE mo_ser_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist

  IMPLICIT NONE
  PUBLIC

  INTEGER ::  ser_output_diag
  INTEGER ::  ser_latbc_data
  INTEGER ::  ser_dynamics
  INTEGER ::  ser_diffusion
  INTEGER ::  ser_step_advection
  INTEGER ::  ser_physics
  INTEGER ::  ser_lhn
  INTEGER ::  ser_nudging
  INTEGER ::  ser_surface
  INTEGER ::  ser_microphysics
  INTEGER ::  ser_convection
  INTEGER ::  ser_cover
  INTEGER ::  ser_radiation
  INTEGER ::  ser_radheat
  INTEGER ::  ser_gwdrag
  INTEGER ::  ser_all_debug                      !! serialize statements using ser_all anywhere for debug purposes
  LOGICAL ::  ser_debug                          !! serialize the debug calls from mo_ser_debug

  NAMELIST /ser_nml/ ser_output_diag, ser_latbc_data, ser_dynamics, ser_debug, &
  &                  ser_diffusion, ser_step_advection, ser_physics, ser_lhn, ser_nudging, ser_all_debug, ser_surface, &
  &                  ser_microphysics, ser_convection, ser_cover, ser_radiation, ser_radheat, &
  &                  ser_gwdrag

  CONTAINS

  SUBROUTINE read_ser_namelist( filename )

   CHARACTER(LEN=*), INTENT(IN) :: filename
   INTEGER :: istat, funit 
   INTEGER :: iunit

   ! turn serialization off by default
   ser_output_diag = 0
   ser_latbc_data = 0
   ser_dynamics = 0
   ser_diffusion = 0
   ser_step_advection = 0
   ser_physics = 0
   ser_lhn = 0
   ser_nudging = 0
   ser_surface = 0
   ser_microphysics = 0
   ser_convection = 0
   ser_cover = 0
   ser_radiation = 0
   ser_radheat = 0
   ser_gwdrag = 0
   ser_all_debug = 0
   ser_debug = .FALSE.

   IF (my_process_is_stdio()) THEN
     iunit = temp_defaults()
     WRITE(iunit, ser_nml)   ! write defaults to temporary text file
   END IF

   CALL open_nml(TRIM(filename))
   CALL position_nml ('ser_nml', status=istat)

   SELECT CASE (istat)
   CASE (POSITIONED)
      READ (nnml, ser_nml)   ! overwrite default settings

      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, ser_nml)   ! write settings to temporary text file
      END IF
   END SELECT

   IF(my_process_is_stdio())  THEN
     funit = open_tmpfile()
     WRITE(funit,NML=ser_nml)                    
     CALL store_and_close_namelist(funit, 'ser_nml') 
   ENDIF

   IF(my_process_is_stdio()) WRITE(nnml_output,nml=ser_nml)

  END SUBROUTINE read_ser_namelist

END MODULE mo_ser_nml
