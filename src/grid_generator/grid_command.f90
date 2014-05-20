PROGRAM grid_command

  USE mo_kind
  USE mo_io_units,  ONLY: filename_max
  USE mo_exception, ONLY: finish
  
  USE mo_global_grid_generator, ONLY: global_graph_generator, global_grid_generator, &
    & prepare_gridref

  USE mo_run_config,            ONLY: profiling_output, TIMER_MODE_DETAILED
 
  USE mo_mpi,                   ONLY: start_mpi, stop_mpi

  IMPLICIT NONE

  CHARACTER(len=filename_max) :: command_file = 'command.grid'
  
  ! drives the old grid generator
  CHARACTER(len=32), PARAMETER :: global_graph_generator_c ='global_graph_generator'
  CHARACTER(len=32), PARAMETER :: global_grid_generator_c ='global_grid_generator'
  CHARACTER(len=32), PARAMETER :: global_grid_refine_c ='global_grid_refine'

  CHARACTER(len=32) :: command =''
  CHARACTER(len=filename_max) :: param_1 = ''
  CHARACTER(len=filename_max) :: param_2 = ''
  CHARACTER(len=filename_max) :: param_3 = ''

  INTEGER :: int_param_1, int_param_2, int_param_3
  REAL(wp):: real_wp_param_1, real_wp_param_2, real_wp_param_3, real_wp_param_4, &
    & real_wp_param_5, real_wp_param_6
  REAL :: real_param_1
  LOGICAL :: logical_flag
!   CALL get_command_argument(1, command)
!   CALL get_command_argument(2, param_1)

  CALL start_mpi()

  profiling_output = TIMER_MODE_DETAILED

  OPEN (500, FILE = command_file,STATUS = 'OLD')
  READ (500, *) command, param_1
  CLOSE(500)
  print *, TRIM(command), ', ', TRIM(param_1)

  SELECT CASE (command)

    CASE (global_graph_generator_c)
      CALL global_graph_generator()

    CASE (global_grid_generator_c)
      CALL global_grid_generator()
      
    CASE (global_grid_refine_c)
      CALL prepare_gridref()

    CASE DEFAULT
      CALL finish('grid_command', 'Unkown command')

  END SELECT
  !----------------------------------------------------------
 
  CALL stop_mpi()

END PROGRAM grid_command
