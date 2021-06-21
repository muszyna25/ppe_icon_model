!

MODULE mo_bench_config

  !USE...

  IMPLICIT NONE

  PRIVATE
  
  ! types
  PUBLIC :: t_bench_config
  ! variables
  PUBLIC :: bench_config

  TYPE t_bench_config
  LOGICAL :: d_unpb, d_ndfo, d_rld, d_n, d_wnlo
  END TYPE t_bench_config

  TYPE(t_bench_config), TARGET :: bench_config
  
  !CONTAINS

END MODULE mo_bench_config