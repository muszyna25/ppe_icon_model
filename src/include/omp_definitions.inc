!-----------------------------------------------
! define default omp parameters
#ifdef __SX__
#define ICON_OMP_DEFAULT_SCHEDULE
#define ICON_OMP_RUNTIME_SCHEDULE SCHEDULE(runtime)
#define ICON_OMP_GUIDED_SCHEDULE  SCHEDULE(guided)
#else
#define ICON_OMP_DEFAULT_SCHEDULE SCHEDULE(dynamic,1)
#define ICON_OMP_RUNTIME_SCHEDULE SCHEDULE(dynamic,1)
#define ICON_OMP_GUIDED_SCHEDULE  SCHEDULE(dynamic,1)
#endif
!-----------------------------------------------

